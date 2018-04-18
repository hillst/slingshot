import sys
import pysam
import os
import cigar as cr
import argparse

def main():
    ## General Tool Information
    desc= "Perform Duplication Consensus Test"
    long_desc="tbd"
    author="Sunil Deochand"
    lastupdated="2018-04-17"

    parser = argparse.ArgumentParser(description='{}'.format(desc), epilog=long_desc)
    parser.add_argument('--bam', required=True, help='Input Bam File')
    parser.add_argument('--vcf', required=True, help='List of Variant Sites (Can work with any format that has the first two columns being chr and position, 1 indexed)')
    args = parser.parse_args()
    print(args)

    ## Function to convert cigar string to pysam list convention 
    do_work(args)

    os.system("samtools sort -o " + os.path.basename( str(args.bam) ) +".sorted.dupetested.bam" + " " +  os.path.basename( str(args.bam) ) + ".dupetested.bam" )
    os.system("samtools index " + os.path.basename( str(args.bam) )  + ".sorted.dupetested.bam" )

def pysam_cigar(x):
    _convertcigar = dict( zip( ['M','I','D', 'N','S', 'H'] , range(6) ) )
    as_list = list(cr.Cigar(x).items())
    rev_list = [ tuple( (i[1],i[0]) ) for i in as_list]
    pysam_ready = [ tuple( (_convertcigar[i[0]], i[1])) for i in rev_list ]
    return(pysam_ready)

def dupeID(read):
    """
    Input assumes...
    ## Function to create a dupe family ID  (Note: U.P. = unclipped 5' position) as a new TAG to each read
    ## Defined as: If on '+' strand : "read's chr" + "read's U.P." + "Mate's chr" + "Mate's U.P."
    ##               If on '-' strand :  "Mate's chr" + "Mate's U.P." + "read's chr" + "read's U.P." 

    :read:

    :returns: ID corresponding to this duplicate family or None (if there is no mate)
    """
    #s# If read on '+' strand
    if not read.is_reverse:
        ## Find read's unclipped 5' position by using its cigar
        unclipped_start = read.reference_start
        for i in range(0,len(read.cigar)):
            if read.cigar[i][0] == 4 and i == 0: ## Subtract softclip at start of read
                unclipped_start -= read.cigar[i][1]
            if read.cigar[i][0] == 5 and i == 0: ## Subtract hardclip at start of read
                unclipped_start -= read.cigar[i][1]

        ## Find the mate's unclipped 5' position using its cigar given in the 'MC' tag
        tag_dict = dict(read.tags) #was a set previously
        if not "MC" in tag_dict: #read mate not present or not mapped -- smarter way to exit out
            return None

        matecigar = pysam_cigar(tag_dict['MC'])
        mate_unclipped_start = read.next_reference_start
        for i in range(0, len(matecigar)):
            if matecigar[i][0] == 0:  ## Add matches
                mate_unclipped_start += matecigar[i][1]
            if matecigar[i][0] == 2:  ## Add deletions
                mate_unclipped_start += matecigar[i][1]            
            if matecigar[i][0] == 4 and i==len(matecigar)-1 :  ## Add soft clip at the end of read
                mate_unclipped_start += matecigar[i][1]
            if matecigar[i][0] == 5 and i==len(matecigar)-1 :  ## Add hard clip at the end of read
                mate_unclipped_start += matecigar[i][1]        
        ID = str(read.reference_id) + str(unclipped_start) + str(read.next_reference_id) + str(mate_unclipped_start)
        return(ID)
    ## If read is on the '-' strand 
    else:
        ## Find read's unclipped 5' position by using its cigar
        unclipped_start = read.reference_start
        for i in range(0, len(read.cigar)):
            if read.cigar[i][0] == 0:  ## Add matches
                unclipped_start += read.cigar[i][1]
            if read.cigar[i][0] == 2:  ## Add deletions
                unclipped_start += read.cigar[i][1]            
            if read.cigar[i][0] == 4 and i==len(read.cigar)-1:  ## Add soft clip at the end of read
                unclipped_start += read.cigar[i][1]
            if read.cigar[i][0] == 5 and i==len(read.cigar)-1:  ## Add hard clip at the end of read
                unclipped_start += read.cigar[i][1]        

        
        ## Find the mate's unclipped 5' position using its cigar given in the 'MC' tag
        tag_dict = dict(read.tags) #was a set previously
        if not "MC" in tag_dict: #read mate not present or not mapped -- smarter way to exit out
            return None

        matecigar = pysam_cigar(tag_dict['MC'])
        mate_unclipped_start = read.next_reference_start
        for i in range(0,len(matecigar)):
            if matecigar[i][0] == 4 and i == 0: ## Subtract softclip at start of read
                mate_unclipped_start -= matecigar[i][1]
            if matecigar[i][0] == 5 and i == 0: ## Subtract hardclip at start of read
                mate_unclipped_start -= matecigar[i][1]
        ID = str(read.next_reference_id) + str(mate_unclipped_start) + str(read.reference_id) + str(unclipped_start)
        return(ID)

def do_work(args):
    """
    Main loop responsible for:
        Argument handling
        I/O (opening and closing files)
        
        
    """
    ## 1st Argument = BAM File
    try:
        BAM=pysam.AlignmentFile( str(args.bam) ,"rb") 
        consensus_bam = pysam.AlignmentFile( os.path.basename( str(args.bam) ) + ".dupetested.bam" , "wb", template=BAM)
    except Exception as e:
        print("Error opening bam: {}".format(args.bam))
        return -1

    try:
        sites = open( str(args.vcf) )
        sites_list = sites.readlines()
        sites.close()
    except Exception as e:
        print("Error opening sites list: {}".format(args.vcf))
        return -1

    ## For each site in the list, get all reads overlapping within 200 bps of that site
    for line in sites_list:
        if line[0][0]=='#':
            continue
        site = line.split()
        read_IDs = []
        read_set = []
        ## Get a List of the Duplicate Family IDs and List with reads corresponding to those IDS
        for read in BAM.fetch( site[0] , int(site[1])-1 , int(site[1]) ):
            read_id = dupeID(read)
            if read_id is not None:
                read_IDs.append(dupeID(read))
                read_set.append(read)

        ## For each Dupe Family ID, Get all reads with that ID and Output Consensus Read
        for ID in read_IDs:
            ## For this particular ID, finds index of all its occurrences in the ID list
            indices = [i for i, x in enumerate(read_IDs) if x == ID]
            ## A set of pair end reads has 2 reads (2 IDs), a duplicate family has 4 or more, >2 is a sufficient test
            if len(indices)>2:
                ## Get the reads corresponding to the ID indices
                dupefamily = [read_set[i] for i in indices]
                ## To Hold the reads that have a base overlapping the site of interest (deletion or N causes this to not be the case)
                dupefamily_overlapping_reads = []
                ## To Hold the base of each of the overlapping reads of the dupe family to determine the majority
                bases = []
                # Loop over the family, recording the base of each read at the site
                for read in dupefamily:
                    ## Proceed only if the read has a base that overlaps the site
                    if read.get_overlap( int(site[1])-1 , int(site[1])) == 1:
                        dupefamily_overlapping_reads.append(read)
                        site_index = read.get_reference_positions(full_length=True ).index( int(site[1])-1 )
                        bases.append( read.get_aligned_pairs( with_seq = True )[site_index][2] )
                ## Choose a read with the majority base in this dupe family and output that
                nucs = ['A', 'T', 'C', 'G']
                base_frequencies = [ float(bases.count(i))/float(len(bases)) for i in nucs ]
                ## Find index of the max base frequency
                maxes = [i for i, x in enumerate(base_frequencies) if x == max(base_frequencies)]
                ## If no ties in frequency, choose first read with the majority base
                if len(maxes)==1:
                    chosen_base = nucs[ base_frequencies.index(max(base_frequencies)) ]
                    chosen_read = dupefamily_overlapping_reads[ bases.index(chosen_base) ]
                    consensus_bam.write(chosen_read)
    consensus_bam.close()


if __name__ == "__main__":
    main()
