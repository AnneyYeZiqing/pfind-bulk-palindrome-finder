#Palindrome summary counts to excel
#Ziqing Ye, Siran Tian
#THIS VERSION IS MODIFIED TO IDENTIFY AT RICH SEQUENCES INSTEAD (E.G. PALINDROMES WITH GC CONTENT LESS THAN (<) the input GC rich threshold)
#updated 9/13: added support for using an absolute address for the input fasta (added line 27 and changed line 28, 29)
#updated 9/13: added cds header format support
#updated 9/24: finished implementing file output format selection command line arguments
#updated 11/19: add path support for both unix and windows
#updat3e 12/13/2021: debugged path support for unix/mac
"""
usage:.py [1]fasta [2]min_length [3]max_length [4]GCthreshold(float out of 1) 
[5]mismatch_tolerance(int, # base pairs on HALF of the palindrome) [6] need_sum (0=no output csv count summary, 1=yes, default yes)
[7]need_palindrome_specifics (txt output with each palindrome's sequence and location data) (0=no, 1=expSSformat, 2=vertical-per-palindrome, default 1)
Sum output format (excel columns):Gene, transcripts, region (3'utr vs cDNA), length (bases) , total palindromes, palindromes >=GC threshold

"""

import sys
import time
import collections
import os
#below are custom modules
from fasta import FASTAReader
import pfind_ATrich #make sure the custom python script pfind and fasta is in the same folder/directory as the current script
from pathlib import Path, PurePath

def main():
    #for command-line args
    #running format: paindrome.py filename.fa min_len max_len gc_thresh mismatch_thresh (whether_count_summary) (palindrome_output)
    reader1= FASTAReader(open(sys.argv[1]))
    min_length=int(sys.argv[2])   #minimum length (nucleotides)
    max_length=int(sys.argv[3]) #maximum length of palindromes 
    gc_thresh = float(sys.argv[4]) #percent gc content
    mismatch_thresh = int(sys.argv[5]) #mismatch tolerance base pairs
    count_summary = 1
    palind_specifics = 1
    if len(sys.argv) >= 7:
        count_summary = int(sys.argv[6])
    if len(sys.argv) >= 8:
        palind_specifics = int(sys.argv[7])
    
    fasta_name = Path(sys.argv[1]).parts[-1] #get the file name from a path
    
    if count_summary != 0:
        os.makedirs(PurePath('csvcounts/'), exist_ok=True)
        filecsv = open(PurePath("csvcounts/counts-"+fasta_name+"-"+sys.argv[2]+"-"+sys.argv[3]+"-"+sys.argv[4]+"-"+sys.argv[5]+"ATrich.csv"), 'w')
        filecsv.write("#param: min-%d, max-%d, gc-%f, mismatch-%d\n" %(min_length, max_length, gc_thresh, mismatch_thresh))
        filecsv.write("Gene, Transcripts, Region, Length, total palindromes, palindromes>=%%GC%.1f \n" %(gc_thresh*100)) #first percent sign to escape and print the second percent sign
    if palind_specifics == 1: #expSSformat, all transcripts in one file
        os.makedirs(PurePath('output-forExpSS/'), exist_ok=True)
        file1 = open(PurePath("output-forExpSS/palind-expSSform-"+fasta_name+"-"+sys.argv[2]+"-"+sys.argv[3]+"-"+sys.argv[4]+"-"+sys.argv[5]+"ATrich.dat"), 'w')
    elif palind_specifics == 2:
        os.makedirs(PurePath('per-palindrome/'), exist_ok=True)
        file2 = open(PurePath("per-palindrome/"+fasta_name+"-"+sys.argv[2]+"-"+sys.argv[3]+"-"+sys.argv[4]+"-"+sys.argv[5]+"ATrich.csv"), 'w')
        
    before_time = time.perf_counter() #used to calculate runtime
    for ident, sequence in reader1:
        fields = ident.rstrip("\n").split()
        #[0]:FBtr0084360 [1]:type=three_prime_untranslated_region; [2]:loc=3R:complement(23049657..23049964);MD5=2315b6abd3f58fa146e32dd60dcb2adc;
        #[3]:length=308; [4]:parent=FBgn0039044; [5]:release=r6.33; [6]:species=Dmel;
        #>FBtr0344352 type=three_prime_UTR; loc=3R:19913085..19913730;MD5=eeed0bffb52c4a9dc45a338e1ff5da8d; length=646; parent=FBgn0038747; release=r6.41; species=Dmel; 
        #for cds: (split by space) 
        #[0]FBpp0088566 [1]eIF4G2-PA [2]type=CDS; [3]loc=3R:complement(join(23810294..23812585, 1..3)); [4]name=eIF4G2-RA;
        #[5]dbxref=FlyBase:FBpp0088566,REFSEQ:NP_651188; [6]MD5=1b87cb4321e6290bffcc5bcb3bb4ab8f; [7]length=6219; [8]parent=FBgn0260634,FBtr0089623;
        #[9]release=r6.41; [10]species=Dmel;
        #for other-way-downloaded cds:
        #>[0]CG9570-PA [1]type=CDS; [2]loc=X:20051456..20052088; [3]name=CG9570-RA; [4]dbxref=xxx1; [5]MD5=d2db82dff0247b192fa77af3fa16358f; [6]length=633;
        #[7]parent=FBgn0031085,FBtr0070002; [8]release=r6.33; [9]species=Dmel; 
        #part 1: get_transcript_info(fields)
        transcript_name = fields[0]
        if transcript_name[:4] == 'FBtr':
            gene_name = fields[4][7:-1] #parent gene in FBgnXXXX
            reg_name = fields[1][5:-1] #region name (e.g. 3'UTR)
            length = fields[3][7:-1] #nucleotide length
            if count_summary != 0:
                filecsv.write("%s, %s, %s, %s, " %(gene_name, transcript_name, reg_name, length))
        elif transcript_name[:4] == 'FBpp': #for cds headers
            gene_name = fields[8][7:-1] #parent=FBgn0260634,FBtr0089623;
            gene_symbol = fields[1] #eIF4G2-PA
            reg_name = fields[2][5:-1]#type=CDS
            length = fields[7][7:-1] #nucleotide length, length=6219;
            if count_summary != 0:
                filecsv.write("%s, %s, %s, %s, %s, " %(gene_name, transcript_name, gene_symbol, reg_name, length)) #add an extra column for CDS, can remove easily
        else:
            #assume it is the other cds but might not work
            gene_name = fields[7][7:-1] #parent=FBgn0260634,FBtr0089623;
            gene_symbol = fields[0] #eIF4G2-PA
            reg_name = fields[1][5:-1]#type=CDS
            length = fields[6][7:-1] #nucleotide length, length=6219;
            #print(gene_name, transcript_name, gene_symbol, reg_name, length)
            if count_summary != 0:
                filecsv.write("%s, %s, %s, %s, %s, " %(gene_name, transcript_name, gene_symbol, reg_name, length)) #add an extra column for CDS, can remove easily
        
        if palind_specifics == 1:
            file1.write(">" + transcript_name + "\n")
        result = pfind_ATrich.find_palindrome(sequence, min_length, max_length, gc_thresh, mismatch_thresh) #int, dict
        #result[0]=count of total palindromes regardless of GC threshhold
        #result[1]=the dict containing all palindromes within the GC threshold

        if count_summary != 0:
            filecsv.write("%d, %d \n" %(result[0], len(result[1])))#total palindromes regardless of GC, palindromes within GC thresh
        if palind_specifics == 1: #expSS-compatible format
            #output format: 
            #>FBtr:FBgn:region:length
            #startInd:ATCGAT:endInd, startInd:AATATATT:endInd
            for value in result[1].values(): #iterate over each palindrome in a dict
                #each palindrome value is a list of (startInd, string sequence, endInd)
                #the ind is the string index, which is 1 less than the nucleotide position
                #we want the actual nucleotide position, ind + 1.
                file1.write("%d:%s:%d," %(value[0]+1, value[1], value[2]+1))
                #e.g. 35:ACCGGT:40 means A is 35, C is 36, ... T is nucleotide position 40
            file1.write("\n")
        if palind_specifics == 2: #vertical rows of all palindromes's sequences in all transcripts in the fasta file, uninterrupted
            for value in result[1].values(): #iterate over each palindrome in a dict
                #each palindrome value is a list of (startInd, string sequence, endInd)
                #the ind is the string index, which is 1 less than the nucleotide position
                #we want the actual nucleotide position, ind + 1.
                #format: starting-nt, sequence, source-transcript, source-gene (the last two are for dedup)
                file2.write("%d,%s,%s,%s\n" %(value[0]+1, value[1], transcript_name, gene_name))

    after_time = time.perf_counter()
    print("Runtime is " + str(after_time-before_time) + " seconds")
    if count_summary != 0:
        filecsv.close()
    if palind_specifics == 1: #expSSformat, all transcripts in one file
        file1.close()
    elif palind_specifics == 2:
        file2.close()

main()
