#Created by Ziqing Ye, March 08, 2021
#Last modified: Nov.4, 2021
##MODIFIED TO KEEP AT-rich sequences < GC threshold
import sys

#functions for palindrome.py for Siran Tian in Trcek Lab at JHU
#cleaned up version, removed inefficient/repetitive functions


#use a dict to save position(int) and palindrome(string) pairs
#Params: string sequence (the DNA sequence), int min_len=minimum palindrome length,
#max_len = maximum palindrome length (in nucleotides)
#gc_thresh: float out of 1, percent minimum gc content
#mis_thresh: tolerance (int in nucleotides) for mismatch in a palindromic sequence
#returns: a dict of position-palindrome pairs
#uses the center of palindrome as the key; returns an int containing total palindrome count (regardless of gc) and a dict of desired gc-rich palindromes
def find_palindrome(sequence, min_len, max_len, gc_thresh, mis_thresh):
    counter = 0 #total palindromes (with distinct centers) regardless of GC content
    seq_len=len(sequence)
    if (max_len > seq_len):
        max_len = seq_len
    results = {}  #initializing an empty dict
    all_palin_centers = []
    for i in range(0, seq_len-min_len+1): #i is the starting position of the palindrome
        j = min_len #j is the length of each palindrome
        while(j<max_len+1 and (i+j)<=seq_len): #start from searching palindrome of len min, increment till searching for seq of len max
            sub_seq=sequence[i:(i+j)] #returns a substring starting from the ith char that is j characters long
            if (j%2 != 0): #skip palindromes of odd length since these will be count as inverted repeats
                j += 0 #do nothing
            elif (DNA_is_imperfect_palindrome(sub_seq, mis_thresh)==True): #check if is palindromic sequence
                gc = perc_gc_content(sub_seq)#calculate gc content
                center = (i+i+j-1)/2
                if (center not in all_palin_centers):
                    counter += 1 #increment total distinct palindrome count if it is a new palindrome w/ a diff center
                    all_palin_centers.append(center);
                #print(i, i+j-1, sub_seq, center)
                #at this point, we could either discard those that do not pass the GC thresh, or keep them
                #if we wanna keep all regardless of gc content, remove the gc part of this if statement
                if (gc < gc_thresh and (not (center in results))):
                    results[(i+i+j-1)/2]=(i, sub_seq, i+j-1) #store center as the key, and starting-index+substring+ending-index as value (Python dict)
            j+=1 #increment palindrome length
        #end of while loop
    #end of for loop
    #print(results)
    return counter, results #return the results dict for output



#calculates the percent gc content of a given sequence
#param seq: a nucleotide sequence in all caps
#returns a float of percent gc content
def perc_gc_content(seq):
    g_count = seq.count('G')
    c_count = seq.count('C')
    return (g_count+c_count)/len(seq)


#Reverse complement string method
#returns the complement nucleotide of the input nucleotide
def comp(nucleotide):
    if (nucleotide=='A'):
        return 'T'
    if (nucleotide=='T'):
        return 'A'
    if (nucleotide=='C'):
        return 'G'
    if (nucleotide=='G'):
        return 'C'


##checks whether a dna sequence is palindromic enough based on
# a given maximum mismatch tolerance, in nucleotides, without
# involving any external library (use only string methods)
#is a different approach I came up with that is more efficient
#than brute force comparing a string with its generated reverse complement
#this should be about twice as efficient cuz only looks at half of the seq
#(compares 1st half with second half)
# param: string sequence, int thresh=mismatch max tolerance
#returns true if the sequence can be considered palindromic, false otherwise
def DNA_is_imperfect_palindrome(seq, tol):
    lens = len(seq)
    half_len = round(lens/2) #will round up to include middle nucleotide if odd length
    max_ind=lens-1
    mis = 0
    for i in range (half_len):
        if seq[i] != comp(seq[max_ind-i]):
            mis+=1
    return (mis<=round(tol/2)) #true if mismatch less than threshold, false otherwise
    
