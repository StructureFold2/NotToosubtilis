#!/usr/bin/env python
#David Tack

'''
Takes RockHopper transcript predictions and applies them to a genomic fasta
therefore making putative transcripts.
Someone decided that commas in names/descriptions are a good thing, so change the <.txt>
transcripts file to semicolon demlimited, and this should get all the fields just fine.
'''

#Imports
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Classes
class RockHopperRecord(object):
    '''Holds Information'''
    def __init__(self,record_line):
        components = record_line.strip().split(';')
        self.tr_start = int(components[0]) if components[0].isdigit() else None
        self.tl_start = int(components[1]) if components[1].isdigit() else None
        self.tl_stop = int(components[2]) if components[2].isdigit() else None
        self.tr_stop = int(components[3]) if components[3].isdigit() else None
        self.strand = components[4]
        self.short_name = components[5]
        self.BSU = components[6]
        self.long_name = components[7]
        self.consistent = all([self.tr_start,self.tl_start,self.tl_stop,self.tr_stop])
        self.full_name = '-'.join([self.short_name,self.BSU])

#Functions
def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def write_fasta(fastadict,outfile,LW=80):
    '''Derp'''
    with open(outfile, 'w') as out:
        for k,v in fastadict.items():
            out.write('>' + k + '\n')
            for i in xrange(0,len(v),LW):
                out.write(v[i:i+LW] + '\n') 

def read_in_rockhopper_transcripts(t_file):
    '''Reads In RockHopper Transcripts'''
    generic_numbs,out_info = iter(range(1,100000)),{}
    with open(t_file,'r') as f:
        for line in f:
            if line.startswith('Transcription Start') or line.strip() =='':
                continue
            else:
                out_info[generic_numbs.next()] = RockHopperRecord(line)
    return out_info

def generate_transcripts(annotation,genome_sequence):
    '''Grab Sequence, Generate Name, Return Information'''
    fasta_dict = {}
    for number, entry in annotation.items():
        
        #For Entries that do have all four fields
        if entry.consistent:
            coords = [entry.tr_start,entry.tl_start,entry.tl_stop,entry.tr_stop]
            left,right = left, right = min(coords), max(coords)
            
            if entry.strand == '+':
                full_sequence = genome_sequence[left-1:right]
                l_brace = (entry.tl_start - entry.tr_start)+1
                r_brace = (entry.tl_stop - entry.tl_start)+l_brace
                numbs = '~'.join([str(l_brace),str(r_brace)])
                entry_name = '.'.join(['ILL',entry.full_name,numbs])
                fasta_dict[entry_name] = full_sequence

            elif entry.strand == '-':
                full_sequence = str(Seq(genome_sequence[left-1:right]).reverse_complement())
                l_brace = (entry.tr_start - entry.tl_start)+1
                r_brace = (entry.tl_start - entry.tl_stop)+l_brace
                numbs = '~'.join([str(l_brace),str(r_brace)])
                entry_name = '.'.join(['ILL',entry.full_name,numbs])
                fasta_dict[entry_name] = full_sequence

        #For Entries that do not have all four fields, and have a translatable sequence
        elif all([entry.tl_start,entry.tl_stop]):
            if entry.strand == '+':
                left_boundary = entry.tr_start if entry.tr_start else entry.tl_start
                right_boundary = entry.tr_stop if entry.tr_stop else entry.tl_stop
                full_sequence = genome_sequence[left_boundary-1:right_boundary]
                l_brace = (entry.tl_start - left_boundary)+1
                r_brace = (entry.tl_stop - entry.tl_start)+l_brace
                numbs = '~'.join([str(l_brace),str(r_brace)])
                entry_name = '.'.join(['ILL',entry.full_name,numbs])
                fasta_dict[entry_name] = full_sequence

            elif entry.strand == '-':
                left_boundary = entry.tr_stop if entry.tr_stop else entry.tl_stop
                right_boundary = entry.tr_start if entry.tr_start else entry.tl_start
                full_sequence = str(Seq(genome_sequence[left_boundary-1:right_boundary]).reverse_complement())
                l_brace = (right_boundary - entry.tl_start)+1
                r_brace = (entry.tl_start - entry.tl_stop)+l_brace
                numbs = '~'.join([str(l_brace),str(r_brace)])
                entry_name = '.'.join(['ILL',entry.full_name,numbs])
                fasta_dict[entry_name] = full_sequence

    #End
    return fasta_dict

#Main Function
def main():
    parser = argparse.ArgumentParser(description='RockHopper to Fasta')
    parser.add_argument('transcripts', type=str, help='Transcripts File')
    parser.add_argument('fasta', type=str, help='Genomic Fasta File')
    parser.add_argument('-name', type=str, default= 'out_transcripts.fa',help='Outfile Name')
    args = parser.parse_args()
    
    #Get Bacterial Chromosome
    in_fasta = read_in_fasta(args.fasta)
    chromosome = in_fasta[in_fasta.keys()[0]]

    #Read In Transcript Info
    transcripts = read_in_rockhopper_transcripts(args.transcripts)
    
    #Generate Transcripts
    out_data = generate_transcripts(transcripts,chromosome)
    
    #Write out
    write_fasta(out_data,args.name)


if __name__ == '__main__': 
    main()
