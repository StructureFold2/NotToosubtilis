#!/usr/bin/env python

'''
Separates <.fasta> or <.react> of any sequences following a given nomenclature into 
Five Prime, Three Prime, CDS, and intronic (or intergenic) segments/files.
'''

#Imports
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse

#Functions
def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def write_out_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for name,seq in sorted(info.items()):
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n')

def read_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) if x!= 'NA' else 'NA' for x in reactivities.split()]
    return information

def write_react(react_info,outfyle='out.react'):
    '''Writes the <.react> back out'''
    with open(outfyle,'w') as g:
        for transcript, data in react_info.items():
            g.write(transcript+'\n')
            g.write('\t'.join([str(i) for i in data])+'\n')

def set_mode(infile):
    '''returns a mode given an input file'''
    fasta_extensions = {'fa':'fasta','fasta':'fasta','fas':'fasta','fna':'fasta','fnn':'fasta','frn':'fasta'}
    react_extensions = {'react':'react'}
    directions = dict(fasta_extensions.items()+react_extensions.items())
    extension = infile.split('.')[-1]
    try:
        fyle_type = directions[extension]
        return fyle_type
    except KeyError:
        return None

def generate_introns(fancy_tuple):
    ''''''
    index,pseudos = 1,[]
    for i in range(0,len(fancy_tuple)-1):
        if fancy_tuple[i+1][0] - fancy_tuple[i][1] > 1:
            pseudo_intron = (fancy_tuple[i][1]+1,fancy_tuple[i+1][0]-1,'_'.join(['INTRON',str(index)]),'INTRON')
            index+=1
            pseudos.append(pseudo_intron)
    return pseudos

def parse_indexes(info_dict):
    '''Generate Indexes based on names and lengths of items'''
    master_indexes = {}
    for k,v in info_dict.items():
        tag,identifier,coords = k.split('.')
        id_bits,coord_bits = identifier.split('_'),[[int(q) for q in x.split('~')] for x in coords.split('_')]
        fused = [tuple(a+[b,'CDS']) for a,b in zip(coord_bits,id_bits)]
        FPUTR = [] if fused[0][0] == 1 else [(1,fused[0][0]-1,'FPUTR','FPUTR')]
        TPUTR = [] if fused[-1][1] == len(v) else [(fused[-1][1]+1,len(v),'TPUTR','TPUTR')]
        introns = generate_introns(fused)
        fused = sorted(FPUTR+fused+TPUTR+introns)
        master_indexes[k] = fused
    return master_indexes

def sunder(fasta_dict,indexes,zbuff):
    '''Return a new dictionary of components based on the indexes'''
    #Buckets
    CDS,FPUTR,TPUTR,INTRONS = {},{},{},{}
    #
    for long_name,sequence in fasta_dict.items():
        seq_indexes = indexes[long_name]
        for index in seq_indexes:
            #
            feature_name = long_name+'.'+''.join([str(index[0]),'~',str(index[1]),'_',index[2],'-',index[3]])
            feature_tag = feature_name.split('-')[-1]
            
            l_index,r_index = index[0]-1-zbuff, index[1]+zbuff
            l_index = l_index if l_index >= 0 else 0
            feature_seq = sequence[l_index:r_index]
            #feature_name = long_name+'.'+'~'.join([str(y) for y in index])
            if feature_tag == 'CDS':
                CDS[feature_name] = feature_seq
            elif feature_tag == 'FPUTR':
                FPUTR[feature_name] = feature_seq
            elif feature_tag == 'TPUTR':
                TPUTR[feature_name] = feature_seq
            elif feature_tag == 'INTRON':
                INTRONS[feature_name] = feature_seq
            else:
                print 'UNRECOGNIZED TYPE'
                print index, feature_seq,feature_name,feature_tag
    return [FPUTR,CDS,TPUTR,INTRONS]

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Splits <.fasta> or <.react> files into bacterial components.')
    parser.add_argument('infile',type=str,help='Input <.react> or <.fasta> file.')
    parser.add_argument('-buff',type=int,default=0,help='[default = 0] Adds n bp to the ends of all features.')
    parser.add_argument('-minimum',type=int,default=10,help='[default = 10] Minimum Feature Length after buffer.')
    parser.add_argument('-restrict',default = None, help = 'Limit analysis to these specific transcripts <.txt> ')
    args = parser.parse_args()
    
    #Get file type and thus mode
    mode = set_mode(args.infile)

    #Fasta Suite
    if mode == 'fasta':
        sequences = read_fasta(args.infile)
        
        if args.restrict != None:
            covered = get_covered_transcripts(args.restrict)
            filter_dictonary(sequences,covered)
        
        indexes = parse_indexes(sequences)
        new_sequences = sunder(sequences,indexes,args.buff)
        new_names = ['_'.join([args.infile.split('.')[0],extension,str(args.buff)+'buff','only.fa']) for extension in ['fp','cds','tp','intron']]
        
        #Filter Information
        for dictionary in new_sequences:
            for name, info in dictionary.items():
                if len(info) < args.minimum:
                    del dictionary[name]
        #Write out
        for pair in zip(new_sequences,new_names):
            write_out_fasta(*pair)

    #React Suite
    elif mode == 'react':
        reactivities = read_reactivities(args.infile)

        if args.restrict != None:
            covered = get_covered_transcripts(args.restrict)
            filter_dictonary(reactivities,covered)
        
        indexes = parse_indexes(reactivities)
        new_reacts = sunder(reactivities,indexes,args.buff)
        new_names = ['_'.join([args.infile.split('.')[0],extension,str(args.buff)+'buff','only.react']) for extension in ['fp','cds','tp','intron']]
        
        #Filter Information
        for dictionary in new_reacts:
            for name, info in dictionary.items():
                if len(info) < args.minimum:
                    del dictionary[name]
        
        #Write out
        for pair in zip(new_reacts,new_names):
            write_react(*pair)

    #Null Suite
    else:
        print 'Invalid input. {} is neither a <.fasta> nor a <.react> file. Try again.'.format(args.infile)



if __name__ == '__main__': 
    main()


