# NotToosubtilis
Utility scripts for generating a newly annotated _B. subtilis_ transcriptome, but very transferable to any prokaryotic genome.

## RockHopper_extract.py
This applies an annotation generated by RockHopper to the genome it was generated against, creating a fasta file of 
transcripts with newly annotated 5' and 3'UTRs appended. New transcripts have specific nomenclature such that it is easy to work 
with them in any capacity. For example, the transcript ILL.yoaE-BSU18570.14~2056:

+ Was generated with Illumina Data
+ Contains the protein 'yoaE', identified as BSU18570
+ The protein resides in bases 14~2056, anything before is 5'UTR, anything after is 3'UTR

### Usage

```
positional arguments:
  transcripts  Transcripts File
  fasta        Genomic Fasta File

optional arguments:
  -h, --help   show this help message and exit
  -name NAME   Outfile Name
```

## transcript_parser.py
This will break apart files generated by **RockHopper_extract.py** into sub-fastas of the respective transcript regions. It will
also work on <.react> files generated via StructureFold2 mapped against a fasta generated by **RockHopper_extract.py**. Basically, 
if you want to separate sequences or reactivity values into separate files for transcript components (5'UTR, CDS, 3'UTR) for any number
of reasons, this script will autmatically do it provided the transcript nomenclature was generated by **RockHopper_extract.py**.

As an additional option, a buffering option is inculuded, i.e., affix n bases to each applicable end of the transcript feature. Thus if you wanted
to study the 5'UTR and have the first start codon, one would assign a buffer of 3. 

### Usage

```
Splits <.fasta> or <.react> files into bacterial components.

positional arguments:
  infile              Input <.react> or <.fasta> file.

optional arguments:
  -h, --help          show this help message and exit
  -buff BUFF          [default = 0] Adds n bp to the ends of all features.
  -minimum MINIMUM    [default = 10] Minimum Feature Length after buffer.
  -restrict RESTRICT  Limit analysis to these specific transcripts <.txt>

```
