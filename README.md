# MyBLAST_SeqAlign
 
This project utilizes two classes designed to process and analyze BLAST (Basic Local Alignment Search Tool) results. Accepting a FASTA-formatted file containing a biological sequence (protein/nucleotide), class BlastData outputs an NCBI BLAST XML file, parses relevant data on a specified number of aligned sequences and stores it in structured formats via JSON, CSV and FASTA (txt). The DataProcessor class then accesses these files and generates summary charts based off the information parsed. Overall, this tool / project aims to facilitate easy access to alignment data for further analysis, visualization techniques and downstream bioinformatics applications.

The default expect value is set to 0.001, but if insufficient matches with this value are received, the value will be automatically adjusted and the program will be rerun with it. Furthermore, if a sequence biotype is not specified by the user, the BlastData class will automatically detect it.

It is recommended that classes are run one at a time, as database queries can take up to 9-10 minutes each run (before I got a timeout error), whereas generating the quality alignment summary chart is much faster.

Dependencies:
- Python 3.8+
- Biopython
- Pandas
- Plotly
