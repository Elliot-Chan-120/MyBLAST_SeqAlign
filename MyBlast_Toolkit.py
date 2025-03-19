from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import json
import pandas as pd
from pathlib import Path


class BlastData:
    def __init__(self, filename, hitlist_size, biotype='auto'):
        """
        This class takes a fasta-formatted file containing a biological sequence, and runs it against NCBI Blast Databases
        depending on the biotype which is provided either by the user or automatically detected by leaving the biotype param blank
        """
        self.filename = filename
        self.biotype = biotype
        self.size = hitlist_size

        self.E_threshold = 0.001  # default expect value
        if self.biotype.lower() == "protein":
            self.program = "blastp"
            self.database = "nr"
        elif self.biotype.lower() in ["dna", "rna", "nucleotide"]:
            self.program = "blastn"
            self.database = "nt"
        #  conditional below sets search parameters based on what it reads from the file presented
        elif self.biotype.lower() == 'auto':
                print("Detecting Biotype...")
                path = Path(self.filename)
                with open(path, 'r') as f:
                    fasta_file = [l.strip() for l in f.readlines()]

                protein_specific = 'ILKMFPSHQRNVYWED'
                nucleotide_list = 'ATCGU'
                full_string = ''

                for line in fasta_file:
                    if '>' in line:
                        label = line
                    else:
                        full_string += line.upper()

                protein_count = 0
                nucleotide_count = 0
                for char in full_string:
                    if char in protein_specific:
                        protein_count += 1
                    elif char in nucleotide_list:
                        nucleotide_count += 1

                if nucleotide_count / len(full_string) > 0.8:  # if more than 80% are valid nucleotides -> allows for ambiguity cases
                    self.biotype = 'Nucleotide'
                    self.program = "blastn"
                    self.database = "nt"
                    print(f"Biotype {self.biotype} Detected -> Search Parameters Set")
                elif protein_count > 0:
                    self.biotype = 'Protein'
                    self.program = "blastp"
                    self.database = "nr"
                    print(f"Biotype {self.biotype} Detected -> Search Parameters Set")
                else:
                    print("Could not confidently determine biotype: please manually set biotype\nOptions: NUCLEOTIDE or PROTEIN")

        self.parameters = {}        # dict

        # BEST MATCH INFO
        self.alignment1_info = {}   # dict
        self.hsp1_info = {}         # dict -> later will combine with alignment1_info
        self.raw_alignment1_hsp_data = {}  # also dict

        # TOP 10 INFO
        self.alignmentmultiple_info = []  # list of dicts
        self.fasta_seqs = []

    def parse(self):
        """
        Initiates function cascade that eventually parses information regarding the sequence hits
        and saves it into class variables
        """
        print(f"[1/4] Searching for sequence matches using {self.program} on the {self.database} database...")

        # E_threshold = 0.001  # expect value -> describes number of hits expected by chance
        # lower values will indicate more significant matches

        # if accepting raw sequence -> convert it into fasta format before sending it to BLAST
        record = SeqIO.read(open(self.filename), format="fasta")
        handle = NCBIWWW.qblast(self.program,
                                self.database,
                                record.format("fasta"),
                                expect = self.E_threshold,
                                hitlist_size=self.size)  # determines max number of sequence matches (hits) returned

        save_file = open("test_output_file.xml", "w")  # outputs information into xml file
        save_file.write(handle.read())
        save_file.close()
        handle.close()

        self.data_check()

    def data_check(self):
        """
        Checks that enough alignments have been gathered to satsify the requested hitsize
        \n if there aren't enough, then adjust the expect value and rerun the search
        :return:
        """
        result_handle = open("test_output_file.xml")
        blast_record = NCBIXML.read(result_handle)

        # results have been parsed, but before you try to process the alignments
        # check in case there aren't enough alignments in the first run and increase E-value
        while len(blast_record.alignments) < self.size and self.E_threshold < 10:
            self.E_threshold = self.E_threshold*10
            print(f"[!] Not enough matches found: running with E threshold at {self.E_threshold}...")
            result_handle.close()
            self.parse()
            if self.E_threshold == 10 and len(blast_record.alignments) < self.size:
                print("[X] There are no meaningful matches available for the provided sequence")
                return None
        else:
            print("[2/4] Obtaining data...")
            self.fetch_info()

    def fetch_info(self):
        """Saves information to class variables"""
        result_handle = open("test_output_file.xml")
        blast_record = NCBIXML.read(result_handle)
        print("[3/4] Parsed XML file obtained --> Saving to class variables...")

        # QUERY ID AND PARAMETER SETTINGS
        self.parameters = {
            "ID": blast_record.query_id,
            "Database": blast_record.database,
            "Matrix": blast_record.matrix,
            "Gap_Penalties": blast_record.gap_penalties
        }


        # BEST ALIGNED SEQUENCE INFORMATION
        best_alignment = blast_record.alignments[0]  # access best aligned sequence
        self.alignment1_info = {
            "Accession": best_alignment.accession,
            "ID": best_alignment.hit_id,
            "Definiton": best_alignment.hit_def,
            "Alignment_Length": best_alignment.length,
            "HSP count": len(best_alignment.hsps)
        }

        # HighestSP INFORMATION DICT
        hsp = best_alignment.hsps[0]  # access highest scoring pair
        self.hsp1_info = {
            "E_value": hsp.expect,
            "Score": hsp.score,
            "Bit_Score": hsp.bits,
            "Identities": hsp.identities,
            "%_Identities": (hsp.identities / hsp.align_length) * 100,
            "Mismatches": hsp.align_length - hsp.identities,
            "Gaps": hsp.gaps,
            "Query": hsp.query,
            "Match": hsp.match,
            "HitSq": hsp.sbjct
        }

        # INFORMATION FOR ALL HIT SEQUENCES
        for alignment in blast_record.alignments:   # basic information obtained here
            align_dict = {
                "Accession": alignment.accession,
                "Definition": alignment.hit_def,
                "Alignment_Length": alignment.length,
            }
            for hsp in alignment.hsps:              # hsp information obtained here
                count = 0
                align_hsp_dict = {"HSP #": count + 1,
                                  "HSP_Bit_Score": hsp.bits,
                                  "HSP_Score": hsp.score,
                                  "HSP_Identities": hsp.identities,
                                  "HSP_%_Identities": (hsp.identities / hsp.align_length) * 100,
                                  "HSP_E_Value": hsp.expect,
                                  "HSP_mismatches": hsp.align_length - hsp.identities,
                                  "HSP_Gaps": hsp.gaps
                                  }
                align_dict.update(align_hsp_dict)
            self.alignmentmultiple_info.append(align_dict)


        # OBTAIN ALL HIT SEQUENCES
        for alignment in blast_record.alignments:
            fasta_info = {"Accession": alignment.accession}
            for hsp in alignment.hsps:
                seq = {"Sequence": hsp.query}
                fasta_info.update(seq)
            self.fasta_seqs.append(fasta_info)

        print("[4/4] Data Acquired")
        return (self.parameters,
                self.alignment1_info,
                self.hsp1_info,
                self.alignmentmultiple_info,
                self.fasta_seqs)

    def all_info_txt_export(self):
        """
        Generates a txt file containing all the information that was parsed
        :return:
        """
        print("Generating TXT file w/ all parsed info")
        path = Path('Data_Storage') / 'Test_Results.txt'
        content = ""

        content += "[1] [ID & Parameters]"
        for key, value in self.parameters.items():
            content += f"\n{key}: {value}"

        content += f"\n{"=" * 100}\n"  # separator
        content += "[2] [Alignment_1 Data]"
        for key, value in self.alignment1_info.items():
            content += f"\n{key}: {value}"

        content += f"\n{"=" * 100}\n"  # separator
        content += "[3] [Alignment_1's Highest HSP Data]"
        for key, value in self.hsp1_info.items():
            content += f"\n{key}: {value}"

        content += f"\n{"=" * 100}\n"  # separator
        content += "[4] [Top 10 Alignments' Data]"
        for align_dict in self.alignmentmultiple_info:
            content += f"\n{"= " * 30}"
            for key, value in align_dict.items():
                content += f"\n{key}: {value}"

        path.write_text(content)
        print("TXT file generated.")

    def database_save(self):
        """Generates multiple files containing their respective info into the data_storage folder"""
        folderpath = Path('Data_Storage')
        print(f"Initiating Data Save in folder [{folderpath}]")

        tmp_list = []  # save everything about alignment 1 into a temporary list which we can then access later
        tmp_list.append(self.alignment1_info)
        tmp_list.append(self.hsp1_info)

        print("[1/4] Storing Parameters as Parameters.json...")
        with open(folderpath / "Parameters.json", 'w') as file:     # generate json file for parameters
            json.dump(self.parameters, file, indent=4)
        print("Parameters.json saved")

        print("[2/4] Storing Alignment 1 Data as Alignment_1.json...")
        with open(folderpath / 'Alignment_1.json', 'w') as file:    # generate json fiel for alignment[0] info
            json.dump(tmp_list, file, indent=4)
        print("Alignment_1.json saved")

        print("[3/4] Storing Info on All Hits as All_hits_info.csv...")
        df = pd.DataFrame(self.alignmentmultiple_info)              # generate csv file for top10alignment info
        df.to_csv(folderpath / 'All_hits_info.csv')
        print("All_hits_info.csv saved")

        print("[4/4] Storing all HSP sequences as All_FASTA_storage...")
        seqs = ""
        for hsp_dict in self.fasta_seqs:
            accession = hsp_dict["Accession"]
            sequence = hsp_dict["Sequence"]
            seqs += f"> {accession}\n{sequence}\n"
        (folderpath / 'All_FASTA_storage').write_text(seqs)
        print("All_FASTA_storage saved")

        # generate fasta compilation file for all seqs
        # this is going to be in order from the greatest similarity to the least

    def test(self):
        print("parameters")
        for key, value in self.parameters.items():
            print(key, value)

        print("alignment1")
        for key, value in self.alignment1_info.items():
            print(key, value)

        print("hsp1")
        for key, value in self.hsp1_info.items():
            print(key, value)

        print("List of all alignments")
        for item in self.alignmentmultiple_info:
            print(item)

        print("raw alignment1 data")
        for key, value in self.raw_alignment1_hsp_data.items():
            print(key, value)
