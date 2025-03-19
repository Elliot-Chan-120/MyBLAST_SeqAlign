import json
import pandas as pd
import plotly.express as px

# NEED TO ADD IN ERROR HANDLING IN CASE THERE ARE NO TARGET FILES

class DataProcessor:
    """This class will work with the data parsed and burned into Data_Storage"""
    def __init__(self, folder_name="Data_Storage"):
        self.folderpath = folder_name

        self.parameters = {}
        self.fasta_seqs = {}  # key = accession, value = seq
        self.alignment1 = []  # list of two dicts: general info and hsp1
        self.alignment10info = [] # list of dicts of top 10 alignments' general info and hsp1 info

    def extract_all_data(self):
        self.get_parameters()
        self.get_fasta()
        self.get_alignment1_json()
        self.get_alignment10_csv()

    def get_parameters(self):
        """
        Acquires the parameters.json information and saves it as a class variable
        """
        parameter_filename = "/Parameters.json"
        with open(self.folderpath + parameter_filename, 'r') as handle:
            parameter_data = json.load(handle)
        self.parameters = parameter_data

        return self.parameters

    def get_fasta(self):
        """
        Saves all HSP (fasta-formatted) hit sequences as a class variable
        """
        fasta_filename = "/ALL_FASTA_storage"
        fastadct = {}
        with open(self.folderpath + fasta_filename, 'r') as fasta_file:
            for line in fasta_file:
                sequence = []
                line = line.strip()
                if line.startswith(">"):
                    label = line
                    fastadct.update({label: sequence})
                else:
                    fastadct[label].append(line)
        self.fasta_seqs = fastadct

        return self.fasta_seqs

    def get_alignment10_csv(self):
        """
        Saves the csv file's information on all hits as a class variable
        """
        multiple_hits_filename = "/All_hits_info.csv"
        with open(self.folderpath + multiple_hits_filename, 'r') as handle:
            all_hit_data = pd.read_csv(handle)
            # you can set data to a class variable in order to save it as a dataframe
        self.alignment10info = all_hit_data.to_string()

        return self.alignment10info

    def get_alignment1_json(self):
        """
        Saves Alignment_1.json file information as a class variable
        """
        alignment1_filename = "/Alignment_1.json"
        with open(self.folderpath + alignment1_filename, 'r') as handle:
            alignment_data = json.load(handle)
        self.alignment1 = alignment_data

        return self.alignment1

    def bit_scores(self):
        """Produces a bar chart that displays each accession ID Hit sequence's bit score"""
        df = pd.read_csv(self.folderpath + '/All_hits_info.csv')
        target_columns = ['Accession', 'HSP_Bit_Score']
        df_filtered = df[target_columns]

        fig = px.bar(df_filtered, y = 'HSP_Bit_Score', x = 'Accession', text_auto='.2s',
                     title = "HSP Bit Score across all Accession ID Hits")

        fig.show()

    def alignment_quality(self):
        """Produces a bubble chart where x is the Bit Score, y is % Identity, the size reflects alignment length and the color intensity indicates the E-value"""
        df = pd.read_csv(self.folderpath + '/All_hits_info.csv')

        fig = px.scatter(df,
                         x='HSP_Bit_Score',
                         y='HSP_%_Identities',
                         size='Alignment_Length',
                         color='HSP_E_Value',
                         opacity=0.6,
                         hover_data='Accession',
                         labels = {'HSP_E_Value': 'E Value', 'HSP_Bit_Score': 'Bit Score', 'HSP_Score': 'Score','HSP_%_Identities': '% Identity', 'Alignment_Length': 'Length'},
                         title="Bubble Chart Plot of HSP Alignment Quality Summary",
                         color_continuous_scale="Viridis")

        fig.show()

    def test(self):
        """Run this to output all the saved data from the files in the terminal: I mainly used this to test if the data was being properly stored"""
        for key, value in self.parameters.items():
            print(f"{key}: {value}")

        for key, value in self.fasta_seqs.items():
            print(f"{key}: {value}")

        for item in self.alignment1:
            for key, value in item.items():
                print(f"{key}: {value}")

        print(self.alignment10info)
