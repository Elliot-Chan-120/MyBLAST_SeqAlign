from MyBlast_Toolkit import BlastData
from Data_Processor import DataProcessor

def test_blast_toolkit():
    test = BlastData("test_file", 20)
    test.parse()
    test.all_info_txt_export()
    test.database_save()

def test_data_procesor():
    demo = DataProcessor()
    demo.extract_all_data()
    demo.bit_scores()
    demo.alignment_quality()
    