from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import sys
import time

#TODO CHANGE MAX TARGET at line 70! Detaul set at MAX 5000line (take very long)


# remember to input email for Entrez
Entrez.email = "kit.liew@outlook.com"


def get_organism_firstrow(html_text):
    """Return a tuple(organism, first_row_of_organism) from blastp taxonomy tabs, cleanup html text"""

    organisms = []
    first_match = []

    contents = html_text

    # seperate into organisms group
    content = contents.split("<tr>")[1:]

    for i in content:
        # seperate organism name from list of accessions
        organism, blast_result = i.partition("</tr>")[0::2]

        # cleanup organism html
        organism = organism.split("</span>")[0].split('<span class="orgTaxHd">')[1].split(">")[1].replace("</a", "").split("(")[0].strip()

        # get the first row after organism name
        first_row = blast_result.split('<tr class="dflLnk">')[1].partition('<td class="c1 l lim">')[-1].split(">")[1].replace("</a", "")

        organisms.append(organism)
        first_match.append(first_row)

    return organisms, first_match


def blastp_scraping(query_protein_id):
    """Automatically Blastp from web browser/chrome, scrap html result and output (organisms, first row accession)"""

    # Setting up web browser
    # remember to download chromedriver and set path
    PATH = "C:\Program Files (x86)\Google\Chrome\Application\chromedriver.exe"
    # init chrome
    driver = webdriver.Chrome(PATH)
    # go to ncbi page
    driver.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins")

    # find query box and input protein_id
    query_box = driver.find_element_by_name("QUERY")
    query_box.clear()
    query_box.send_keys(query_protein_id)
    time.sleep(2)

    # algorithm parameters
    # search only Bacteria
    organism_box = driver.find_element_by_id("qorganism")
    organism_box.clear()
    organism_box.send_keys("Bacteria (taxid:2) ")
    # number of max target sequence
    driver.find_element_by_id("algPar").click()
    driver.find_element_by_xpath("//*[@id='NUM_SEQ']/option[7]").click()    # {1:10, 2:50, 3:100, 4:250, 5:500, 6:1000, 7:5000} default at 3

    # submit blastp
    driver.find_element_by_id("blastButton1").click()
    driver.back()
    driver.find_element_by_id("b1").click()

    try:
        element = WebDriverWait(driver, 600).until(EC.presence_of_element_located((By.ID, "btnTaxn")))
        element.click()
        element = WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.ID, "btnTaxOrg")))
        element.click()
        time.sleep(5)

        # get HTML format
        table = driver.find_element_by_xpath('//*[@id="orgTable"]/tbody')
        html_text = table.get_attribute("innerHTML")

        return get_organism_firstrow(html_text)

    except IOError:
        print("Problem connecting to NCBI")


def get_nucleotideid_accession(accession_list):
    """fetch refseq accession from protein db using protein id"""
    result = []
    handle = Entrez.efetch(db="protein", id=accession_list, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, format="gb")

    for record in records:

        # get the last feature, CDS
        cds = record.features[-1]

        if str(cds.type) == "CDS":

            # cleanup. Only get the cds nucleotide id
            cds_feature = str(cds).split("Key: coded_by, Value: ")[1].split("\n")[0].replace("[\'", "").replace("complement(", "").replace(r"join(", "").replace(")", "").replace("\']", "")

            # cds
            n_id, cds_range = cds_feature.split(":")
            result.append(n_id)

            #TODO for future development
            # # cleanup cds_range
            # start, end = cds_range.split("..")
            # cds_range = str(start) + ":" + str(int(end) + 1)
            # nucleotide_ids.append(n_id)
            # cds_regions.append(cds_range)
            # print(n_id, cds_range)

        else:

            result.append([])

    return result


def get_nucleotide_fasta(accession_list):
    """fetch fasta sequence from protein db using protein id"""
    handle = Entrez.efetch(db="nucleotide", id=accession_list, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    return [str(i.seq) for i in records]


def build_df(input_csv, output_csv, organisms_list, fasta_sequences, gc_fasta):

    # read the input.csv file
    df = pd.read_csv(input_csv)

    # compile data together
    data = zip(organisms_list, fasta_sequences, gc_fasta)

    # make dataframe from blast results
    df1 = pd.DataFrame(data, columns= ["Organism Name", "Fasta", "GC_content"])

    # cleanup ["Organism Name"] column data, case sensitive
    df["Organism Name"] = df["Organism Name"].str.lower()
    df1["Organism Name"] = df1["Organism Name"].str.lower()

    # merged two dataframe on Organism Name
    df_merged = pd.merge(df, df1, how="left", on="Organism Name")

    return df_merged


if __name__ == "__main__":
    query_protein_id = str(sys.argv[1])
    input_csv = sys.argv[2]
    output_csv = sys.argv[3]

    # start with Blastp the protein
    organisms_list, first_match_list = blastp_scraping(query_protein_id)
    print("Number of Organisms and high score Blastp: " + str(len(organisms_list)))

    # find the nucleotide/cds reference accession
    nucleotide_accession_list = get_nucleotideid_accession(first_match_list)
    print("Number of Nucleotide Accession: " + str(len(nucleotide_accession_list)))

    # get the fasta sequences
    print("Getting Fasta sequence...")
    fasta_sequences = get_nucleotide_fasta(nucleotide_accession_list)

    # calculate the GC content
    print("Calculating GC content...")
    gc_fasta = [GC(seq) for seq in fasta_sequences]

    # build the dataframe
    print("Comparing to {}  ...".format(input_csv))
    df = build_df(input_csv, output_csv, organisms_list, fasta_sequences, gc_fasta)
    print("A list of DataFrame Columns: ", str(df.columns))
    print("Showing few rows that have matched result...")
    print(df[df.Fasta.notnull()])

    # write to output_csv
    print("Writing to {}...".format(output_csv))
    df.to_csv(output_csv, sep=",", index = False)
    print("Write complete. Script Ends!")
