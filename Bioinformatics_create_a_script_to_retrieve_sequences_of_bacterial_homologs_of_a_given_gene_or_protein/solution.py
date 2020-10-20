from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import sys
import re
import time

# remember to input email for Entrez
Entrez.email = "kit.liew@outlook.com"

def search_term(gene, organism_list):
    """ combine gene symbol and search organism
        e.g iles AND Coxiella endosymbiont of Amblyomma americanum[ORGN] """
    result = []
    for organism in organism_list:
        term = str(gene + " AND " + organism + "[ORGN]")
        result.append(term)
    return result


def get_id(term_list):
    """Return first GeneID from search"""

    result = []

    for term in term_list:

        # search GENE AND ORGN from gene database
        handle = Entrez.esearch(db="gene", term=term)
        records = Entrez.read(handle)

        # check if any records
        if records["IdList"]:

            current_ID = records["IdList"][0]

            result.append(current_ID)

        # if search no result then add empty1
        else:
            print(term)
            result.append("")

    return result


def nucleotide_id_and_seq_range(gene_id_list):
    """Return (gene_id, nucleotide_id, seq_range)"""

    result = []

    # fetch summary for annotation(nucleotide id) #shows seq range and complement?
    handle = Entrez.efetch(db="gene", id=gene_id_list, rettype="gb", retmode="text")
    records = handle.read()

    # records are separated by newline
    test = records.strip().split("\n\n")
    for i in test:

        # Get ID
        gene_id = i.split("ID: ")[-1]

        if "Annotation:" in i:

            annotation_range = i.split("Annotation: ")[1].split("\n")[0].strip()

            if "Chromosome " in annotation_range:

                annotation_range = re.sub("Chromosome [0-9a-zA-Z] ", "", annotation_range)

            annotation_id, seq_range = annotation_range.partition(" ")[::2]
            seq_range = seq_range_clean(seq_range)
            result.append((gene_id, annotation_id, seq_range))

    return result


def efetch_fasta_seq(nucleotide_ids_seq_range):
    """Return fasta/nucleotide sequence"""

    results = []

    # first element is Gene ID
    # second element in tuple is Nucleotide Id
    nucleotide_ids = [item[1] for item in nucleotide_ids_seq_range]
    # second element in tuple is sequence range
    seq_range = [item[2] for item in nucleotide_ids_seq_range]

    handle = Entrez.efetch(db="nucleotide", id=nucleotide_ids, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")

    test = zip(seq_range, records)

    for element in test:
        start = element[0][0]
        end = element[0][1]
        condition = element[0][2]

        record = element[1]

        # get exact fasta seq
        sequence = record.seq[start:end]

        # reverse complement if True
        if condition:
            sequence = sequence.reverse_complement()
            results.append((record.id, str(sequence), GC(sequence)))

        else:
            sequence = sequence
            results.append((record.id, str(sequence), GC(sequence)))


    return results


def seq_range_clean(seq_range_list):
    """ clean up this format '(109382..111916, complement)' """

    start, end = seq_range_list.split("..")
    start = int(re.sub("[^0-9]", "", start))
    end = int(re.sub("[^0-9]", "", end))

    if "complement" in seq_range_list:
        return (start, end, True)

    else:
        return (start, end, False)


def get_organisms_csv(filename):
    df = pd.read_csv(filename)
    organism_list = df["Organism Name"].tolist()
    return organism_list


def merge_3_data(df, data1, data2, data3, output_csv):
    # create 3 dataframe each have at least one common column
    df1 = df
    df1["Gene_Id"] = data1
    df1["_count"] = df1.groupby("Gene_Id").cumcount()

    df2 = pd.DataFrame(data2, columns=["Gene_Id", "Nucleotide_Id"])
    df2["_count"] = df2.groupby("Gene_Id").cumcount()

    df3 = pd.DataFrame(data3, columns=["Nucleotide_Id", "Fasta_Seq", "GC_Fasta_Seq"])
    df3["_ncount"] = df3.groupby("Nucleotide_Id").cumcount()


    # merge the 3 dataframe
    df_merged = pd.merge(df1, df2, how="left", on=["Gene_Id", "_count"]).drop(columns=["_count"])
    df_merged["Nucleotide_Id"] = df_merged["Nucleotide_Id"]
    df_merged["_ncount"] = df_merged.groupby("Nucleotide_Id").cumcount()
    df_merged1 = pd.merge(df_merged, df3, how="left", on=["Nucleotide_Id", "_ncount"]).drop(columns=["_ncount"])

    return df_merged1


if __name__ == "__main__":
    # test_target = ["LeuS", "IleS", "ValS"]
    #
    # LeuS[624 rows x 4 columns]
    # Summary: Run
    # time = 3110.285308122635
    #
    # Iles
    # [620 rows x 4 columns]
    # Summary: Run
    # time = 2988.8152072429657
    #
    # #Vlas
    # # time = 3139.697003364563
    # [595 rows x 4 columns]
    # Summary: Run

    gene_name = sys.argv[1]
    input_csv = sys.argv[2]
    output_csv = sys.argv[3]

    start = time.time()

    # get organism list
    print("Retrieve organism name from " + str(input_csv))
    organism_list = get_organisms_csv(input_csv)
    # get search term list
    print("Number of Organism: " + str(len(organism_list)))
    term_list = search_term(gene_name, organism_list)       #TODO can run partially if needed

    # get id list
    # this is most time consuming part. Average 1 sec for 1 search
    print("This is most time consuming part. Average 1 sec for 1 search...")
    print("Getting Gene Id...")
    print("Loading.... ............Displaying organisms without result......................")
    list_gene_ids = get_id(term_list)
    print("Search DONE! NOW THINGS are fast!")


    # get annotation_id and sequence_range
    print("Getting Nucleotide Id...")
    annotation = nucleotide_id_and_seq_range(list_gene_ids)

    # get annotation_id, fasta_seq, gc_fasta
    print("Fetching Fasta Sequence and Calculate GC...")
    nucleotide_seq = efetch_fasta_seq(annotation)

    # compile data together
    print("Compile data...")
    df = pd.read_csv(input_csv)
    data1 = list_gene_ids
    data2 = [i[0:2] for i in annotation]
    data3 = nucleotide_seq

    # get final dataframe
    print("Merging data...")
    df_result = merge_3_data(df, data1, data2, data3, output_csv)

    # show number of hits
    print("Whole DataFrame")
    print(df_result)

    print("Added columns and values....")
    print(df_result[df_result.GC_Fasta_Seq.notnull()][["Organism Name", "Gene_Id", "Nucleotide_Id", "Fasta_Seq", "GC_Fasta_Seq"]])

    test = df_result[df_result.GC_Fasta_Seq.notnull()]
    print(str(test.shape) + " hits!")


    # Write to output csv
    df_result.to_csv(output_csv, sep=",", index=False)

    end = time.time()
    print("Summary: Run time = " + str(end - start))
    print("Write complete. Script Ends!")
