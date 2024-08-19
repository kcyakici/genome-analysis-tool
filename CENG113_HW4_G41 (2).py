# CENG113 HW4 TEMPLATE
# AUTHOR: SAMET TENEKECI
# DATE: 27/12/2021

# MODIFY AND RENAME THIS FILE FOR YOUR ASSIGNMENT
# EXAMPLE FILE NAMING: CENG113_HW4_G01.py OR CENG_113_HW4_G25.py
# WRITE STUDENT IDS, NAMES & SURNAMES OF THE GROUP MEMBERS AT THE TOP
# SUBMIT ONLY THIS FILE AND ONLY ONCE FOR YOUR GROUP

#Deniz Kaya 280201033
#Kürşat Çağrı Yakıcı 290201098



def read_genes(file_path):
    genes_dict = {}
    lines = []
    with open(file_path) as f:
        for line in f.readlines():  #Stripping \n and > s from lines so it will be in the desired format.
            line = line.strip() 
            line = line.strip(">")
            lines.append(line)
    for i in range(len(lines)):  #Since the data contains key in even lines values in odd lines, adding to dictionary according to their index is easier.
        if i % 2 == 0:
            genes_dict[lines[i]] = lines[i+1]
    return genes_dict

def get_fragments(genes_dict, frag_len = 50):
    genes_dict_keys = genes_dict.keys()     #Getting the keys of genes to assign them to their gene count.
    gene_counts = {}
    frag_dict = {}

    for key in genes_dict_keys:  #Spilitting the key to its components to find their gene counts and assign them back to its key in a new dict named gene_counts.
        numbers = key.split("|") 
        numbers = numbers[1].split("-")
        number1 = int(numbers[0])
        number2 = int(numbers[1])
        gene_counts[key] = number2 - number1

    for key in genes_dict_keys:
        subset_count = 0
        if gene_counts[key] >= frag_len:
            subset_count = gene_counts[key] // frag_len #Finding how many subset the code will create from that gene which we will use later in a for loop.
            components_of_key = key.split("|")  #Spilitting the components of key because we will rename each key to contain 50 difference between numbers.
            name = components_of_key[0]
            numbers = components_of_key[1].split("-")
            number1 = int(numbers[0])
            genes_of_key = genes_dict[key]
            for subset in range(subset_count):
                gene_part = ""
                for gene_index in range(subset*50,subset * 50 + 50): #Adding first 50 gene one by one to new key.
                    gene_part += genes_of_key[gene_index]
                key_part = name + "|" + str(number1) + "-" + str(number1+50) #Creating the name of the key by using first number and last number as 50 more of it.
                frag_dict[key_part] = gene_part
                number1 += 50  #Incrementing the number by 50 to start the other subset if needed.
    return frag_dict

def filter_frags(frag_dict, threshold = 0.7):
    def get_similarity(s1,s2):
        similar_chr_count = 0
        total_chr_count = len(s1)
        for i in range(len(s1)): #Finding how many similar character the both genes.
            if s1[i] == s2[i]:
                similar_chr_count += 1
        
        match_percentage = similar_chr_count/total_chr_count  #Dividing the similar to total to find percentage.
        return match_percentage
    frag_dict_keys = list(frag_dict.keys())  #Turning the keys to a list and then sorting it for it to become lexicographicly sorted.
    frag_dict_keys = sorted(frag_dict_keys)  
    dissimilar_frag_dict = dict()
    for elem in frag_dict_keys:  #Adding every element of frag_dict in order of sorted version.
        dissimilar_frag_dict[elem] = frag_dict[elem]
    for key_index in range(len(frag_dict_keys)): #Working with their indexes to only compare the first element with the elements after it to not compare the two element two times and delete them both.
        key = frag_dict_keys[key_index]          #Impelemting an algorithm that works like bubble sort. First key index is the element that will be compared with all other key values in the dictionary.
        key_value = frag_dict[key]               #Other key index is the variable that goes through all the elements after the key index to compare it with them.
        for other_key_index in range(key_index,len(frag_dict_keys)):
            other_key = frag_dict_keys[other_key_index]
            other_key_value = frag_dict[other_key]
            if key == other_key:
                continue
            elif get_similarity(key_value,other_key_value) >= threshold: #If the similarity is bigger than the threshold, delete the element from dissimilag frag dict.
                if other_key in dissimilar_frag_dict:
                    del dissimilar_frag_dict[other_key]
    return dissimilar_frag_dict

def get_sentences(dissimilar_frag_dict):
    def generate_kmers(seq,k):
        kmer = ""
        for i in range(len(seq)-k+1): #Number of k length genes.
            if i != 0: #Create a space after each kmer.
                kmer += " "
            for a in range(i,i+k): #Add each element k times starting from i.
                kmer += seq[a]
        return kmer
    
    dissimilar_frag_dict_keys = dissimilar_frag_dict.keys()
    sentences_dict = {}
    for key in dissimilar_frag_dict_keys:  #Take the key and value of an element.
        key_value = dissimilar_frag_dict[key]
        _4mer = generate_kmers(key_value,4) #Generate a kmer from genes.
        sentences_dict[key] = _4mer #Assign the value with the same key to another dictionary named sentences_dict.
    return sentences_dict

def clean_dict(sentences_dict):
    def clean_sentence(sentence):
        split_vers = sentence.split() #Split the sentence to compare them with each other.
        unique = ""
        for elem in split_vers: #Adding the 4mers to the unique if unique does not have the same elem.
            if elem not in unique:
                unique += elem
                unique += " "
        unique = unique[:-1] #To delete the last space in each unique sentence.
        return unique
    sentences_dict_keys = sentences_dict.keys()
    clean_sentences_dict = {}

    for key in sentences_dict_keys: #First find the key and then clean the value of the key and reassing the value to the same key in a new dictionary.
        key_value = sentences_dict[key]
        clean_vers = clean_sentence(key_value)
        clean_sentences_dict[key] = clean_vers
    return clean_sentences_dict


def write_genes(file_path,clean_sentences_dict):
    import csv
    header = ["fragment_id","sentence","sentence_length","number_of_words"]
    with open(file_path, "w",newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header) #Write header to the top.
        clean_sentences_dict_keys = clean_sentences_dict.keys() #Get the keys of dictionary.
        for key in clean_sentences_dict_keys: #Give the value to the variables for all dictionary elements and then write them to the csv file one by one.
            fragment_id = key
            sentence = clean_sentences_dict[key]
            sentence_lenght = len(sentence)
            number_of_words = len(sentence.split())
            row = [fragment_id,sentence,sentence_lenght,number_of_words]
            writer.writerow(row)

def main():
    genes_dict = read_genes("input.txt")
    print(len(genes_dict))

    frag_dict = get_fragments(genes_dict, frag_len = 50)
    print(len(frag_dict))

    dissimilar_frag_dict = filter_frags(frag_dict, threshold = 0.7)
    print(len(dissimilar_frag_dict))

    sentences_dict = get_sentences(dissimilar_frag_dict)
    sentence = sentences_dict["chrX|105727750-105727800"] #Just a random element in sentence to measure its lenght.
    split_vers = sentence.split()
    print(len(split_vers)) #Printing the length of split version of the sentence.


    #Finding the values before clean procedure and after to compare them two and print the result.
    clean_sentences_dict = clean_dict(sentences_dict)
    clean_sentences_dict_values = clean_sentences_dict.values()

    sentences_dict_values = sentences_dict.values()
    all_values_before = ""

    for value in sentences_dict_values:
        all_values_before += value
        all_values_before += " "

    all_values_after = ""

    for value in clean_sentences_dict_values:
        all_values_after += value
        all_values_after += " "
    
    all_values_before = all_values_before.split() #Split them to only find how many kmers in them.
    all_values_after = all_values_after.split()
    print(len(all_values_before) - len(all_values_after))

    write_genes("output.csv",clean_sentences_dict)


main()
