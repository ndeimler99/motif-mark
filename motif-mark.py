#!/usr/bin/env python

#import required modules (argparse for user input, cairo for graphics)
import argparse
import cairo
import re
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt

#########################################################################################################################################
#argparse function to take user input
def get_args():
    """argparse function that takes user input"""
    parser = argparse.ArgumentParser(description = "motif-mark: for easy visualization of the presence of motifs in your sequences")
    #action = append will take any duplicated -f flag and append it to the end of a list, allows for multiple input fasta files
    parser.add_argument("-f", "--files", help = "what is your input fasta files?", required = True, action="append")
    #a txt file containing motifs, each motif must be on its own line
    parser.add_argument("-m", "--motifs", help="What are your motifs? Each motif should appear on a new line", required=True)
    #must be a matplotlib color palette, default contains ten colors
    parser.add_argument("-c", "--color", help="What matplotlib color palette would you like to use", required=False, default="tab10")
    return parser.parse_args()

#assign user input to global script variables
args = get_args()
fasta_files = args.files
motif_file = args.motifs
color = args.color

#########################################################################################################################################
#Dictionariy containing keys of possible letters that could appear in an RNA or DNA sequence with values of their meaning
#In the case of ambiguous codes the value is a list of possible nucleotides
iupac_conversions = {"A":"A", "T":"T", "U":"T", "C":"C", "G":"G", "R":["A","G"], "Y":["C", "T"], "S":["G", "C"], "W":["A", "T"], "K":["G", "T"], "M":["A", "C"],"B":["C", "G", "T"],"D":["A", "G", "T"], "H":["A", "C", "T"], "V":["A", "C", "G"],"N":["A", "C", "G", "T"]}

#########################################################################################################################################
#motif_maker takes the input of a .txt file that contains one motif per line (while this program will technically work with as many motifs as you want, the default color scheme allows for ten colors)
def motif_maker(motif_file):
    """This function will take a file argument and will create a dictionary of DNA motifs that we are trying to analyze, \ 
    the keys will be the original motif, while the values are all possible motifs that correspond to it.  A second dictionary will be made in case the fasta sequence is RNA"""
    #create two empty dictionaries, one for DNA and one for RNA motifs
    #the keys will be the actual motif sequence with the value of all possible motif sequences (after adjusting them for ambiguous nucleotides)    
    DNA = {}
    RNA = {}

    #open the motif file
    with open(motif_file, "r") as fh:
        #go through the file line by line
        for line in fh:
            #remove the newline character
            line = line.strip()
            #make the motif sequence uppercase
            line = line.upper()
            #create new dictionary entries for each motif
            DNA[line] = []
            RNA[line] = []

    #for every motif in the dictionary
    for motif in DNA:
        #intialize empty temporary lists
        new_dna_motif=[]
        new_rna_motif=[]
        
        #use some fancy list comprehension
        #for letter in motif, create an itertools.product permutation replacing that letter with each other possible nucleotide
        for i in itertools.product(*[iupac_conversions[j] for j in motif]):
            #by default iupac_conversions dictionary creates DNA motif sequences
            new_dna_motif.append("".join(i))
            #replace every thymine with uracil for the RNA motif sequences
            new_rna_motif.append("".join(i).replace("T", "U"))

        #add all possible motif sequences into respective dictionaries
        DNA[motif] = new_dna_motif
        RNA[motif] = new_rna_motif
    
    #return dictionaries
    return DNA, RNA
#########################################################################################################################################



#########################################################################################################################################
#function that parses through files extracting header lines and sequences
def parser(fh):
    """This function will parse through the file passed to it and return a dictionary containing the keys of sequence headers, and the value of a list containing the sequence itself.  \
    This list will contain the sequence split into three components; the region prior to the exon, the exon, and the region following the exon. Note that this only works if the exon sequences are capitalized."""
    
    #initialize a new dictionary
    sequences = {}
    #open the file passed to the method
    with open(fh, "r") as file:
        #create empty header and sequence strings
        header = ""
        sequence = ""
        #for each line in the file
        for line in file:
            #remove the newline character
            line = line.strip()
            #if line is a header line it would start with a ">"
            if line.startswith(">"):
                #if the len(header) > 0 then we know that we are not on the first entry so a sequence exists
                #the very first line would be a header line, but we do not want to split the sequence, as we have not yet reached the sequence
                if len(header) > 0:
                    #split the sequence by lowercase and uppercase letters and assign the resulting list to the sequences dictionary with the respective header line as the key
                    sequences[header] = re.split(r"([A-Z]+)", sequence)
                #assign the header line to header
                header = line
                #reset the sequence
                sequence = ""
            #else the line is not a header
            else:
                #concatenate the new line to the last line
                sequence += line 
        #split and assign the very last sequence in the file
        #this is necessary because the last sequence split will not be triggered by a new header
        sequences[header] = re.split(r"([A-Z]+)", sequence)
    #return the dictionary containing the sequences
    return sequences
#########################################################################################################################################


#########################################################################################################################################
#function that determines if the given gene sequences is DNA
def determine_type(sequence):
    """This function determines if the given gene sequence is DNA. If a T is present we know the sequence is DNA, while if a U is present we know the sequence is RNA"""
    #in case one of the segments of sequence does not contain a U for loop through each one
    for item in sequence:
        #if U is in the sequence return False
        if "U" in item.upper():
            return False
        #else if U is not in the sequence return True
        else:
            return True
#########################################################################################################################################




#########################################################################################################################################
#this function will take the sequence as well as a dictionary of motifs and return a dictionary containing the motif as a key, and a list of the starting positions as values
def determine_coords(sequence, motif):
    """This function will take the sequence dictionary entry (a list of length three) as well as the dictionary of DNA motifs ifDNA is true or RNA motifs ifDNA is false. \
    This function will return a dictionary containing the given motif sequence as a key, with a list of all the starting positions determined as values"""
    #return all coords for given motif in sequence
    #concatenate the three parts of the sequence together into one string
    sequence = sequence[0] + sequence[1] + sequence[2]
    #since the motifs are all uppercase, change the sequence to all uppercase
    sequence = sequence.upper()
    #create new, empty dictionary
    coords = {}
    
    #for motif in our dictionary
    for key in motif:
        #create a new dictionary entry for that motif
        coords[key] = []
        
        #determine the length of the motif
        motif_length = len(key)
        
        #for possible motif sequence per given motif
        for value in motif[key]:
            #for every possible sliding window in the sequence for that motif (if motif sequence is of length 4 this will take position 0-3 then 1-4, 2-5, etc.)
            for i in range(0, len(sequence)-motif_length + 1):
                #what is the sequence of the sliding window
                check = sequence[i: i + motif_length]
                #is the sliding window sequence the same sequence as the motif we are checking
                if check == value:
                    #if it is append the start position of the motif
                    coords[key].append(i)
        #sort the coordinates for ease of plotting
        coords[key].sort()
    #return the new dictionary containing the coordinates    
    return coords
#########################################################################################################################################



#########################################################################################################################################
#this function will return the required size of images based onthe height of the key (number of motifs), longest gene, and number of genes
def determine_fig_size(sequences):
    """This function takes in dictionary containing all sequences and determines which sequence is the longest to set the width of the image.  In addition it determines the number of motifs \
    and number of genes to determine the required height of the image"""
    #number of genes is simply the number of entries in the sequences dictionary
    number_of_genes = len(sequences)
    #set temporary int to 0
    max_length = 0
    #iterate through dictionary
    for item in sequences:
        #concatenate strings together to determine length of total sequence
        sequence = sequences[item][0] + sequences[item][1] + sequences[item][2]
        #if the length of the sequence is greater then the max length assign a new max lenght, else do nothing
        if len(sequence) > max_length:
            max_length = len(sequence)
        
    #add 100 to maximum sequence length for borders of image
    figure_width = max_length + 100
    #multiply number of genes by 150 for each gene + 50 for border + 30 times the number of motifs for the size of key
    figure_height = number_of_genes * 150 + 50 + 30 * len(DNA)

    #return dimensions
    return figure_width, figure_height
#########################################################################################################################################



#########################################################################################################################################

def make_legend(cr, motif):
    """This function takes the pycairo context surface as well as the dictionary of motifs and creates a key for the image"""
    #draw box for key
    cr.set_line_width(2)
    cr.rectangle(25, 25, 180, 30 * len(motif))
    cr.stroke()
    cr.set_line_width(5)
    color = 0
    motif_num = 1
    #for each motif in the dictionary
    for item in motif:
        #get color from cmap object created
        rgb = cmap(color)
        #set color
        cr.set_source_rgb(rgb[0], rgb[1], rgb[2])
        #set and move to position based on motif number
        pos = motif_num * 30 + 10
        cr.move_to(50, pos)
        cr.line_to(80, pos)
        cr.move_to(85, pos+3)
        #print line with color (above) and motif sequence (below)
        cr.show_text(item)
        cr.stroke()

        #increment color counter indentically to how it was incremented below
        if(len(DNA) < 10):
            color += 0.1
        else:
            color += 1/num_colors_needed

        #increment motif number which is used to calculate y position
        motif_num += 1

#########################################################################################################################################


#########################################################################################################################################

def plot_gene(context, sequence, gene_number, title):
    """This function takes the pycairo context surface, the sequence dictionary entry, the gene_number, and the gene title and will draw the genes rectangle outline as well as the intron \
    before and after the exon as well as the exon. The length of the three segments of the sequence correlate directly to the length of the sequence"""
    #determine lengths of each portion of sequence
    pre_exon_length = len(sequence[0])
    exon_length = len(sequence[1])
    post_exon_length = len(sequence[2])

    #set color to black
    context.set_source_rgb(0,0,0)
    context.set_line_width(2)

    #determine what y coordinates should be used for top left corner of gene box outline
    ypos=(gene_number-1) * 150 + 30*len(DNA) + 25
    
    #paste gene header in top left corner of gene box
    context.move_to(50, ypos + 25)
    context.show_text(title)

    #draw rectangle for gene box outline
    context.rectangle(25, ypos, 50 + pre_exon_length+exon_length+post_exon_length, 150)

    #add 115 to yposition to indicate where the intron should be started to drawa
    ypos = ypos + 115
    #draw first intron (pre-exon)
    context.move_to(50, ypos)
    context.line_to(50+pre_exon_length, ypos)
    #draw exon
    context.rectangle(50+pre_exon_length, ypos-20, exon_length, 40)
    #draw sepond intron (post-exon)
    context.move_to(50+pre_exon_length+exon_length, ypos)
    context.line_to(50+pre_exon_length+exon_length+post_exon_length, ypos)
    #add lines to canvas
    context.stroke()
    #return ypos
    return ypos 
#########################################################################################################################################


#########################################################################################################################################

#create DNA and RNA dictionaries using motif_maker function described above
DNA, RNA = motif_maker(motif_file)
#determine the number of motifs being used
num_colors_needed = len(DNA)

#using matplotlib color palettes get cmap color object
cmap = plt.get_cmap(name=color)


#for each fasta file run the program on it
for file in fasta_files:
    #extract file name for use in naming plot
    figname = re.findall('(.*)\.fa', file)[0] + ".svg"
    #sequences = parser(file) will set sequences equal to a dictionary containing each header and sequence in the fasta file
    sequences = parser(file)

    #use the determine_fig_size function to return the figure width and height
    figure_width, figure_height = determine_fig_size(sequences)

    #create a pycairo surface
    surface = cairo.SVGSurface(figname, figure_width, figure_height)

    #create context on the pycairo surface
    context = cairo.Context(surface)

    #set gene incrementer variable
    gene_count = 1

    #iterate through each sequence in the sequences dictionary
    for sequence in sequences:
        #plot intron and exon lengths
        ypos = plot_gene(context, sequences[sequence], gene_count, sequence)

        #determine if the sequence is DNA, if it is use the DNA motifs, else use the RNA motifs
        isDNA =  determine_type(sequences[sequence])
        if isDNA:
            coords = determine_coords(sequences[sequence], DNA)
        else:
            coords = determine_coords(sequences[sequence], RNA)

        #create an empty dictionary and tracker integer
        #used_pos keys are indexes on sequences while values will be the number of motifs that overlap
        used_pos = {}
        height = 0

        #determine what segment of color palette to use first
        color = 0

        #coords dictionary contains each motif as well as there starting positions
        #iterate through this dictionary one motif at a time
        for motif in coords:
            #extract color from color palette
            rgb = cmap(color)

            #iterate through each startpos for each motif
            for startpos in coords[motif]:
                #set color
                context.set_source_rgb(rgb[0], rgb[1], rgb[2])

                max_height = 0

                #for the numbers from 0 to length of motif
                for i in range(0, len(motif)):
                    #x would be every coordinate the motif overlaps with
                    x = i + startpos

                    #if x is in the used_pos dictionary
                    if x in used_pos:
                        #increment the height the next motif will be plotted at by one
                        height_increment = used_pos[x] + 1
                    else:
                        #no other motif is found here so set the height increment to one
                        height_increment = 1
                    
                    #if height increment is greater then max_height
                    if height_increment > max_height:
                        #max heigh equals height increment
                        max_height = height_increment
                #go back through range of motif
                for i in range(0, len(motif)):
                    #set x again
                    x = i + startpos
                    #change height in dictionary
                    used_pos[x] = max_height

                #determine height needed
                height += 6 * used_pos[startpos]
                context.set_line_width(0.1)
                #draw motif
                context.rectangle(50+startpos, ypos-2-height, len(motif), 4)
                context.fill_preserve()
                #change color to black so you can get a black outline
                context.set_source_rgb(0,0,0)
                context.stroke()
                #reset height to 0
                height = 0

            #increment color variable (ten colors in default palette. if more than ten motifs new palette must be chosen)
            if(len(DNA) < 10):
                color += 0.1
            else:
                color += 1/num_colors_needed
        #increment gene_count by one
        gene_count +=1


    #make and add legends
    make_legend(context, DNA)
    #double check to make sure everything shows up on surface
    context.stroke()
    #finish surface
    surface.finish()

#####################################################################