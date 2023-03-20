#!/usr/bin/env python
#
#.-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#/ / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
#`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
# motif-mark-OOP
#
# Bi625 2023
# sj kim
#
# script to visualize motifs and exons on sequences, using pycairo
#
# script overview: Takes a FASTA file (seqs are less than 1000 bases)
# and a motif file (less than 10 bases each, a different motif per line)
# and generates a png image of each sequence marked with exons and motifs.

#### [ port ] ####
import cairo
import math
import argparse
import re
import os

def get_args():
    parser = argparse.ArgumentParser(description= "Graphical representation of motifs, introns, and exons on sequences.")
    parser.add_argument("-f", "--fasta", help = "input fasta file", required = True)
    parser.add_argument("-m", "--motif", help = "file of valid motifs",  required = True)
    return parser.parse_args()

#############################################################
#                                                           #
#                         Classes                           #
#                                                           #
#############################################################

class Sequence:
    ''' 
    Sequence class describes a single record of a fasta file. 
    It looks for exons in a given sequence by looking for capitalized letters in the input sequence. 
    NOTE: This describes our given fasta file's convention - could be different depending on how input fasta was generated.
    Includes methods to draw the sequence and the exons.
    '''
    def __init__(self, sequence_, y_pos_, context_):
        """ defines Sequence objects with the following characteristics """
        self.sequence = sequence_
        self.seq_line_width = 5 #arbitrary width of line for sequences
        self.length = len(self.sequence)
        self.exon_seq = ""
        self.y_pos = y_pos_
        self.context = context_

    def draw_seq(self):
        """ draws a sequence as a line, using the length of the sequence """
        self.context.set_source_rgba(0, 0, 0, 1) #black
        self.context.set_line_width(self.seq_line_width)
        self.context.move_to(0, self.y_pos)
        self.context.line_to(self.length, self.y_pos)
        self.context.stroke()

    def draw_exon(self):
        """ draws an exon as a thick black line, over the sequence at the exon's location """
        self.exon_seq = re.search(r'([A-Z]+)', self.sequence)
        exon_start = self.exon_seq.span()[0]
        exon_stop = self.exon_seq.span()[1]

        self.context.set_source_rgba(0, 0, 0, 1) #black
        self.context.set_line_width(10*self.seq_line_width) # exon should be thickest
        self.context.move_to(exon_start, self.y_pos)
        self.context.line_to(exon_stop, self.y_pos)
        self.context.stroke()

class Motif(Sequence):
    ''' 
    Motif class describes individual motifs and draws motifs from input file. Different motifs are 
    represented as different colors and are drawn to scale on a sequence. 
    This motif class can handle ambiguous bases, as defined by the IUPAC degenerate base symbols.
    This class inherits features from the Sequence class.
    '''
    def __init__(self, sequence_, y_pos_, context_, color_):
        super().__init__(sequence_, y_pos_, context_)
        self.color = color_
        self.start_pos = list()
        self.translated = ""
        self.universalMotifs()

    def universalMotifs(self):
        '''
        motifs are in RNA, need to translate to DNA
        Have to account for degenerate bases:
        IUPAC degenerate base symbols
        '''
        translateMotif = self.sequence.replace("U","T")
        translateMotif = translateMotif.replace("W", "[AT]")
        translateMotif = translateMotif.replace("S", "[CG]")
        translateMotif = translateMotif.replace("M", "[AC]")
        translateMotif = translateMotif.replace("K", "[GT]")
        translateMotif = translateMotif.replace("R", "[AG]")
        translateMotif = translateMotif.replace("Y", "[CT]")
        translateMotif = translateMotif.replace("B", "[C, G, T]")
        translateMotif = translateMotif.replace("D", "[A, G, T]")
        translateMotif = translateMotif.replace("H", "[A, C, T]")
        translateMotif = translateMotif.replace("V", "[A, C, G]")
        translateMotif = translateMotif.replace("N", "[A, C, G, T]")

        # use lookahead for finding overlapping motifs
        translateMotif = "(?=(" + translateMotif + "))"
        self.translated = translateMotif

    def motifLocation(self, sequence):
        '''
        Locates all the start positions of the motif in the sequence.
        Populates list of start positions relative to the sequence.
        '''
        motif_locations = re.finditer(rf"{self.translated}", sequence.upper())
        self.start_pos = list(m.start() for m in motif_locations)

    def drawMotifs(self):
        '''
        This function draws motifs as colorful rectangles
        '''
        for start in self.start_pos:
            self.context.set_source_rgba(self.color[0], self.color[1], self.color[2], self.color[3])
            #self.context.set_source_rgba(0.83529, 0.36863, 0.00000, .8)
            self.context.rectangle(start, (self.y_pos - 50/2), self.length, 50)
            self.context.fill()
            self.context.stroke()

#############################################################
#                                                           #
#                    Global Functions                       #
#                                                           #
#############################################################

def onelineFasta(originalFaFile, intermediaryFaFile : str):
    '''This function takes a .fasta file and removes newline breaks in the protein sequence and
    creates an intermediary .fasta file with 2 lines per record: the header and the sequence'''
    with open(originalFaFile, "r") as ogInfile:
        with open(intermediaryFaFile, "w") as onelineOutfile:
            for i, line in enumerate(ogInfile):
                header = ""
                record = ""
                if line.startswith(">"):
                    header = line
                    if i == 0:
                        onelineOutfile.write(f"{header}")
                    else:
                        onelineOutfile.write(f"\n{header}")
                else:
                    line = line.strip("\n")
                    record += line
                    onelineOutfile.write(f"{record}")

def importMotifs(motifInfile:str):
    ''' 
    Takes motif file with one motif per line
    returns list of motifs, all upper-case
    '''
    motifs = []
    with open(motifInfile, "r") as motifFile:
        for line in motifFile:
            line = line.strip()
            motifs.append(line.upper())
    return motifs

def fastaDimensions(oneline_fasta_file):
    '''
    This function takes transformed one-line fasta file and returns the 
    number of records, the length of the longest sequence and the dictionary
    that holds header:sequence from the fasta file
    The dimensions are mostly to set the canvas size but also for printing
    the header.
    '''
    recordCount = 0
    longestSeq = 0
    faDictionary = {}
    header = ""
    record = ""
    with open(oneline_fasta_file, "r") as faOnelineInfile:
        for line in faOnelineInfile:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else: #line doesn't start with >, indicating it's a record
                record = line
                if len(record) > longestSeq:
                    longestSeq = len(record)
            faDictionary[header] = record
            recordCount = len(faDictionary)
    return recordCount, longestSeq, faDictionary

#############################################################
#                                                           #
#                          Main                             #
#                                                           #
#############################################################

def main():
    # first, process the input files through argparse
    args = get_args()
    fasta_file = args.fasta
    motif_file = args.motif
    png_file = re.sub(r'(.fasta$)', ".png", fasta_file)
    fasta_oneline = re.sub(r'(.fasta$)', "_oneline.fasta", fasta_file)

    # process intake files
    onelineFasta(fasta_file, fasta_oneline)
    motifs = importMotifs(motif_file)

    # rgbaMotifs is a list of global colors for motifs
    # transparency applied to rgba colors here so overlaps can show through
    rgbaMotifs = [(0.83529, 0.36863, 0.00000, .8), # vermillion
                (0.90196, 0.62353, 0.00000, .8), # orange
                (0.33725, 0.70588, 0.91373, .8), # sky blue
                (0.00000, 0.61961, 0.45098, .8), # green
                (0.94118, 0.89412, 0.25882, .8), # yellow
                (0.00000, 0.44706, 0.69804, .8), # blue
                (1, 0, 0, .8), # red
                (0.80000, 0.47451, 0.65490, .8)] # pink

    # get dimensions for canvas size
    numRecords, longestRecord, fastaDictionary = fastaDimensions(fasta_oneline)

    # set canvas sizing
    WIDTH, HEIGHT = (100+longestRecord, (100*numRecords)+200) #size of canvas, +400 to leave room for legend
    recordHeight = 50  # arbitrary height for exons and motifs
    start_y = range(0, 100*numRecords, 100) # range from 0 to 100*numRecords, increment by 100

    # make canvas
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    context = cairo.Context(surface)
    context.set_source_rgba(1, 1, 1, 1) # white background
    context.rectangle(0, 0, WIDTH, HEIGHT)  
    context.fill()
    context.stroke()
    context.translate(50, 50)    # shift by 50 to the right for all following code, to leave 50 for a margin

    # make legend outline
    legendWidth = 175 # arbitrary size
    legendHeight = 16*(len(motifs)+2) #leave room for all motifs plus an exon and intron entry

    context.move_to(0, 0)
    context.set_source_rgba(0, 0, 0, 1)
    context.set_font_size(20)
    context.show_text("Legend")
    context.set_line_width(1)
    context.rectangle(0, 5, legendWidth, legendHeight)
    context.stroke()

    # fill in legend
    ## start with motifs
    for i in range(len(motifs)):
        # color labels
        context.set_source_rgba(rgbaMotifs[i][0], rgbaMotifs[i][1], rgbaMotifs[i][2], 1)
        context.rectangle(5, 10 + (i*15), 10, 10)
        context.fill()
        context.stroke()
        # motifs 
        context.move_to(20, 20 + i*15)
        context.set_source_rgba(0, 0, 0, 1)
        context.set_font_size(15)
        context.show_text(f"Motif {motifs[i]}")

    ## add exon to legend
    context.set_source_rgba(0, 0, 0, 1)
    context.rectangle(5, 10 + len(motifs)*15, 10, 10)
    context.fill()
    context.stroke()
    context.move_to(20, 20 + len(motifs)*15)
    context.show_text(f"Exon") 

    ## add intron to legend
    context.move_to(5, 13 + (len(motifs)+1)*15)
    context.set_line_width(2)
    context.line_to(15, 13 + (len(motifs)+1)*15)
    context.stroke()
    context.set_line_width(1)
    context.move_to(20, 20 + (len(motifs)+1)*15)
    context.show_text("Intron")

    context.translate(0,150) # move context away from legend

    # draw the rest: sequences, exons, and motifs 
    # parse through previously made fastaDictionary from dimensions function
    recordNum = 0
    for entry in fastaDictionary:
        seq = Sequence(fastaDictionary[entry], start_y[recordNum], context)
        seq.draw_seq()
        seq.draw_exon()
        for m in range(len(motifs)):
            motifObj = Motif(motifs[m], start_y[recordNum], context, rgbaMotifs[m])
            motifObj.motifLocation(fastaDictionary[entry])
            motifObj.drawMotifs()
        # draw headers after drawing the sequences
        context.move_to(0, start_y[recordNum]-30) # give each sequence a header of -30 pixels
        context.set_source_rgba(0, 0, 0, 1)
        context.show_text(f"{entry}")
        recordNum += 1

    surface.write_to_png(png_file)
    surface.finish()
    
    # remove temp file
    if os.path.exists(fasta_oneline):
        os.remove(fasta_oneline)

if __name__ == "__main__":
    main()
