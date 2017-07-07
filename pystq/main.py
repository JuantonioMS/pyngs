 #!/usr/bin/python
import argparse
import time
import subprocess
import shlex
import sys
import collections
from zipfile import ZipFile
import numpy as np
from itertools import islice
from pystq.plots import BasicPlots, SamPlots
from pystq import SamReader, FastqReader, FastQ, Sam, Fasta, FastaReader

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

def kmer_detection(seq, n):
    iterator = iter(seq)
    result = tuple(islice(iterator, n))
    if len(result) == n:
        yield ''.join(result)
    for element in iterator:
        result = result[1:] + (element,)
        yield ''.join(result)

def debug(VERBOSE_LEVEL, *args, **kwargs):

    if VERBOSE_LEVEL == DEBUG_LEVEL:
        print DEBUG_LEVEL + "! " + "".join(str(x) for x in args)

def run(arguments):

    start_time = time.time()


    arguments['input'] = argparse.FileType('r')(arguments['input'])
    arguments['text'] = argparse.FileType('w')(arguments['text'])

    args = Bunch(arguments)
    debug("Info", "Starting the program")
    debug("Flag", "Start time: ", start_time)



    name, ext = args.input.name.split(".")
    sys.stdout.write("File %s loaded, extension .%s detected\n" % (args.input.name, ext))

    order = "wc -l " + args.input.name
    comamand_line = shlex.split(order)

    try:
        total_lines, nothing = subprocess.check_output(comamand_line).split(" ")
        total_lines = int(total_lines)
        sys.stdout.write("Detected %i lines\n" % total_lines)
    except:
        pass

    if ext not in ("fq", "fastq", "fasta", "sam"):
        sys.exit("Wrong extension!\n" + "Please, check your file extension is (fq, fastq, fasta, sam).")

    if ext in ("fq", "fastq"):
        debug("Flag", "FastQ section enter")

        try:
            estimate_sequences = total_lines / 4
            sys.stdout.write("Estimated %i sequences\n" % estimate_sequences)

        except:
            with FastqReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1


        store_nuc = collections.defaultdict(lambda: collections.defaultdict(int))
        store_qual = collections.defaultdict(lambda: collections.defaultdict(int))
        store_kmers = collections.defaultdict(lambda: collections.defaultdict(int))
        store_duplicants = collections.defaultdict(int)
        store_singlekmer = collections.defaultdict(int)
        store_len = collections.defaultdict(int)
        store_gc = collections.defaultdict(int)
        store_qualitysec = collections.defaultdict(int)
        store_secn = collections.defaultdict(int)


        infile = FastqReader(args.input)
        reads = infile.sampling()

        skip = False
        for read in reads:
            if not skip:
                phred_converter = read.asciidetection(read)
            skip = True

            #Store sequences lengths
            store_len[len(read.sequence)] += 1
            #Store sequences mean quality
            store_qualitysec[read.qualitysec(phred_converter)] += 1
            #Store sequence GC proportion
            store_gc[read.gc_proportion()] += 1
            #Store sequence N content
            store_secn[read.sequence.count("N")] += 1
            #Store sequences duplicated
            store_duplicants[read.sequence] += 1

            #Store nucleotides and quality per cycle
            for position, (seq, qual) in enumerate(zip(read.sequence, read.quality)):
                store_nuc[position][seq] += 1
                store_qual[position][ord(qual) - phred_converter] += 1

            #Store overexpressed kmers and cycle
            for cycle, kmer in enumerate(kmer_detection(read.sequence, args.kmerlen)):
                store_kmers[cycle][kmer] += 1
                store_singlekmer[kmer] += 1

        #Creating secn array
        array_secn = [0] * (max(store_secn) + 1)
        for element in store_secn:
            array_secn[element] = store_secn[element]
        array_secn = array_secn[1:]

        #Choosing overexpressed kmers
        sorted_store_singlekmer = sorted(store_singlekmer.items(), key=lambda(k,v): v)
        selected_kmer = sorted_store_singlekmer[-args.reprkmer:]
        selected_kmer.reverse()

        #Representing start cycle kmers
        kmer_array = np.zeros((max(store_len), len(selected_kmer)), dtype=int)
        for row, (_,cycle) in enumerate(store_kmers.items()):
            for column, i in enumerate(selected_kmer):
                key = i[0]
                kmer_array[row][column] = cycle[key]

        #Choosing duplicants
        sorted_duplicants = sorted(store_duplicants.items(), key=lambda(k,v): v)
        array_duplicants = [0] * args.reprduplicate
        for sequence, repeats in sorted_duplicants:
            try:
                array_duplicants[repeats] += 1
            except IndexError:
                continue


        #Building uality boxplot array
        qual_array = []
        for _, cycle in store_qual.items():
            column = []
            for value, times in cycle.items():
                for i in range(times):
                    column.append(value)
            qual_array.append(column)

        #Building nuc_array, gc_arrays and cyclen_array
        nuc_array = np.zeros((max(store_len.keys()), 5), dtype=float)
        array_cyclen = []
        array_gccycle = []

        for row, (_, cycle) in enumerate(store_nuc.items()):

            total = sum(cycle.values())
            nuc_array[row][0] = float(cycle["A"]) / total * 100
            nuc_array[row][1] = float(cycle["T"]) / total * 100
            C_proportion = float(cycle["C"]) / total * 100
            nuc_array[row][2] = C_proportion
            G_proportion = float(cycle["G"]) / total * 100
            nuc_array[row][3] = G_proportion
            nuc_array[row][4] = float(cycle["N"]) / total * 100
            array_cyclen.append(cycle["N"])
            array_gccycle.append(C_proportion + G_proportion)

        with ZipFile(args.zipname + '.zip', mode='w') as zip_archive:

            plotter = BasicPlots(args.input.name)
            plotter.uninuc_plot(nuc_array, zip_archive, args.grid, args.output)
            plotter.uniqual_boxplot(qual_array, zip_archive, args.grid, args.color, args.output)
            plotter.unigc_plot(array_gccycle, zip_archive, args.grid, args.output)
            plotter.uniqualitysec_bar(store_qualitysec, zip_archive, args.grid, args.color, args.output)
            plotter.unilenght_bar(store_len, zip_archive, args.grid, args.color, args.output)
            plotter.unigcproportion_scatter(store_gc, zip_archive, args.grid, args.color, args.output)
            plotter.unioverkmer_bar(selected_kmer, zip_archive, args.grid, args.color, args.output)
            plotter.unikmer_plot(kmer_array, selected_kmer, zip_archive, args.grid, args.output)
            plotter.uniduplicants_hist(array_duplicants, zip_archive, args.grid, args.color, args.output)
            plotter.unisecn_bar(array_secn, zip_archive, args.grid, args.color, args.output)
            plotter.unicyclen_bar(array_cyclen, zip_archive, args.grid, args.color, args.output)

    if ext in ("fa", "fasta"):
        debug("Flag", "Fasta section enter")

        try:
            sys.stdout.write("Estimated %i sequences\n" % 10)

        except:
            with FastaReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1


        store_nuc = collections.defaultdict(lambda: collections.defaultdict(int))
        store_kmers = collections.defaultdict(lambda: collections.defaultdict(int))
        store_duplicants = collections.defaultdict(int)
        store_singlekmer = collections.defaultdict(int)
        store_len = collections.defaultdict(int)
        store_gc = collections.defaultdict(int)
        store_secn = collections.defaultdict(int)


        infile = FastaReader(args.input)
        reads = infile.sampling()

        for read in reads:

            #Store sequences lengths
            store_len[len(read.sequence)] += 1
            #Store sequence GC proportion
            store_gc[read.gc_proportion()] += 1
            #Store sequence N content
            store_secn[read.sequence.count("N")] += 1
            #Store sequences duplicated
            store_duplicants[read.sequence] += 1

            #Store nucleotides and quality per cycle
            for position, seq in enumerate(read.sequence):
                store_nuc[position][seq] += 1

            #Store overexpressed kmers and cycle
            for cycle, kmer in enumerate(kmer_detection(read.sequence, args.kmerlen)):
                store_kmers[cycle][kmer] += 1
                store_singlekmer[kmer] += 1

        #Creating secn array
        array_secn = [0] * (max(store_secn) + 1)
        for element in store_secn:
            array_secn[element] = store_secn[element]
        array_secn = array_secn[1:]

        #Choosing overexpressed kmers
        sorted_store_singlekmer = sorted(store_singlekmer.items(), key=lambda(k,v): v)
        selected_kmer = sorted_store_singlekmer[-args.reprkmer:]
        selected_kmer.reverse()

        #Representing start cycle kmers
        kmer_array = np.zeros((max(store_len), len(selected_kmer)), dtype=int)
        for row, (_,cycle) in enumerate(store_kmers.items()):
            for column, i in enumerate(selected_kmer):
                key = i[0]
                kmer_array[row][column] = cycle[key]

        #Choosing duplicants
        sorted_duplicants = sorted(store_duplicants.items(), key=lambda(k,v): v)
        array_duplicants = [0] * args.reprduplicate
        for sequence, repeats in sorted_duplicants:
            try:
                array_duplicants[repeats] += 1
            except IndexError:
                continue

        #Building nuc_array, gc_arrays and cyclen_array
        nuc_array = np.zeros((max(store_len.keys()), 5), dtype=float)
        array_cyclen = []
        array_gccycle = []

        for row, (_, cycle) in enumerate(store_nuc.items()):

            total = sum(cycle.values())
            nuc_array[row][0] = float(cycle["A"]) / total * 100
            nuc_array[row][1] = float(cycle["T"]) / total * 100
            C_proportion = float(cycle["C"]) / total * 100
            nuc_array[row][2] = C_proportion
            G_proportion = float(cycle["G"]) / total * 100
            nuc_array[row][3] = G_proportion
            nuc_array[row][4] = float(cycle["N"]) / total * 100
            array_cyclen.append(cycle["N"])
            array_gccycle.append(C_proportion + G_proportion)

        with ZipFile(args.zipname + '.zip', mode='w') as zip_archive:

            plotter = BasicPlots(args.input.name)
            plotter.uninuc_plot(nuc_array, zip_archive, args.grid, args.output)
            plotter.unigc_plot(array_gccycle, zip_archive, args.grid, args.output)
            plotter.unilenght_bar(store_len, zip_archive, args.grid, args.color, args.output)
            plotter.unigcproportion_scatter(store_gc, zip_archive, args.grid, args.color, args.output)
            plotter.unioverkmer_bar(selected_kmer, zip_archive, args.grid, args.color, args.output)
            plotter.unikmer_plot(kmer_array, selected_kmer, zip_archive, args.grid, args.output)
            plotter.uniduplicants_hist(array_duplicants, zip_archive, args.grid, args.color, args.output)
            plotter.unisecn_bar(array_secn, zip_archive, args.grid, args.color, args.output)
            plotter.unicyclen_bar(array_cyclen, zip_archive, args.grid, args.color, args.output)

    elif ext in ["sam"]:

        try:
            estimate_sequences = total_lines
            sys.stdout.write("Estimated %i sequences\n" % estimate_sequences)

        except:
            with SamReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1

        store_flag = [0] * 12
        store_cigar = [0, "M", 0, "I", 0, "D", 0, "N", 0, "S", 0, "H", 0, "P", 0, "=", 0, "X"]
        store_nuc = collections.defaultdict(lambda: collections.defaultdict(int))
        store_qual = collections.defaultdict(lambda: collections.defaultdict(int))
        store_kmers = collections.defaultdict(lambda: collections.defaultdict(int))
        store_duplicants = collections.defaultdict(int)
        store_singlekmer = collections.defaultdict(int)
        store_len = collections.defaultdict(int)
        store_gc = collections.defaultdict(int)
        store_qualitysec = collections.defaultdict(int)
        store_secn = collections.defaultdict(int)
        store_position = collections.defaultdict(int)
        store_mapq = collections.defaultdict(int)


        infile = SamReader(args.input)
        reads = infile.sampling()

        skip = False
        for read in reads:
            if not skip:
                phred_converter = read.asciidetection(read)
            skip = True

            #Store sequences lengths
            store_len[len(read.sequence)] += 1
            #Store sequences mean quality
            store_qualitysec[read.qualitysec(phred_converter)] += 1
            #Store sequence GC proportion
            store_gc[read.gc_proportion()] += 1
            #Store sequence N content
            store_secn[read.sequence.count("N")] += 1
            #Store sequences duplicated
            store_duplicants[read.sequence] += 1
            #Store mapping position
            store_position[read.position] += 1
            #Store Flags
            store_flag = read.flaglist(store_flag)
            #Store mapping quality
            store_mapq[read.mapq] += 1
            #Store cigar
            store_cigar = read.cigarlist(store_cigar)

            #Store nucleotides and quality per cycle
            for position, (seq, qual) in enumerate(zip(read.sequence, read.quality)):
                store_nuc[position][seq] += 1
                store_qual[position][ord(qual) - phred_converter] += 1

            #Store overexpressed kmers and cycle
            for cycle, kmer in enumerate(kmer_detection(read.sequence, args.kmerlen)):
                store_kmers[cycle][kmer] += 1
                store_singlekmer[kmer] += 1

        #Creating cigar array
        array_cigar = []
        for n, number in enumerate(store_cigar):
            if n % 2 == 0:
                array_cigar.append(number)
            elif n % 2 == 1:
                continue

        #Creating secn array
        array_secn = [0] * (max(store_secn) + 1)
        for element in store_secn:
            array_secn[element] = store_secn[element]
        array_secn = array_secn[1:]

        #Creating position array
        array_position = [0] * (max(store_position) + 1)
        for element in store_position:
            array_position[element] = store_position[element]
        array_position = array_position[1:]

        #Creating mapq array
        array_mapq = [0] * (max(store_mapq) + 1)
        for element in store_mapq:
            array_mapq[element] = store_mapq[element]


        #Choosing overexpressed kmers
        sorted_store_singlekmer = sorted(store_singlekmer.items(), key=lambda(k,v): v)
        selected_kmer = sorted_store_singlekmer[-args.reprkmer:]
        selected_kmer.reverse()

        #Representing start cycle kmers
        kmer_array = np.zeros((max(store_len), len(selected_kmer)), dtype=int)
        for row, (_,cycle) in enumerate(store_kmers.items()):
            for column, i in enumerate(selected_kmer):
                key = i[0]
                kmer_array[row][column] = cycle[key]

        #Choosing duplicants
        sorted_duplicants = sorted(store_duplicants.items(), key=lambda(k,v): v)
        array_duplicants = [0] * args.reprduplicate
        for sequence, repeats in sorted_duplicants:
            try:
                array_duplicants[repeats] += 1
            except IndexError:
                continue

        #Building uality boxplot array
        qual_array = []
        for _, cycle in store_qual.items():
            column = []
            for value, times in cycle.items():
                for i in range(times):
                    column.append(value)
            qual_array.append(column)

        #Building nuc_array, gc_arrays and cyclen_array
        nuc_array = np.zeros((max(store_len.keys()), 5), dtype=float)
        array_cyclen = []
        array_gccycle = []

        for row, (_, cycle) in enumerate(store_nuc.items()):

            total = sum(cycle.values())
            nuc_array[row][0] = float(cycle["A"]) / total * 100
            nuc_array[row][1] = float(cycle["T"]) / total * 100
            C_proportion = float(cycle["C"]) / total * 100
            nuc_array[row][2] = C_proportion
            G_proportion = float(cycle["G"]) / total * 100
            nuc_array[row][3] = G_proportion
            nuc_array[row][4] = float(cycle["N"]) / total * 100
            array_cyclen.append(cycle["N"])
            array_gccycle.append(C_proportion + G_proportion)

        with ZipFile(args.zipname + '.zip', mode='w') as zip_archive:

            plotter = SamPlots(args.input.name)
            plotter.uninuc_plot(nuc_array, zip_archive, args.grid, args.output)
            plotter.uniqual_boxplot(qual_array, zip_archive, args.grid, args.color, args.output)
            plotter.unigc_plot(array_gccycle, zip_archive, args.grid, args.output)
            plotter.uniqualitysec_bar(store_qualitysec, zip_archive, args.grid, args.color, args.output)
            plotter.unilenght_bar(store_len, zip_archive, args.grid, args.color, args.output)
            plotter.unigcproportion_scatter(store_gc, zip_archive, args.grid, args.color, args.output)
            plotter.unioverkmer_bar(selected_kmer, zip_archive, args.grid, args.color, args.output)
            plotter.unikmer_plot(kmer_array, selected_kmer, zip_archive, args.grid, args.output)
            plotter.uniduplicants_hist(array_duplicants, zip_archive, args.grid, args.color, args.output)
            plotter.unisecn_bar(array_secn, zip_archive, args.grid, args.color, args.output)
            plotter.unicyclen_bar(array_cyclen, zip_archive, args.grid, args.color, args.output)
            plotter.uniflag_bar(store_flag, zip_archive, args.grid, args.color, args.output)
            plotter.uniposition_bar(array_position, zip_archive, args.grid, args.color, args.output)
            plotter.unimapq_bar(array_mapq, zip_archive, args.grid, args.color, args.output)
            plotter.unicigar_bar(array_cigar, zip_archive, args.grid, args.color, args.output)

    elif ext in ["fasta", "mpfa", "fna", "fsa", "fas"]:
        debug("Flag", "Fasta section enter")

        try:
            estimate_sequences = total_lines / 2
            sys.stdout.write("Estimated %i sequences\n" % estimate_sequences)

        except:
            with FastaReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1


        store_nuc = collections.defaultdict(lambda: collections.defaultdict(int))
        store_kmers = collections.defaultdict(lambda: collections.defaultdict(int))
        store_duplicants = collections.defaultdict(int)
        store_singlekmer = collections.defaultdict(int)
        store_len = collections.defaultdict(int)
        store_gc = collections.defaultdict(int)
        store_secn = collections.defaultdict(int)

        infile = FastaReader(args.input)
        reads = infile.sampling()

        skip = False
        for read in reads:
            if not skip:
                phred_converter = read.asciidetection(read)
            skip = True

            #Store sequences lengths
            store_len[len(read.sequence)] += 1
            #Store sequences mean quality
            store_qualitysec[read.qualitysec(phred_converter)] += 1
            #Store sequence GC proportion
            store_gc[read.gc_proportion()] += 1
            #Store sequence N content
            store_secn[read.sequence.count("N")] += 1
            #Store sequences duplicated
            store_duplicants[read.sequence] += 1

            #Store nucleotides and quality per cycle
            for position, (seq, qual) in enumerate(zip(read.sequence, read.quality)):
                store_nuc[position][seq] += 1
                store_qual[position][ord(qual) - phred_converter] += 1

            #Store overexpressed kmers and cycle
            for cycle, kmer in enumerate(kmer_detection(read.sequence, args.kmerlen)):
                store_kmers[cycle][kmer] += 1
                store_singlekmer[kmer] += 1

        #Creating secn array
        array_secn = [0] * (max(store_secn) + 1)
        for element in store_secn:
            array_secn[element] = store_secn[element]
        array_secn = array_secn[1:]

        #Choosing overexpressed kmers
        sorted_store_singlekmer = sorted(store_singlekmer.items(), key=lambda(k,v): v)
        selected_kmer = sorted_store_singlekmer[-args.reprkmer:]
        selected_kmer.reverse()

        #Representing start cycle kmers
        kmer_array = np.zeros((max(store_len), len(selected_kmer)), dtype=int)
        for row, (_,cycle) in enumerate(store_kmers.items()):
            for column, i in enumerate(selected_kmer):
                key = i[0]
                kmer_array[row][column] = cycle[key]

        #Choosing duplicants
        sorted_duplicants = sorted(store_duplicants.items(), key=lambda(k,v): v)
        array_duplicants = [0] * args.reprduplicate
        for sequence, repeats in sorted_duplicants:
            try:
                array_duplicants[repeats] += 1
            except IndexError:
                continue


        #Building uality boxplot array
        qual_array = []
        for _, cycle in store_qual.items():
            column = []
            for value, times in cycle.items():
                for i in range(times):
                    column.append(value)
            qual_array.append(column)

        #Building nuc_array, gc_arrays and cyclen_array
        nuc_array = np.zeros((max(store_len.keys()), 5), dtype=float)
        array_cyclen = []
        array_gccycle = []

        for row, (_, cycle) in enumerate(store_nuc.items()):

            total = sum(cycle.values())
            nuc_array[row][0] = float(cycle["A"]) / total * 100
            nuc_array[row][1] = float(cycle["T"]) / total * 100
            C_proportion = float(cycle["C"]) / total * 100
            nuc_array[row][2] = C_proportion
            G_proportion = float(cycle["G"]) / total * 100
            nuc_array[row][3] = G_proportion
            nuc_array[row][4] = float(cycle["N"]) / total * 100
            array_cyclen.append(cycle["N"])
            array_gccycle.append(C_proportion + G_proportion)

        with ZipFile(args.zipname + '.zip', mode='w') as zip_archive:

            plotter = BasicPlots(args.input.name)
            plotter.uniqual_boxplot(qual_array, zip_archive, args.grid, args.color, args.output)
            plotter.uniqualitysec_bar(store_qualitysec, zip_archive, args.grid, args.color, args.output)
            plotter.unilenght_bar(store_len, zip_archive, args.grid, args.color, args.output)
            plotter.unigcproportion_scatter(store_gc, zip_archive, args.grid, args.color, args.output)
            plotter.unioverkmer_bar(selected_kmer, zip_archive, args.grid, args.color, args.output)
            plotter.uniduplicants_hist(array_duplicants, zip_archive, args.grid, args.color, args.output)
            plotter.unisecn_bar(array_secn, zip_archive, args.grid, args.color, args.output)
            plotter.unicyclen_bar(array_cyclen, zip_archive, args.grid, args.color, args.output)









def main():
    parser = argparse.ArgumentParser(prog='analyzer', description="Bioinformatics files analyzer")
    parser.add_argument('input', type=str, help="input file (Fastq, Fasta, Sam")
    parser.add_argument('-v', '--verbose', type=str, dest='DEBUG_LEVEL', action='store', default=None, help="Level screen output (default: %(default)s)")

    parser.add_argument('-z', '--zipname', type=str, default='pystq-figures', help="name of the zip (default: %(default)s)")
    parser.add_argument('-k', '--kmerlen', type=int, default=5, choices=range(2, 10), help='length of kmer for over-repesented kmer counts (default: %(default)s)')
    parser.add_argument('-K', '--reprkmer', type=int, default=7, choices=range(2, 10), help='number of represented kmers (default: %(default)s)')
    parser.add_argument('-d', '--reprduplicate', type=int, default=20, choices=range(2, 80), help='number of represented duplicated sequences (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str, default="", help="name in figures (default: %(default)s)")
    parser.add_argument('-e', '--text', type=str, default='-', help="file name for text output (default: %(default)s)")
    parser.add_argument('-g', '--grid', action="store_true", default=False, help="show grid on graphics (default: %(default)s)")
    parser.add_argument('-c', '--color', type=str, default="cornflowerblue", help="graphics color, boxplot only blue (default: %(default)s)")


    args = parser.parse_args()
    globals().update(args.__dict__)
    arguments = vars(args)
    run(arguments)

if __name__ == "__main__":
    main()
