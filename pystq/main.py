#!/usr/bin/python
import argparse
import time
import subprocess
import shlex
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from itertools import islice
from pystq.plots import BasicPlots


class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)


class FastQ(object):

    def __init__(self, name, sequence, extra, quality):
        self.name = name
        self.sequence = sequence
        self.extra = extra
        self.quality = quality

    def __str__(self):
        return "\n".join([self.name, self.sequence, self.extra, self.quality]) + "\n"

    def __len__(self):
        return len(self.sequence)

    def asciidetection(self, phredconverter = None):
        low_list = [chr(character) for character in range(33,59)]
        high_list = [chr(character) for character in range(74,127)]
        phred_value = int()
        for ascii_letter in self.quality:

            if ascii_letter in low_list:
                phred_value = 33

            elif ascii_letter in high_list:
                phred_value = 64

            if phred_value:
                return phred_value
                break

    def qualitysec(self, phred_value):
        qualities = [(ord(character) - phred_value) for character in self.quality]
        return sum(qualities) / (len(qualities))


    def gc_proportion(self):
        g_number = self.sequence.count("G")
        c_number = self.sequence.count("C")
        proportion = float(g_number + c_number) / len(self.sequence)
        return int(proportion * 100)



class FastqReader:

    def __init__(self, filename):
        self.filename = filename

    def __iter__(self):
        return self

    def ___next__(self):
        return self.next()

    def next(self):
        try:
            name = next(self.filename).strip()
            sequence = next(self.filename).strip()
            extra = next(self.filename).strip()
            quality = next(self.filename).strip()

            return FastQ(name, sequence, extra, quality)

        except StopIteration:
            raise StopIteration

    def sampling(self):
        n = 4
        for i, line in enumerate(self.filename):

            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                yield FastQ(name=name, sequence=seq, extra=strand, quality=qual)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.filename.close()

class Sam(object):

    def __init__(self, line):

        self.qname = str(line[0])
        self.flag = int(line[1])
        self.rname = str(line[2])
        self.position = int(line[3])
        self.mapq = int(line[4])
        self.cigar = str(line[5])
        self.rnext = str(line[6])
        self.pnext = int(line[7])
        self.length = int(line[8])
        self.sequence = str(line[9])
        self.quality = str(line[10])

    def __str__(self):
        return "\t".join([self.qname, str(self.flag), self.rname, str(self.position), str(self.mapq), str(self.cigar), self.rnext, str(self.pnext), str(self.length), self.sequence, self.quality]) + "\n"

    def flaglist(self, store):
        str_flag = bin(self.flag)
        str_flag = str_flag[2:]
        list_flag = [int(b) for b in str_flag]

        for position, value in enumerate(list_flag):
            store[position] += value

        return store

class SamReader():

    def __init__(self, filename):
        print "Iniciando"
        self.file = filename

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            line = next(self.filename),strip()
        except StopIteration:
            raise StopIteration

    def sampling(self):
        for line in self.file:
            if line[0] == "@":
                pass
            else:
                yield Sam(tuple(line.rstrip().split("\t")))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.filename.close()








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
    debug("Flags", "Start time: ", start_time)

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
        debug("Flags", "FastQ section enter")

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

            store_duplicants[read.sequence] += 1
            store_qualitysec[read.qualitysec(phred_converter)] += 1
            store_len[len(read.sequence)] += 1
            store_gc[read.gc_proportion()] += 1
            store_secn[read.sequence.count("N")] += 1

            for position, (seq, qual) in enumerate(zip(read.sequence, read.quality)):
                store_nuc[position][seq] += 1
                store_qual[position][ord(qual) - phred_converter] += 1

            for cycle, kmer in enumerate(kmer_detection(read.sequence, args.kmer)):
                store_kmers[cycle][kmer] += 1
                store_singlekmer[kmer] += 1


        array_secn = []
        for value in store_secn.values():
            array_secn.append(value)

        array_secn = array_secn[1:]



        sorted_store_singlekmer = sorted(store_singlekmer.items(), key=lambda(k,v): v)
        #definir un argumento por aqui
        selected_kmer = sorted_store_singlekmer[-10:]
        selected_kmer.reverse()

        sorted_duplicants = sorted(store_duplicants.items(), key=lambda(k,v): v)
        sorted_duplicants.reverse()
        """Montaje del array para la calidad (BOXPLOT)"""

        qual_array = []
        for _, cycle in store_qual.items():
            column = []
            for value, times in cycle.items():
                for i in range(times):
                    column.append(value)
            qual_array.append(column)

        """MONTAJE DEL ARRAY PARA PORCENTAJE DE NUC (PLOT)"""


        nuc_array = np.zeros((max(store_len.keys()), 4), dtype=float)
        array_gccycle = []
        cyclen_array = np.zeros((max(store_len), 1), dtype=int)
        for row, (_, cycle) in enumerate(store_nuc.items()):

            total = sum(cycle.values())
            nuc_array[row][0] = float(cycle["A"]) / total * 100
            nuc_array[row][1] = float(cycle["T"]) / total * 100
            C_proportion = float(cycle["C"]) / total * 100
            nuc_array[row][2] = C_proportion
            G_proportion = float(cycle["G"]) / total * 100
            nuc_array[row][3] = G_proportion
            cyclen_array[row] = cycle["N"]
            array_gccycle.append(C_proportion + G_proportion)

        """MONTAJE DEL ARRAY PARA PORCENTAJE DE NUC"""
        kmer_array = np.zeros((max(store_len), len(selected_kmer)), dtype=int)

        for row, (_,cycle) in enumerate(store_kmers.items()):
            for column, i in enumerate(selected_kmer):
                key = i[0]
                kmer_array[row][column] = cycle[key]

        plotter = BasicPlots(args.input.name)
        plotter.uniqual_boxplot(qual_array)
        plotter.uninuc_plot(nuc_array)
        plotter.unigc_plot(array_gccycle)
        plotter.uniqualitysec_bar(store_qualitysec)
        plotter.unilenght_bar(store_len)
        plotter.unigcproportion_scatter(store_gc)
        plotter.unioverkmer_bar(selected_kmer)
        plotter.unikmer_plot(kmer_array, selected_kmer)
        plotter.uniduplicants_hist(store_duplicants)
        plotter.unisecn_plot(array_secn)
        plotter.unicyclen_plot(cyclen_array)

    elif ext in ["sam"]:
        print "Extension SAM"

        infile = SamReader(args.input)
        reads = infile.sampling()

        flag_store = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        for read in reads:
             flag_store = read.flaglist(flag_store)


        plotsinstance = Plots()
        plotsinstance.uniflag_bar(flag_store)










def main():
    parser = argparse.ArgumentParser(prog='analyzer', description="Bioinformatics files analyzer")
    parser.add_argument('input', type=str, help="input file (Fastq, Fasta, Sam")
    parser.add_argument('-v', '--verbose', type=str, dest='DEBUG_LEVEL', action='store', default=None, help="Level screen output (default: %(default)s)")
    parser.add_argument('-a', '--name', type=str, help='sample name identifier for text and graphics output (default: input file name)')
    parser.add_argument('-n', '--nreads', type=int, default=2000000, help='number of reads sample from input (default: %(default)s)')
    parser.add_argument('-p', '--base-probs', type=str, default='0.25,0.25,0.25,0.25,0.1', help='probabilites for observing A,T,C,G,N in reads (default: %(default)s)')
    parser.add_argument('-k', '--kmer', type=int, default=5, choices=range(2, 8), help='length of kmer for over-repesented kmer counts (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str, default='fastqp_figures', help="base name for output figures (default: %(default)s)")
    parser.add_argument('-e', '--text', type=str, default='-', help="file name for text output (default: %(default)s)")
    parser.add_argument('-ll', '--leftlimit', type=int, default=1, help="leftmost cycle limit (default: %(default)s)")
    parser.add_argument('-rl', '--rightlimit', type=int, default=-1, help="rightmost cycle limit (-1 for none) (default: %(default)s)")
    parser.add_argument('-mq', '--median-qual', type=int, default=30, help="median quality threshold for failing QC (default: %(default)s)")

    align_group = parser.add_mutually_exclusive_group()
    align_group.add_argument('--aligned-only', action="store_true", default=False, help="only aligned reads (default: %(default)s)")
    align_group.add_argument('--unaligned-only', action="store_true", default=False, help="only unaligned reads (default: %(default)s)")
    parser.add_argument('-d', '--count-duplicates', action="store_true", default=False, help="calculate sequence duplication rate (default: %(default)s)")

    args = parser.parse_args()
    globals().update(args.__dict__)
    arguments = vars(args)
    run(arguments)

if __name__ == "__main__":
    main()
