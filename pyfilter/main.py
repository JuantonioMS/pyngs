import argparse
import sys
import time
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import tempfile
import collections

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

class FastQ(object):

    def __init__(self, name, sequence, extra, quality, phred_value=None):
        self.name = name
        self.sequence = sequence
        self.extra = extra
        self.quality = quality
        self.N_chain = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        self.phred = phred_value
        if self.phred is not None:
            self.asciiconversion()



    def __str__(self):
        return "\n".join([self.name, self.sequence, self.extra, self.quality]) + "\n"

    def __len__(self, qualen = False):
        if qualen:
            return len(self.quality)

        else:
            return len(self.sequence)

    def qualitysec(self):
        qualities = [(ord(character) - self.phred) for character in self.quality]
        return sum(qualities) / (len(qualities))


    def gc_proportion(self):
        g_number = self.sequence.count("G")
        c_number = self.sequence.count("C")
        proportion = float(g_number + c_number) / len(self.sequence)
        return int(proportion * 100)

    def nucproportion(self, nuc):
        proportion = self.sequence.count(nuc) / float(len(self.sequence)) * 100
        return int(proportion)

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

    def asciiconversion(self):
        self.converted = [(ord(character) - self.phred) for character in self.quality]
        return self.converted

    def minquality(self, min_quality):
        position = None
        while position == None:
            try:
                position = self.converted.index(min_quality - 1)
            except ValueError:
                min_quality = min_quality + 1
                if min_quality > 41:
                    sys.stdout.write("Min set over 41")
                    break
        transformed_seq = self.sequence[:position]
        transformed_qual = self.quality[:position]
        if len(transformed_seq) == 0:
            return False, transformed_seq, transformed_qual
        else:
            return True, transformed_seq, transformed_qual

    def maxquality(self, max_quality):
        position = None
        while position == None:
            print max_quality
            try:
                position = self.converted.index(max_quality)
            except ValueError:
                max_quality = max_quality - 1
                if max_quality < 0:
                    sys.stdout.write("Max set to 0")
                    break
        transformed_seq = self.sequence[position:]
        transformed_qual = self.quality[position:]
        if len(transformed_seq) == 0:
            return False, transformed_seq, transformed_qual
        else:
            return True, transformed_seq, transformed_qual


    def nucnumber(self, nuc_number):
        if nuc_number <= len(self.sequence):
            transformed_seq = self.sequence[:nuc_number]
            transformed_qual = self.quality[:nuc_number]
            return True, transformed_seq, transformed_qual

        else:
            diference = nuc_number - len(self.sequence)
            transformed_seq = self.sequence + self.N_chain[:diference]
            transformed_qual = self.quality

            while len(transformed_seq) != len(transformed_qual):
                transformed_qual = self.quality + self.quality[-1]

            return True, transformed_seq, transformed_qual

    def gcpercentage(self, gc_percentage):
        if self.gc_proportion() <= gc_percentage:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def GCpercentage(self, gc_percentage):
        if self.gc_proportion() >= gc_percentage:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def meanquality(self, mean_quality):
        if self.qualitysec() <= mean_quality:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def MEANquality(self, mean_quality):
        if self.qualitysec() >= mean_quality:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def nucpercentage(self, nuc, nuc_percentage):
        if int(nuc_percentage) >= self.nucproportion(nuc):
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def NUCpercentage(self, nuc, nuc_percentage):
        if nucproportion(nuc) >= int(nuc_percentage):
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def removeduplicants(self,data):
        if data[self.sequence] >= 2:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def POLYXfilter(self, size):
        counter_a = 0
        counter_t = 0

        for nuc in self.sequence:

            if nuc == "T":
                counter_a = 0
                counter_t += 1

            elif nuc == "A":
                counter_a += 1
                counter_t = 0

            elif nuc != "A" and nuc != "T":
                counter_a = 0
                counter_t = 0

            if counter_a == size or counter_t == size:
                return False, self.sequence, self.quality
                break

        return True, self.sequence, self.quality

    def polyxfilter(self, size):
        counter_a = 0
        counter_t = 0

        for nuc in self.sequence:

            if nuc == "T":
                counter_a = 0
                counter_t += 1

            elif nuc == "A":
                counter_a += 1
                counter_t = 0

            elif nuc != "A" and nuc != "T":
                counter_a = 0
                counter_t = 0

            if counter_a == size:
                add_chain = ""
                while len(add_chain) != size:
                    add_chain = add_chain + "A"
                return True, self.sequence.replace(add_chain,self.N_chain[:size]), self.quality
                break

            if counter_t == size:
                add_chain = ""
                while len(add_chain) != size:
                    add_chain = add_chain + "T"
                return True, self.sequence.replace(add_chain,self.N_chain[:size]), self.quality
                break

        return True, self.sequence, self.quality


    def gcregion(self, size):
        counter = 0
        for nuc in self.sequence:
            if nuc == "G" or nuc == "C":
                counter += 1
            else:
                counter = 0
            if counter == size:
                return False, self.sequence, self.quality
                break
        return True, self.sequence, self.quality

    def qualitynchanger(self, min_quality):
        transformed_seq = ""
        for position, (nuc, qual) in enumerate(zip(self.sequence, self.converted)):
            if qual <= min_quality:
                transformed_seq = transformed_seq + "N"
            else:
                transformed_seq = transformed_seq + self.sequence[position]
        return True, transformed_seq, self.quality
    def Qualitynchanger(self, max_quality):
        transformed_seq = ""
        for position, (nuc, qual) in enumerate(zip(self.sequence, self.converted)):
            if qual >= max_quality:
                transformed_seq = transformed_seq + "N"
            else:
                transformed_seq = transformed_seq + self.sequence[position]
        return True, transformed_seq, self.quality

    def complement(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        complement = str(seq.complement())
        return True, complement, self.quality

    def reversecomplement(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        reverse_complement = str(seq.reverse_complement())
        return True, reverse_complement, reverse(self.quality)

    def transcription(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        transcript = str(seq.transcribe())
        return True, transcript, self.quality

    def retrotranscription(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_rna)
        retro_transcript = str(seq.back_transcribe())
        return True, retro_transcript, self.quality

    def righttrimmer(self, size):
        return True, self.sequence[:-size], self.quality[:-size]

    def lefttrimmer(self, size):
        return True, self.sequence[size:], self.quality[size:]

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

    def sampling(self, phred_value):
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
                yield FastQ(name=name, sequence=seq, extra=strand, quality=qual, phred_value=phred_value)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.filename.close()

class Fasta(object):

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return "\n".join([self.name, self.sequence]) + "\n"

    def nucproportion(self, nuc):
        proportion = self.sequence.count(nuc) / float(len(self.sequence)) * 100
        return int(proportion)

    def righttrimmer(self, size):
        return True, self.sequence[:-size]

    def lefttrimmer(self, size):
        return True, self.sequence[size:]

    def removenucpercentage(self, nuc, nuc_percentage):
        if self.nucproportion(nuc) <= int(nuc_percentage):
            return False, self.sequence
        else:
            return True, self.sequence

    def removeNUCpercentage(self, nuc, nuc_percentage):
        if self.nucproportion(nuc) >= int(nuc_percentage):
            return False, self.sequence
        else:
            return True, self.sequence

    def removeambiguous(self):
        for character in self.sequence:
            if caracter in "RYKMSWBDHVNX-":
                return False, self.sequence
            else:
                return True, self.sequence

    def removeduplicants(self, data):
        if data[self.sequence] >= 2:
            return False, self.sequence, self.quality

        else:
            return True, self.sequence, self.quality

    def complement(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        complement = str(seq.complement())
        return True, complement, self.quality

    def reversecomplement(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        reverse_complement = str(seq.reverse_complement())
        return True, reverse_complement

    def transcription(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        transcript = str(seq.transcribe())
        return True, transcript

    def retrotranscription(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_rna)
        retro_transcript = str(seq.back_transcribe())
        return True, retro_transcript


class FastaReader:

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

            return FastQ(name, sequence)

        except StopIteration:
            raise StopIteration

    def sampling(self, phred_value):
        n = 2
        for i, line in enumerate(self.filename):

            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
                yield Fasta(name=name, sequence=seq)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.filename.close()


def asciidetection(quality):
    low_list = [chr(character) for character in range(33,59)]
    high_list = [chr(character) for character in range(74,127)]
    phred_value = int()
    for ascii_letter in quality:

        if ascii_letter in low_list:
            phred_value = 33

        elif ascii_letter in high_list:
            phred_value = 64

        if phred_value:
            return phred_value
            break

def reverse(list):
	if len(list)==1:
		return list
	else:
		return list[-1]+reverse(list[:-1])

def switch(firstflag, secondflag):
    if firstflag == True and secondflag == False:
        firstflag = False
        secondflag = True
    elif firstflag == False and secondflag == True:
        firstflag = True
        secondflag = False
    else:
        sys.exit("Warning! Switch failure!")
    return firstflag, secondflag

def debug(VERBOSE_LEVEL, *args, **kwargs):

    if VERBOSE_LEVEL == DEBUG_LEVEL:
        print DEBUG_LEVEL + "! " + "".join(str(x) for x in args)


def run(arguments):
    print "Iniciando el programa"

    start_time = time.time()
    debug("Flags", "arguments", arguments)

    arguments['input'] = argparse.FileType('r')(arguments['input'])
    arguments['text'] = argparse.FileType('w')(arguments['text'])
    args = Bunch(arguments)
    debug("Flags", "args", args)
    debug("Flags", "argument list", sys.argv)


    name, ext = args.input.name.split(".")
    sys.stdout.write("File %s loaded, extension .%s detected\n" % (args.input.name, ext))

    if ext not in ("fq", "fastq", "fasta", "sam"):
        sys.exit("Wrong extension!\n" + "Please, check your file extension is (fq, fastq, fasta, sam).")

    if ext in ("fq", "fastq") and not args.testpool:
        debug("Flags", "FastQ section enter")

        try:
            estimate_sequences = total_lines / 4
            sys.stdout.write("Estimated %i sequences\n" % estimate_sequences)

        except:
            with FastqReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1
            if args.duplicants:
                store_duplicants = collections.defaultdict(int)

        with FastqReader(open(args.input.name)) as loaded_file:
            phred_value = 0
            for sequences in loaded_file:
                if phred_value == 0:
                    phred_value = asciidetection(sequences.quality)
                if args.duplicants:
                    store_duplicants[sequences] += 1





        temporalfileone = tempfile.TemporaryFile()
        temporalfiletwo = tempfile.TemporaryFile()

        flag_one = True
        flag_two = False
        flag_original = True
        flag_operation = False

        for argument in sys.argv:
            temporalfileone.seek(0)
            temporalfiletwo.seek(0)

            if not flag_operation:
                debug("Flag", "Original file input")
                infile = FastqReader(args.input)
                reads = infile.sampling(phred_value)
                flag_original = False
            elif flag_one == False and flag_original == False:
                debug("Flag", "First temporal file input")
                infile = FastqReader(temporalfileone)
                reads = infile.sampling(phred_value)
            elif flag_two == False and flag_original == False:
                debug("Flag", "Second temporal file input")
                infile = FastqReader(temporalfiletwo)
                reads = infile.sampling(phred_value)


            if argument in ["-q", "--minquality"]:
                debug("Flag", "Minquality treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.minquality(args.minquality)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-Q", "--maxquality"]:
                debug("Flag", "Maxquality treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.maxquality(args.maxquality)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-n", "--nucnumber"]:
                debug("Flag", "Nucnumber treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.nucnumber(args.nucnumber)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-gc", "--gcpercentage"]:
                debug("Flag", "gcpercentage treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.gcpercentage(args.gcpercentage)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-GC", "--GCpercentage"]:
                debug("Flag", "GCpercentage treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.GCpercentage(args.GCpercentage)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-mq", "--meanquality"]:
                debug("Flag", "meanquality treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.meanquality(args.meanquality)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-MQ", "--MEANquality"]:
                debug("Flag", "MEANquality treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.MEANquality(args.MEANquality)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-nuc", "--nucpercentage"]:
                debug("Flag", "nucpercentage treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.nucpercentage(args.nucpercentage[0], args.nucpercentage[1])
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-NUC", "--NUCpercentage"]:
                debug("Flag", "nucpercentage treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.NUCpercentage(args.NUCpercentage[0], args.NUCpercentage[1])
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-d", "--duplicants"]:
                debug("Flag", "duplicants treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.removeduplicants(store_duplicants)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-P", "--POLYXfilter"]:
                debug("Flag", "POLYXfilter treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.POLYXfilter(args.POLYXfilter)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-p", "--polyxfilter"]:
                debug("Flag", "polyxfilter treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.polyxfilter(args.polyxfilter)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-gr", "--gcregion"]:
                debug("Flag", "gcregion treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.gcregion(args.gcregion)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-qc", "--qualitynchanger"]:
                debug("Flag", "qualitynchanger treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.qualitynchanger(args.qualitynchanger)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-QC", "--QUALITYnchanger"]:
                debug("Flag", "QUALITYnchanger treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.Qualitynchanger(args.QUALITYnchanger)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-c", "--complement"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.complement()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-t", "--transcription"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.transcription()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-r", "--reversecomplement"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.reversecomplement()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-lt", "--lefttrimmer"]:
                debug("Flag", "lefttrimmer treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.lefttrimmer(args.lefttrimmer)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-rt", "--righttrimmer"]:
                debug("Flag", "righttrimmer treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence, quality = read.righttrimmer(args.righttrimmer)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n" + read.extra + "\n" + quality + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

    if ext in ("fasta") and not args.testpool:
        debug("Flags", "Fasta section enter")

        try:
            estimate_sequences = total_lines / 4
            sys.stdout.write("Estimated %i sequences\n" % estimate_sequences)

        except:
            with FastaReader(open(args.input.name)) as loaded_file:
                estimate_sequences = 0
                for sequences in loaded_file:
                    estimate_sequences += 1
        if args.duplicants:
            store_duplicants = collections.defaultdict(int)
            with FastaReader(open(args.input.name)) as loaded_file:
                for sequences in loaded_file:
                    store_duplicants[sequences] += 1

        temporalfileone = tempfile.TemporaryFile()
        temporalfiletwo = tempfile.TemporaryFile()

        flag_one = True
        flag_two = False
        flag_original = True
        flag_operation = False

        for argument in sys.argv:
            temporalfileone.seek(0)
            temporalfiletwo.seek(0)

            if not flag_operation:
                debug("Flag", "Original file input")
                infile = FastqReader(args.input)
                reads = infile.sampling(phred_value)
                flag_original = False
            elif flag_one == False and flag_original == False:
                debug("Flag", "First temporal file input")
                infile = FastqReader(temporalfileone)
                reads = infile.sampling(phred_value)
            elif flag_two == False and flag_original == False:
                debug("Flag", "Second temporal file input")
                infile = FastqReader(temporalfiletwo)
                reads = infile.sampling(phred_value)

            if argument in ["-lt", "--lefttrimmer"]:
                debug("Flag", "lefttrimmer treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.lefttrimmer(args.lefttrimmer)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-rt", "--righttrimmer"]:
                debug("Flag", "righttrimmer treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.righttrimmer(args.righttrimmer)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-d", "--duplicants"]:
                debug("Flag", "duplicants treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.removeduplicants(store_duplicants)
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-c", "--complement"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.complement()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-t", "--transcription"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.transcription()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-ra", "--removeambiguous"]:
                debug("Flag", "removeambiguous treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.removeambiguous()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)

            if argument in ["-r", "--reversecomplement"]:
                debug("Flag", "complement treatment")
                flag_operation = True
                for read in reads:
                    writable, sequence = read.reversecomplement()
                    if writable:
                        if flag_one:
                            debug("Flag", "Writing over first tamporal file")
                            temporalfileone.write(read.name + "\n" + sequence + "\n")
                        elif flag_two:
                            debug("Flag", "Writing over second tamporal file")
                            temporalfiletwo.write(read.name + "\n" + sequence + "\n")
                if flag_one:
                    debug("Flag", "Cleaning second temporal file")
                    temporalfiletwo.close()
                    temporalfiletwo = tempfile.TemporaryFile()
                elif flag_two:
                    debug("Flag", "Cleaning first temporal file")
                    temporalfileone.close()
                    temporalfileone = tempfile.TemporaryFile()

                flag_one, flag_two = switch(flag_one, flag_two)



        temporalfileone.seek(0)
        temporalfiletwo.seek(0)

        for i in range(20):
            lecture1 = temporalfileone.read()
            lecture2 = temporalfiletwo.read()
            print lecture1, "/", lecture2

def main():
    parser = argparse.ArgumentParser(prog='filter', description="Bioinformatics NGS filter")
    parser.add_argument('input', type=str, help="input file (Fastq, Fasta, Sam")
    parser.add_argument('-v', '--verbose', type=str, dest='DEBUG_LEVEL', action='store', default=None, help='Level screen output (default: %(default)s)')
    parser.add_argument('-e', '--text', type=str, default='-', help="file name for text output (default: %(default)s)")

    parser.add_argument('-q', '--minquality', type=int, help='minimum quality accepted')
    parser.add_argument('-Q', '--maxquality', type=int, help='minimum quality accepted')
    parser.add_argument('-n', '--nucnumber', type=int, help='number of nuc per sequence (lenght)')
    parser.add_argument('-gc', '--gcpercentage', type=int, help='min gc content per sequence')
    parser.add_argument('-GC', '--GCpercentage', type=int, help='max gc content per sequence')
    parser.add_argument('-mq', '--meanquality', type=int, help='min quality mean per sequence')
    parser.add_argument('-MQ', '--MEANquality', type=int, help='max quality mean per sequence')
    parser.add_argument('-nuc', '--nucpercentage', nargs=2, help='min nuc content per sequence')
    parser.add_argument('-NUC', '--NUCpercentage', type=tuple, help='max nuc content per sequence')
    parser.add_argument('-d', '--duplicants', action="store_true", default=False, help='remove duplicated sequences')
    parser.add_argument('-P', '--POLYXfilter', type=int, help='remove sequences with polyx region')
    parser.add_argument('-p', '--polyxfilter', type=int, help='change polyx region with N')
    parser.add_argument('-gr', '--gcregion', type=int, help='remove sequences with gcregions(size)')
    parser.add_argument('-qc', '--qualitynchanger', type=int, help='change to N if not minimal quality')
    parser.add_argument('-QC', '--QUALITYnchanger', type=int, help='change to N if higher than max quality')
    parser.add_argument('-c', '--complement', action="store_true", default=False, help='change sequence to complement')
    parser.add_argument('-t', '--transcription', action="store_true", default=False, help='change DNA to RNA')
    parser.add_argument('-r', '--reversecomplement', action="store_true", default=False, help='change sequence to reverse complement')
    parser.add_argument('-rt', '--righttrimmer', type=int, help='number of nuc right trimmed')
    parser.add_argument('-lt', '--lefttrimmer', type=int, help='number of nuc left trimmed')
    parser.add_argument('-ra', '--removeambiguous', action="store_true", default=False, help='remove non standar bases (A,C,G,T)')
    parser.add_argument('-tp', '--testpool', action="store_true", default=False, help='testing diferents traitments, make .csv file')

    args = parser.parse_args()
    globals().update(args.__dict__)
    arguments = vars(args)
    run(arguments)

if __name__ == "__main__":
    main()