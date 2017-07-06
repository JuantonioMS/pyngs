import argparse
import sys
import time
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import tempfile
import collections

"""CLASE FASTQ"""

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

"""CLASE FASTA"""

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

"""CLASE SAM"""

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
        self.N_chain = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        self.phred = phred_value
        if self.phred is not None:
            self.asciiconversion()

    def __str__(self):
        return "\t".join([self.qname, str(self.flag), self.rname, str(self.position), str(self.mapq), str(self.cigar), self.rnext, str(self.pnext), str(self.length), self.sequence, self.quality]) + "\n"

    def __len__(self):
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
            return Fasta(name, sequence)

        except StopIteration:
            raise StopIteration

    def sampling(self):
        n = 4
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
