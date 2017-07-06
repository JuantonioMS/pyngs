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

    def __len__(self):
        return len(self.sequence)

    def flaglist(self, store):
        str_flag = bin(self.flag)
        str_flag = str_flag[2:]
        list_flag = [int(b) for b in str_flag]

        for position, value in enumerate(list_flag):
            store[position] += value

        return store

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

    def cigarlist(self, store):
        if self.cigar == "*":
            return store
        else:
            scafold = ""
            for character in self.cigar:
                try:
                    int(character)
                    scafold = "".join([scafold, character])
                except ValueError:
                    number = int(scafold)
                    reserve = store.index(character) - 1
                    store[reserve] += number
                    scafold = ""
                    continue
            return store


class SamReader:

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

class Fasta(object):

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return "\n".join([self.name, self.sequence]) + "\n"

    def __len__(self):
        return len(self.sequence)

    def gc_proportion(self):
        g_number = self.sequence.count("G")
        c_number = self.sequence.count("C")
        proportion = float(g_number + c_number) / len(self.sequence)
        return int(proportion * 100)
