import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
from six import BytesIO

def add_figure_to_archive(fig, zipfile, filename):
    bytes_buf = BytesIO()
    plt.savefig(bytes_buf, format='png')
    bytes_buf.seek(0)
    zipfile.writestr(filename, bytes_buf.read())
    bytes_buf.close()

class BasicPlots():

    def __init__(self, name):
        self.name = name

    def uniqual_boxplot(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots(figsize = (10,8))
        if color != "blue":
            white = False
        else:
            white = True
        axes.boxplot(data, sym = "", patch_artist = white, meanline=True, showmeans=True)
        if not grid:
            axes.grid()
        axes.set_title("Quality per Cycle " + sample)
        axes.set_xlabel("Cycle")
        axes.set_ylabel("Quality level")
        add_figure_to_archive(fig, filename, "cycle-quality.png")

    def uniqualitysec_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        max_len = max(data.keys())
        x_axis = tuple(range(max_len + 1))
        axes.bar(x_axis, [data[n] for n in x_axis], color=color)
        if not grid:
            axes.grid()
        axes.set_title('Mean Sequence Quality ' + sample)
        axes.set_xlabel('Quality level')
        axes.set_ylabel('Number of sequences')
        add_figure_to_archive(fig, filename, "sequence-quality.png")

    def uninuc_plot(self, data, filename, grid, sample):
        fig, axes = plt.subplots()
        axes.plot(data)
        if not grid:
            axes.grid()
        axes.set_title("Nucleotide per Cycle " + sample)
        axes.set_ylabel("Nucleotide percentage (%)")
        axes.set_xlabel("Cycle")
        axes.axis([0, data.shape[0], 0, 100])
        axes.legend(["A", "T", "C", "G", "N"])
        add_figure_to_archive(fig, filename, "nucleotide-percentage.png")

    def unigc_plot(self, data, filename, grid, sample):
        fig, axes = plt.subplots()
        axes.plot(data)
        if  not grid:
            axes.grid()
        axes.set_title("GC per Cycle " + sample)
        axes.axis([0, len(data), 0, 100])
        axes.set_ylabel("GC percentage (%)")
        axes.set_xlabel("Cycle")
        add_figure_to_archive(fig, filename, "gc-percentage.png")

    def unigcproportion_scatter(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        gc_list = []

        for key in data.keys():
            for counter in range(data[key]):
                gc_list.append(key)

        mean = np.mean(gc_list)
        sigma = np.std(gc_list)
        second_x = np.linspace(0,100,100)
        x_axis = tuple(range(101))
        y_axis = [data[height] for height in x_axis]
        axes.scatter(x_axis, y_axis, color=color)
        if not grid:
            axes.grid()
        x1,x2,y1,y2 = axes.axis()
        axes.axis((x1,x2,0,y2))
        axes2 = axes.twinx()
        axes2.plot(second_x,mlab.normpdf(second_x,mean,sigma), color='red')
        axes2.get_yaxis().set_visible(False)
        handles, labels = axes.get_legend_handles_labels()
        display = (0,1,2)
        a = plt.Line2D((0,1),(0,0), color=(0.1,0.6,0.8))
        b = plt.Line2D((0,1),(0,0), color='red')
        axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
        axes.set_axisbelow(True)
        axes.set_title('Sequence GC Content Distribution ' + sample)
        axes.set_xlabel('Mean sequence GC content (%)')
        axes.set_ylabel('Number of sequences')
        axes.legend([handle for i,handle in enumerate(handles) if i in display]+[a,b],
                             [label for i,label in enumerate(labels) if i in display]+['Observed','Theoretical'])
        add_figure_to_archive(fig, filename, "gc-sequence-distribution.png")

    def unilenght_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        max_len = max(data.keys())
        x_axis = tuple(range(max_len + 1))
        axes.bar(x_axis, [data[n] for n in x_axis], color=color)
        if not grid:
            axes.grid()
        axes.set_title('Sequence Lengths ' + sample)
        axes.set_xlabel('Lenght')
        axes.set_ylabel('Number of sequences')
        add_figure_to_archive(fig, filename, "sequence-lengths.png")

    def unioverkmer_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots(figsize = (10,8))
        x_axis = range(len(data))
        y_axis = []
        x_proaxis = []

        for i in data:
            x_proaxis.append(i[0])
            y_axis.append(i[1])

        axes.bar(x_axis, y_axis, color=color)
        if not grid:
            axes.grid()
        plt.xticks(x_axis, x_proaxis, rotation = "vertical")
        axes.set_title("Over-Represented Kmers " + sample)
        axes.set_ylabel("Number of repeats")
        axes.set_xlabel("Kmer")
        axes.set_xticks(x_axis)
        axes.set_xticklabels(x_proaxis)
        add_figure_to_archive(fig, filename, "overrepresented-kmer.png")

    def unikmer_plot(self, data, over, filename, grid, sample):
        fig, axes = plt.subplots()
        x_proaxis = []
        for i in over:
            x_proaxis.append(i[0])
        axes.plot(data)
        if not grid:
            axes.grid()
        axes.legend(x_proaxis)
        axes.grid()
        axes.set_ylabel("Number of Repeats")
        axes.set_xlabel("Cycle")
        axes.set_title("Over-represented Kmers Per Cycle " + sample)
        add_figure_to_archive(fig, filename, "kmer-cycle.png")

    def uniduplicants_hist(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        label = data[2]
        data = data[3:]
        x_axis = range(3, len(data) + 3)
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid()
        axes.set_title("Duplicated Sequences " + sample)
        axes.set_xlabel("Repeats" + "     " + "Sequences with 2 repeats: " + str(label))
        axes.set_ylabel("Number of sequences")
        add_figure_to_archive(fig, filename, "sequences-duplicated.png")

    def unisecn_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        x_axis = range(1,(len(data)+1))
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid()
        axes.set_title("Ns per sequence " + sample)
        axes.set_ylabel("Number of sequences")
        axes.set_xlabel("Number of Ns")
        add_figure_to_archive(fig, filename, "sequence-n.png")

    def unicyclen_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        x_axis = range(1,len(data) + 1)
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid()
        axes.set_title("Ns per cycle " + sample)
        axes.set_ylabel("Number of Ns")
        axes.set_xlabel("Cycle")
        add_figure_to_archive(fig, filename, "cycle-n.png")

class SamPlots(BasicPlots):

    def __init__(self, name):
        self.name = name

    def uniflag_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots(figsize = (10,8))
        x_axis = range(1,13)
        x_proaxis = [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000, 100000000000]
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid
        plt.xticks(x_axis, x_proaxis, rotation = "vertical")
        axes.set_title("Flags bits " + sample)
        axes.set_ylabel("Number of sequences")
        axes.set_xlabel("Flags bits")
        add_figure_to_archive(fig, filename, "sequence-flag.png")

    def unicigar_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        x_axis = range(8)
        x_proaxis = ["I", "D", "N", "S", "H", "P", "=", "X"]
        axes.bar(x_axis, data[1:], color=color)
        if not grid:
            axes.grid
        plt.xticks(x_axis, x_proaxis)
        axes.set_title("Cigars " + sample)
        axes.set_ylabel("Number of cigars")
        axes.set_xlabel("Cigars" + "  " + "Number of cigars M: " + str(data[0]))
        add_figure_to_archive(fig, filename, "cigars.png")

    def uniposition_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        x_axis = range(1,(len(data)+1))
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid()
        axes.set_title("Mapping positions " + sample)
        axes.set_ylabel("Number of sequences")
        axes.set_xlabel("Position")
        add_figure_to_archive(fig, filename, "position.png")

    def unimapq_bar(self, data, filename, grid, color, sample):
        fig, axes = plt.subplots()
        x_axis = range(len(data))
        axes.bar(x_axis, data, color=color)
        if not grid:
            axes.grid()
        axes.set_title("Mapping Quality " + sample)
        axes.set_ylabel("Number of sequences")
        axes.set_xlabel("Mapping quality")
        add_figure_to_archive(fig, filename, "mapping-quality.png")
