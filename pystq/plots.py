import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class BasicPlots1():

    def __init__(self, name):
        self.name = name

    def uniqual_boxplot(self, data):
        fig, axes = plt.subplots(figsize = (10,8))
        axes.boxplot(data, sym = "", patch_artist = True)
        axes.set_title("Quality Per Cycle")
        axes.set_xlabel("Cycle")
        axes.grid()
        axes.set_ylabel("Quality level")
        fig.savefig("qualcycle.png")

    def uniqualitysec_bar(self, data):
        fig, axes = plt.subplots()
        max_len = max(data.keys())
        x_axis = tuple(range(max_len + 1))
        axes.bar(x_axis, [data[n] for n in x_axis], color="red")
        axes.set_title('Mean Sequence Quality')
        axes.set_xlabel('Quality level')
        axes.set_ylabel('Number of sequences')
        fig.savefig("qualitysec.png")

    def uninuc_plot(self, data):
        fig, axes = plt.subplots()
        axes.plot(data)
        axes.grid()
        axes.set_title("NUC per Cycle")
        axes.set_ylabel("NUC percentage (%)")
        axes.set_xlabel("Cycle")
        axes.axis([0, data.shape[0], 0, 100])
        axes.legend(["A", "T", "C", "G"])
        fig.savefig("nuc.png")

    def unigc_plot(self, data):
        fig, axes = plt.subplots()
        axes.plot(data)
        axes.grid()
        axes.set_title("GC per cycle")
        axes.axis([0, len(data), 0, 100])
        axes.set_ylabel("GC percentage (%)")
        axes.set_xlabel("Cycle")
        fig.savefig("gcnuc.png")

    def unigcproportion_scatter(self, data):
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
        axes.scatter(x_axis, y_axis)
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
        axes.set_title('Sequence GC content distribution')
        axes.set_xlabel('Mean sequence GC content (%)')
        axes.set_ylabel('Number of sequences')
        axes.legend([handle for i,handle in enumerate(handles) if i in display]+[a,b],
                             [label for i,label in enumerate(labels) if i in display]+['Observed','Theoretical'])
        fig.savefig("gcdistribution.png")

    def unilenght_bar(self, data):
        fig, axes = plt.subplots()
        max_len = max(data.keys())
        x_axis = tuple(range(max_len + 1))
        axes.bar(x_axis, [data[n] for n in x_axis])
        axes.set_title('Sequence lengths')
        axes.set_xlabel('Lenght')
        axes.set_ylabel('Number of sequences')
        axes.grid()
        fig.savefig("lenghts.png")

    def unioverkmer_bar(self, data):
        fig, axes = plt.subplots(figsize = (10,8))
        x_axis = range(len(data))
        y_axis = []
        x_proaxis = []

        for i in data:
            x_proaxis.append(i[0])
            y_axis.append(i[1])

        axes.bar(x_axis, y_axis)
        plt.xticks(x_axis, x_proaxis, rotation = "vertical")
        axes.set_title("Over-represented Kmers")
        axes.set_ylabel("Number of repeats")
        axes.set_xlabel("Kmer")
        axes.set_xticks(x_axis)
        axes.set_xticklabels(x_proaxis)
        fig.savefig("overkmer.png")

    def unikmer_plot(self, data, over):
        fig, axes = plt.subplots()
        x_proaxis = []
        for i in over:
            x_proaxis.append(i[0])
        axes.plot(data)
        axes.legend(x_proaxis)
        axes.grid()
        axes.set_ylabel("Number of repeats")
        axes.set_xlabel("Cycle")
        axes.set_title("Over-represented Kmers Per Cycle")
        fig.savefig("kmer.png")

    def uniduplicants_hist(self, data):
        fig, axes = plt.subplots()
        new_list = []
        values = []
        for value in data.values():
            if value != 1:
                values.append(value)
        for i in values:
            if i < 20:
                new_list.append(i)
        x_axis = range(2,max(new_list)+1)
        axes.hist(new_list)
        axes.set_xticks(x_axis)
        axes.set_title("Duplicated Sequences")
        axes.set_xlabel("Repeats")
        axes.set_ylabel("Number of sequences")
        fig.savefig("duplicants.png")

    def unisecn_plot(self, data):
        fig, axes = plt.subplots()
        x_axis = range(1,(len(data)+1))
        axes.plot(x_axis, data)
        axes.grid()
        axes.set_title("Ns per sequence")
        axes.set_ylabel("Number of sequences")
        axes.set_xlabel("Number of Ns")
        fig.savefig("secn.png")

    def unicyclen_plot(self, data):
        fig, axes = plt.subplots()
        axes.plot(data)
        axes.grid()
        axes.set_title("Ns per cycle")
        axes.set_ylabel("number of Ns")
        axes.set_xlabel("Cycle")
        fig.savefig("cyclen.png")

class SamPlots(BasicPlots):

    def __init__(self, name):
        self.name = name

    def uniflag_bar(self, data):
        fig, axes = plt.subplots()
        x_axis = range(1,13)
        axes.bar(x_axis, data)
        fig.savefig("prueba.png")
