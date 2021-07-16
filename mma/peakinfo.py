import os
import xml.etree.ElementTree

class Adduct:
    def __init__(self, mw, name, fragment, ms2specdata, metabolites_dict, ms2_files):

        self.name = name
        self.fragment = fragment
        self.mw = mw
        self.hmdb = self.get_hmdb(metabolites_dict)
        if self.hmdb:
            self.score, self.ms2spec = self.get_score(ms2specdata, ms2_files, metabolites_dict)
        else:
            self.score, self.ms2spec = 0, ''

    def get_hmdb(self, metabolites_dict, tolerance = 0.3):

        for metabolite in metabolites_dict:

            mmw = metabolites_dict[metabolite][0]
            if  mmw != None:
                if self.mw >= float(mmw) - tolerance and self.mw <= float(mmw) + tolerance:

                    return metabolite

    def get_score(self, ms2specdata, ms2_files, metabolites):

        import glob

        best_score = 0.0
        best_ms2spec = ''
        files = glob.glob(ms2_files + self.hmdb +'*')

        if len(files)>0:

            for file in files:
                ms2spec = self.get_ms2_spec_from_hmdb(file, float(metabolites[self.hmdb][0]))

                if ms2spec != None: #msms is positive

                        for spectra in ms2specdata:
                            ms2 = ms2specdata[spectra]
                            if ms2 != None:
                                score, match = self.fast_cosine(ms2spec, ms2, 0.2, 2)

                                if score > best_score:
                                    best_score = score
                                    best_ms2spec = ms2spec

        return best_score, best_ms2spec

    def get_ms2_spec_from_hmdb(self, file, parent_mz):
        import sys
        sys.path.append('/Users/anamaria/git/molnet/code/')
        import mnet
        path_to_hmdbfile = file
        et = xml.etree.ElementTree.parse(path_to_hmdbfile)
        element = et.getroot()
        mode = element.find('ionization-mode').text
        if mode == 'positive':
            instrtype = element.find('instrument-type').text
            filename = element.find('database-id').text
            np = element.find('peak-counter').text
            peaks = []

            for msms in element.find('ms-ms-peaks'):
                mz = float(msms.find('mass-charge').text)
                intensity = float(msms.find('intensity').text)

                peaks.append((mz, intensity))


            ms2_spectrum = mnet.Spectrum(peaks, filename, None, None, parent_mz, parent_mz, metadata = (instrtype, mode))
            return ms2_spectrum


    def find_pairs(self, spec1,spec2,tol,shift=0):
        #from Simon's molnet code
        matching_pairs = []
        spec2lowpos = 0
        spec2length = len(spec2)

        for idx,(mz,intensity) in enumerate(spec1):
            # do we need to increase the lower idx?
            while spec2lowpos < spec2length and spec2[spec2lowpos][0] + shift < mz - tol:
                spec2lowpos += 1
            if spec2lowpos == spec2length:
                break
            spec2pos = spec2lowpos
            while(spec2pos < spec2length and spec2[spec2pos][0] + shift < mz + tol):
                matching_pairs.append((idx,spec2pos,intensity*spec2[spec2pos][1]))
                spec2pos += 1

        return matching_pairs

    def fast_cosine(self, spectrum1,spectrum2,tol,min_match):
        #from Simon's molnet code
        # spec 1 and spec 2 have to be sorted by mz
        if spectrum1.n_peaks == 0 or spectrum2.n_peaks == 0:
            return 0.0,[]
        # find all the matching pairs

        spec1 = spectrum1.normalised_peaks
        spec2 = spectrum2.normalised_peaks

        matching_pairs = self.find_pairs(spec1,spec2,tol,shift = 0.0)
        matching_pairs = sorted(matching_pairs,key = lambda x:x[2],reverse = True)
        used1 = set()
        used2 = set()
        score = 0.0
        used_matches = []
        for m in matching_pairs:
            if not m[0] in used1 and not m[1] in used2:
                score += m[2]
                used1.add(m[0])
                used2.add(m[1])
                used_matches.append(m)
        if len(used_matches) < min_match:
            score = 0.0
        return score,used_matches



class AdductList:

    def __init__(self, mz, ms2specdata, metabolites_dict, ms2_files, positive):

        self.adducts = []
        adduct_list = self.compute_adducts(mz, positive)
        for i in adduct_list:
            adduct = Adduct(i[0], i[1], i[2], ms2specdata, metabolites_dict, ms2_files)
            self.adducts.append(adduct)

    def print_adducts(self):
        for adduct in self.adducts:
            print(adduct.name, adduct.mw, adduct.fragment, adduct.hmdb, adduct.score)

    def get_adduct_with_best_score(self):
        score = 0.0
        final_adduct = ''
        for adduct in self.adducts:
            if adduct.score > score:
                score = adduct.score
                final_adduct = adduct
        return final_adduct

    def compute_adducts(self, mz, positive = True):

        PROTON = 1.00727646677
        if positive:
            addList = [(mz - PROTON, 'M+H[1+]', ''),
                 ((mz - PROTON)*2, 'M+2H[2+]', ''),
                 ((mz - PROTON)*3, 'M+3H[3+]', ''),
                 (mz - 1.0034 - PROTON, 'M(C13)+H[1+]', 'C'),
                 ((mz - 0.5017 - PROTON)*2, 'M(C13)+2H[2+]', 'C'),
                 ((mz - 0.3344 - PROTON)*3, 'M(C13)+3H[3+]', 'C'),
                 (mz -1.9958 - PROTON, 'M(S34)+H[1+]', 'S'),
                 (mz -1.9972 - PROTON, 'M(Cl37)+H[1+]', 'Cl'),
                 (mz - 21.9820 - PROTON, 'M+Na[1+]', ''),
                 ((mz - 10.991 - PROTON)*2, 'M+H+Na[2+]', ''),
                 (mz - 37.9555 - PROTON, 'M+K[1+]', ''),
                 (mz - 18.0106 - PROTON, 'M+H2O+H[1+]', ''),
                 (mz + 18.0106 - PROTON, 'M-H2O+H[1+]', 'H2O'),
                 (mz + 36.0212 - PROTON, 'M-H4O2+H[1+]', 'H4O2'),
                 (mz + 17.0265 - PROTON, 'M-NH3+H[1+]', 'NH3'),
                 (mz + 27.9950 - PROTON, 'M-CO+H[1+]', 'CO'),
                 (mz + 43.9898 - PROTON, 'M-CO2+H[1+]', 'CO2'),
                 (mz + 46.0054 - PROTON, 'M-HCOOH+H[1+]', 'H2CO2'),
                 (mz - 67.9874 - PROTON, 'M+HCOONa[1+]', ''),
                 (mz + 67.9874 - PROTON, 'M-HCOONa+H[1+]', 'HCO2Na'),
                 (mz - 57.9586 - PROTON, 'M+NaCl[1+]', ''),
                 (mz + 72.0211 - PROTON, 'M-C3H4O2+H[1+]', 'C3H4O2'),
                 (mz - 83.9613 - PROTON, 'M+HCOOK[1+]', ''),
                 (mz + 83.9613 - PROTON, 'M-HCOOK+H[1+]', 'HCO2K')]

        else:
            addList = [ ((mz + PROTON)*3, 'M-3H[3-]' ,''),
                       ((mz + PROTON)*2, 'M-2H[2-]' ,''),
                       (mz + 19.01839, 'M-H2O-H' ,''),
                        (mz + PROTON, 'M-H' ,''),
                        (mz - 20.974666, 'M+Na-2H' ,''),
                       (mz - 34.969402, 'M+Cl' ,''),
                       (mz - 36.948606, 'M+K-2H' ,''),
                       (mz - 44.998201, 'M+FA-H' ,''),
                       (mz - 59.013851, 'M+Hac-H' ,''),
                       (mz - 78.918885, 'M+Br' ,''),
                       (mz - 112.985586, 'M+TFA-H' ,''),
                       ((mz + PROTON)/2, '2M-H' ,''),
                       ((mz - 44.998201)/2, '2M+FA-H' ,''),
                       ((mz -59.013851)/2, '2M+Hac-H' ,''),
                       ((mz + PROTON)/3, '3M-H' ,'')]

        return addList


class PeakInfo:

    def __init__(self, cid, mz, rt, mode):
        self.cid = mode+str(cid)
        self.mz = mz
        self.rt = rt
        self.pval = 0.0
        self.tval = 0.0
        self.logfc = 0.0
        self.mm_annotation = ''
        self.mm_pathway = ''
        self.mm_kegg_id = ''
        self.std_annotation = ''
        self.std_kegg_id = ''
        self.spectra = ''
        self.adducts = ''
        self.best_ms2_match_adduct = ''
        self.ms2_annotation = ''
        self.ms2_kegg_id = ''
        self.final_annotation = ''
        self.intensities = {}

    def add_pval(self, pval):
        self.pval = pval

    def add_tval(self, tval):
        self.tval = tval

    def add_logfc(self, logfc):
        self.logfc = logfc

    def add_mm_annotation(self, annotation):
        self.mm_annotation = annotation
    def add_mm_pathway(self, pathway):
        self.mm_pathway = pathway
    def add_mm_kegg_id(self, kegg_id):
        self.mm_kegg_id = kegg_id
    def change_cid(self, mode):
        new_cid = int(self.cid.split(mode)[1].split('.')[0])
        self.cid = new_cid

    def add_intensities(self, intensities):
        self.intensities = intensities


    def get_intensities(self):
        return self.intensities

    def calculate_tolerance(self, x, ppm):

        tolerance = x*ppm/1000000

        return tolerance

    def add_std_match(self, standards_dict, metabolites_dict, ppm = 3, ppm_rt = 30, positive=True):
        import numpy as np
        PROTON = 1.00727646677
        std_annotations = []

        for key in standards_dict.keys():

            for adduct in self.adducts:

                tolerance = self.calculate_tolerance(adduct[0], ppm)
                if positive:
                    mz = float(key) - PROTON
                else:
                    mz = float(key) + PROTON
                if mz - tolerance <= float(adduct[0]) <= mz + tolerance:


                    if standards_dict[key][1]-ppm_rt <= self.rt*60 <= standards_dict[key][1]+ppm_rt:
                        diff = np.abs((float(key)- PROTON) - float(adduct[0]))
                        std_annotations.append(((standards_dict[key][0], adduct[1]), diff))
                        print(self.cid, standards_dict[key][0], adduct[1], diff)
        if len(std_annotations) >1:
             self.std_annotation = sorted(std_annotations, key=lambda x: x[1])[0][0]

        elif len(std_annotations) == 1:
            self.std_annotation =  std_annotations[0][0]
       # for metabolite in metabolites_dict:
       #     if str(self.std_annotation[0]).strip() == str(metabolites_dict[metabolite][2]).lower().strip():
       #         if metabolites_dict[metabolite][5] != None:
       #             self.std_kegg_id = metabolites_dict[metabolite][5]
       #         else:
       #             self.hmdb_kegg_id = metabolite

    def add_spectra(self, aligner_gp):

        #new_peakid = int((self.cid).split('p')[1].split('.')[0]) - 1
        new_peakid = int((self.cid)) - 1

        num_peaks = aligner_gp.peaksets[new_peakid].n_peaks
        spectra = {}
        if num_peaks > 1:
            for i in range(num_peaks-1):
                    source = aligner_gp.peaksets[new_peakid].peaks[i+1].source_file
                    msms = aligner_gp.peaksets[new_peakid].peaks[i+1].ms2_spectrum
                    spectra[source] = msms

        self.spectra = spectra

    def add_adduct_list(self, mode = True):
        self.adducts = AdductList.compute_adducts(self, self.mz, mode)

    def add_adducts(self, metabolites_dict, ms2_files, positive = True):
        self.adducts = AdductList(self.mz, self.spectra, metabolites_dict, ms2_files, positive)
        self.best_ms2_match_adduct = self.adducts.get_adduct_with_best_score()

    def add_ms2_match(self, metabolites_dict, score_theshold = 0.1):
        if self.best_ms2_match_adduct:
            if self.best_ms2_match_adduct.score > score_theshold:
                self.ms2_annotation = metabolites_dict[self.best_ms2_match_adduct.hmdb][2]
                self.ms2_kegg_id = metabolites_dict[self.best_ms2_match_adduct.hmdb][5]

    def get_possible_kegg_ids(self):
        import numpy as np
        kegg_ids = []
        if self.mm_kegg_id != '':
            kegg_list = self.mm_kegg_id.split(';')
            for i in kegg_list:
                kegg_ids.append(i)

        if self.std_kegg_id != '' or self.std_kegg_id != None:
            kegg_ids.append(self.std_kegg_id)

        if self.ms2_kegg_id != '':
            kegg_ids.append(self.ms2_kegg_id)

        kegg_ids = np.array(kegg_ids)
        if len(kegg_ids) >0:
            return kegg_ids
        else:
            return None

    def get_taxonomy(self, metabolites_dict):
        taxonomy = []
        if self.get_possible_kegg_ids() is not None:
            for kegg in self.get_possible_kegg_ids():
                for metabolite in metabolites_dict:
                    if metabolites_dict[metabolite][5] == kegg:
                        taxonomy.append((metabolites_dict[metabolite][2], metabolites_dict[metabolite][8]))
            return taxonomy
        else:
            return None

    def is_standard(self):
        if self.std_kegg_id != '':
            return True
        else:
            return False

    def plot_boxplots(self, data, path):
        import seaborn as sns
        import matplotlib.pyplot as plt


        fig, ax = plt.subplots()


        nice_fonts = {"text.usetex": True,
            "font.family": "serif",
            "font.serif" : "Times New Roman"}
        plt.rcParams.update(nice_fonts)

        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize=19)
        plt.rc('ytick', labelsize=19)



        light_blue = '#34BBE3'
        light_red = '#E3688F'
        pal = {'control': light_blue, 'infected': light_red}

        lighter_red = '#FF75A3'
        lighter_blue = '#39D1FF'
        face_pal = {'control': lighter_blue, 'infected': lighter_red}
        hue_order = ['control', 'infected']

        boxprops = {'edgecolor': 'k', 'linewidth': 2, 'alpha':.5}
        lineprops = {'color': 'k', 'linewidth': 2}

        boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                          'whiskerprops': lineprops, 'capprops': lineprops,
                          'width': 0.75, 'palette': face_pal,
                          'hue_order': hue_order}

        stripplot_kwargs = {'linewidth': 0.6, 'size': 6, 'alpha': 0.7,
                            'palette': pal, 'hue_order': hue_order}

        cid = str(self.cid)

        ax = sns.boxplot(y=data[cid], x='Dataset',hue = 'Condition', data = data, fliersize=0, **boxplot_kwargs)
        ax = sns.stripplot(y=data[cid], x='Dataset',hue = 'Condition', data = data, split=True, jitter=0.2, **stripplot_kwargs)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles[0:2], labels[0:2], loc='best', fontsize='medium',
               handletextpad=0.5)
        lgd.legendHandles[0]._sizes = [40]
        lgd.legendHandles[1]._sizes = [40]
        ax.get_legend().remove()
        plt.ylabel("Peak Intensity (log2)")
        title = 'Peak ID '+str(self.cid) + ':m/z=' + str(round(self.mz, 4))
        plt.title(title)



        #plt.show()

        fig.savefig(path+str(cid)+'.png', format='png', dpi=400, bbox_inches='tight')
        plt.close(fig)


    def __str__(self):
        return """Peak {self.cid}: m/z = {self.mz}, rt = {self.rt} \n
                p-val: {self.pval}, logFC = {self.logfc} \n
                mm annotation: {self.mm_annotation} \n
                std annotation: {self.std_annotation} \n
                hmdb annotation: {self.ms2_annotation} \n

                {self.spectra} \n""".format(self = self)




#make a peakinfolist class
def get_peak_by_cid(peakinfolist, cid):
    for peak in peakinfolist:
        if peak.cid == cid:
            return peak
