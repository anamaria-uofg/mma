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

    def __init__(self, mz, ms2specdata, metabolites_dict, ms2_files):

        self.adducts = []
        adduct_list = self.compute_adducts(mz)
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
                 (mz + 83.9613 - PROTON, 'M-HCOOK+H[1+]', 'HCO2K'),
                 ]

        return addList


class PeakInfo:

    def __init__(self, cid, mz, rt):
        self.cid = cid
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

    def add_std_match(self, standards_dict, metabolites_dict, ppm = 0.1, ppm_rt = 60):

        for key in standards_dict.keys():
            if self.mz >= key-ppm and self.mz <= key + ppm:
                if self.rt*60 >= standards_dict[key][1]-ppm_rt and self.rt*60 <= standards_dict[key][1]+ppm_rt:
                    self.std_annotation = standards_dict[key][0]

        for metabolite in metabolites_dict:
            if str(self.std_annotation).strip() == str(metabolites_dict[metabolite][2]).lower().strip():
                self.std_kegg_id = metabolites_dict[metabolite][5]

    def add_spectra(self, aligner_gp):
        new_peakid = self.cid - 1
        num_peaks = aligner_gp.peaksets[new_peakid].n_peaks
        spectra = {}
        if num_peaks > 1:
            for i in range(num_peaks-1):
                    source = aligner_gp.peaksets[new_peakid].peaks[i+1].source_file
                    msms = aligner_gp.peaksets[new_peakid].peaks[i+1].ms2_spectrum
                    spectra[source] = msms

        self.spectra = spectra

    def add_adducts(self, metabolites_dict, ms2_files):
        self.adducts = AdductList(self.mz, self.spectra, metabolites_dict, ms2_files)
        self.best_ms2_match_adduct = self.adducts.get_adduct_with_best_score()

    def add_ms2_match(self, metabolites_dict, score_theshold = 0.3):
        if self.best_ms2_match_adduct:
            if self.best_ms2_match_adduct.score > score_theshold:
                self.ms2_annotation = metabolites_dict[self.best_ms2_match_adduct.hmdb][2]
                self.ms2_kegg_id = metabolites_dict[self.best_ms2_match_adduct.hmdb][5]

    def get_possible_kegg_ids(self):
        kegg_ids = []
        if self.mm_kegg_id != '':
            kegg_list = self.mm_kegg_id.split(';')
            for i in kegg_list:
                kegg_ids.append(i)

        if self.std_kegg_id != '':
            kegg_ids.append(self.std_kegg_id)

        if self.ms2_kegg_id != '':
            kegg_ids.append(self.ms2_kegg_id)

        return kegg_ids

    def plot_boxplots(self, data, order = ['controlVL','infectedVL','controlMalaria', 'infectedMalaria', 'controlZika','infectedZika']):
        import seaborn as sns
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,5))

        ax = sns.boxplot(y=data[str(self.cid)], x='ConditionDataset', data = data, order = order, palette = 'Set2' )
        ax = sns.swarmplot(y=data[str(self.cid)], x='ConditionDataset', data = data, order = order, color="black")

        plt.ylabel("Intensity (log2)")
        plt.title(self.mz)
        plt.show()

    def __str__(self):
        return """peak {self.cid}: {self.mz}, {self.rt} \n
                p-val: {self.pval}, {self.tval} \n
                {self.mm_annotation} \n
                {self.mm_pathway} \n
                {self.mm_kegg_id} \n
                {self.spectra} \n
                {self.adducts} \n
                {self.best_ms2_match_adduct}""".format(self = self)

    def get_peak_by_id(self, cid):
        return self



def get_peak_by_cid(peakinfolist, cid):
    for peak in peakinfolist:
        if peak.cid == cid:
            return peak
