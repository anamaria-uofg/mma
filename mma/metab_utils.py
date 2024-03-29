###methods useful for analysis of metabolomics files
import pandas as pd
import numpy as np
from scipy import stats
import csv
import matplotlib.pyplot as plt
import glob
import os
import re
import itertools

import sys
import mass_spec_utils


WIN = 'D:/'
MAC = '/Volumes/Transcend2/17_20_PhD/19_20_PhD_Metabolomics/'

osp = MAC
#MNET_PATH = '/Users/anamaria/git/molnet/code/'
#sys.path.append(MNET_PATH)

sys.path.append('/Users/anamaria/git/pymzmine/')
from mzmine import align, load_aligned_peaks, match_aligned_to_original
MZMINE_COMMAND = '/Users/anamaria/git/MZmine-2.40.1/startMZmine_MacOSX.command'

def calculate_tolerance(x, ppm):
    tolerance = x*ppm/1000000
    return tolerance

def create_std_dict(input_std_csv_file, polarity, include_actual_rt= True):

    """ Creates and returns a dictionary for a standards csv file obtained from ToxID."""

    standards = {}
    for _,row in input_std_csv_file.iterrows():
        compound_name = row['Compound Name'].lower()
        detected_mz = row['Detected m/z']
        expected_rt = row['Expected RT']
        actual_rt = row['Actual RT']
        cmpd_polarity = row['Polarity']
        if cmpd_polarity == polarity:
            if include_actual_rt:
                standards[compound_name] = (detected_mz, expected_rt, actual_rt, polarity)
            else:
                standards[compound_name] = (detected_mz, expected_rt, polarity)
    return standards

def create_std_dict_by_mz(input_std_csv_file, polarity):

    """ Creates and returns a dictionary for a standards csv file obtained from ToxID with m/z as the key."""
    standards = {}
    for _,row in input_std_csv_file.iterrows():
        compound_name = row['Compound Name'].lower()
        detected_mz = row['Detected m/z']
        expected_rt = row['Expected RT']
        actual_rt = row['Actual RT']
        polar = row['Polarity']
        adduct_int = row['Adduct1 H+']
        if polarity == 'both':
            standards[detected_mz] = (compound_name, expected_rt*60, actual_rt*60, polar, adduct_int)
        elif polarity == '+':
            if polar == '+':
                standards[detected_mz] = (compound_name, expected_rt*60, polar, adduct_int)
        elif polarity == '-':
            if polar == '-':
                standards[detected_mz] = (compound_name, expected_rt*60, polar, adduct_int)
    return standards

def createMetaboliteSetList(original_files, output_dir, stds_csvs, original_csvs, rt_range = 0.5, mz_range = 0.0003):

    from standardsMetabolite import MetaboliteSet, Metabolite
    listMetabolites = []

    for file_pos,o in enumerate(original_files):

        with open(stds_csvs[file_pos],'r') as s:
             for line in s:
                #to skip the comments read only the rows which start with a number
                if re.match(r"^\d+.*$",line):
                    name = line.split(',')[2].lower()
                    polarity = line.split(',')[4]
                    mz = float(line.split(',')[6])
                    rt = float(line.split(',')[9])
                    metaboliteSet = MetaboliteSet(name)

                    with open(original_csvs[file_pos], 'r') as c:
                        reader = csv.reader(c)
                        header = next(reader)
                        for line in reader:

                            id_o = int(line[0])
                            mz_o = float(line[1])
                            rt_o = float(line[2])
                            int_o = float(line[3])

                            if (mz_o <= mz+mz_range and mz_o >= mz-mz_range):
                                if (rt_o <= rt+rt_range and rt_o >= rt-rt_range):
                                    metabolite = Metabolite(o,id_o,mz_o,rt_o,int_o)
                                    metaboliteSet.add_file(o)
                                    metaboliteSet.add_metabolite(metabolite)
                    listMetabolites.append(metaboliteSet)
    return listMetabolites

def get_std_matches_for_files(original_files, output_dir, stds_csvs, original_csvs, rt_range = 0.5, mz_ppm = 3):

    matches = {}
    for file_pos,o in enumerate(original_files):


        with open(stds_csvs[file_pos],'r') as s:
             for line in s:
                #to skip the comments read only the rows which start with a number
                if re.match(r"^\d+.*$",line):
                    name = line.split(',')[2].lower()
                    polarity = line.split(',')[4]
                    if line.split(',')[6] != '':
                        mz = float(line.split(',')[6])
                    if line.split(',')[9] != '':
                        rt = float(line.split(',')[9])

                    with open(original_csvs[file_pos], 'r') as c:
                        reader = csv.reader(c)
                        header = next(reader)
                        for line in reader:

                            id_o = int(line[0])
                            mz_o = float(line[1])
                            rt_o = float(line[2])
                            int_o = float(line[3])

                            if not name in matches:
                                matches[name] = {}
                            mz_range = calculate_tolerance(mz,3)
                            if (mz-mz_range <= mz_o <= mz+mz_range):
                                if (rt-rt_range <= rt_o <= rt+rt_range):
                                    matches[name][o] = (id_o,mz_o,rt_o,int_o)


    output_file = os.path.join(output_dir,'stds_match_links.csv')
    with open(output_file,'w') as f:
        heads = ['stds_name'] + original_files
        writer = csv.writer(f)
        writer.writerow(heads)
        for matched_peak in matches:

            new_row = [matched_peak]

            for o in original_files:

                val = matches[matched_peak].get(o,None)

                if val:
                    new_row.append(val[0])
                else:
                    new_row.append('null')
            writer.writerow(new_row)

    return matches

def get_rts_between_datasets(dataset1, dataset2, matches, replicate = 1):
    if replicate == 1:
        replicate = '_1_'
    elif replicate == 2:
        replicate = '_2_'
    dataset1_rts = []
    dataset2_rts = []
    dataset1_mzs = []
    dataset2_mzs = []
    rts_diff = []
    for metabolite in matches:
        count = 0
        for file in matches[metabolite]:

            if replicate in file:
                if dataset1 in file.lower():
                    dataset1_rt = matches[metabolite][file][2]
                    dataset1_mz = matches[metabolite][file][1]
                    count += 1
                if dataset2 in file.lower():
                    dataset2_rt = matches[metabolite][file][2]
                    dataset2_mz = matches[metabolite][file][1]
                    count += 1
        if count == 2:
            dataset1_rts.append(dataset1_rt)
            dataset2_rts.append(dataset2_rt)
            dataset1_mzs.append(dataset1_mz)
            dataset2_mzs.append(dataset1_mz)
            rts_diff.append(dataset1_rt-dataset2_rt)
    return dataset1_rts, dataset1_mzs, dataset2_rts, dataset2_mzs

def get_stds_between_datasets(dataset1, dataset2, matches, replicate = 1):
    if replicate == 1:
        replicate = '_1_'
    elif replicate == 2:
        replicate = '_2_'
    dataset1_info = []
    dataset2_info = []

    for metabolite in matches:
        count = 0
        for file in matches[metabolite]:

            if replicate in file:
                if dataset1 in file.lower():
                    dataset1_in = (matches[metabolite][file][1], matches[metabolite][file][2], matches[metabolite][file][3])
                    count += 1
                if dataset2 in file.lower():
                    dataset2_in = (matches[metabolite][file][1], matches[metabolite][file][2], matches[metabolite][file][3])
                    count += 1
        if count == 2:
            dataset1_info.append(dataset1_in)
            dataset2_info.append(dataset2_in)

    return dataset1_info, dataset2_info

def get_rts_within_datasets(dataset, matches):

    replicate1_rts = []
    replicate2_rts = []
    rts_diff = []
    for metabolite in matches:
        count = 0

        for file in matches[metabolite]:

            if dataset in file.lower():

                if '_1_'  in file:

                    replicate1_rt = matches[metabolite][file][2]
                    count = count + 1
                if '_2_'  in file:

                    replicate2_rt = matches[metabolite][file][2]
                    count = count + 1

        if count == 2:
            replicate1_rts.append(replicate1_rt)
            replicate2_rts.append(replicate2_rt)
            rts_diff.append(replicate1_rt-replicate2_rt)
    return replicate1_rts, replicate2_rts, rts_diff

def get_rt_diff_between_datasets(dataset1_mzmine_stds, dataset2_mzmine_stds, file_number = 3):

    dataset1_rt = []
    dataset2_rt = []
    stds = {}
    diff = []

    for i in range(file_number):
        dataset1 = set(dataset1_mzmine_stds[i])
        dataset2 = set(dataset2_mzmine_stds[i])

        for name in dataset1.intersection(dataset2):

            diff.append((dataset1_mzmine_stds[i][name][2]-dataset2_mzmine_stds[i][name][2])*60)
            dataset1_rt.append(dataset1_mzmine_stds[i][name][2])
            dataset2_rt.append(dataset2_mzmine_stds[i][name][2])

            stds[name] = (i+1, dataset1_mzmine_stds[i][name], dataset2_mzmine_stds[i][name])



    return (dataset1_rt, dataset2_rt, stds, diff)

def get_stats_on_diff(dataset1_rt,name1, dataset2_rt,name2,diff, setbinwidth, zcut = 1):


    zscore = np.abs(stats.zscore(diff))
    d1_rt_mod = []
    d2_rt_mod = []

    for (d1, d2, z) in zip(dataset1_rt, dataset2_rt, zscore):
         if (z<zcut):
                d1_rt_mod.append(d1)
                d2_rt_mod.append(d2)

    d1_d2_diff_no_outliers = (np.subtract([x * 60 for x in d1_rt_mod],[x*60 for x in d2_rt_mod]))

    d1_d2_diff_std = np.std(d1_d2_diff_no_outliers)
    d1_d2_diff_mean = np.mean(d1_d2_diff_no_outliers)
    d1_d2_diff_max = np.max(abs(d1_d2_diff_no_outliers))

    print('Mean: {0:.2f} \nSD: {1:.2f}\nMaximum difference: {2:.2f}'.format(d1_d2_diff_mean, d1_d2_diff_std,  d1_d2_diff_max))

    binwidth = setbinwidth
    plt.hist(d1_d2_diff_no_outliers, bins=np.arange(min(d1_d2_diff_no_outliers), max(d1_d2_diff_no_outliers) + binwidth, binwidth), alpha =0.5)
    plt.xlabel('RT '+ name1+ ' - RT ' +name2 + '(s)')
    plt.show()

def try_gp_regressions(t_dataset_mod, t_reference_mod, plot= False):
    import GPy
    from sklearn import metrics

    k_lin = GPy.kern.Linear(1)
    k_exp = GPy.kern.Exponential(1)
    k_nn = GPy.kern.MLP(1)

    k_rbf = GPy.kern.RBF(input_dim=1, variance=0.1, lengthscale=0.5)
    k_mat32 = GPy.kern.Matern32(input_dim=1, variance=2., lengthscale=0.2)
    k_mat52 = GPy.kern.Matern52(1)

    k_per = GPy.kern.StdPeriodic(1, period=3.)
    k_cos = GPy.kern.Cosine(1)
    k_brwn = GPy.kern.Brownian(1)

    t_reference = np.array(t_reference_mod)
    t_dataset = np.array(t_dataset_mod)

    REF = t_reference[:,None]
    TOPRED = t_dataset[:,None]


    #Cross-validation
    from sklearn.model_selection import train_test_split

    training_data, test_data, training_data_target, test_data_target = train_test_split(TOPRED, REF, test_size=0.3, random_state=0)

    #training_data = TOPRED[:len(TOPRED)//2] #Need to put double / so it could work
    #training_data_target = REF[:len(REF)//2]
    #test_data = TOPRED[len(TOPRED)//2+1:]
    #test_data_target = REF[len(REF)//2+1:]

    #write a vector with kernels which might work
    ks = [k_nn, k_nn+k_rbf,k_nn*k_rbf]
    ks_names = ['MLP','MLP+RBF','MLP*RBF']
    names = ['RBF','MLP','MLP+RBF','MLP*RBF']

    accuracy_list = []
    mae_list = []
    mse_list = []

    mtest = GPy.models.GPRegression(training_data,training_data_target,k_rbf)
    mtest.optimize()
    if plot:
        mtest.plot()
        plt.title('RBF')
        plt.show()

    gpr_predicted_data,_ = mtest.predict(test_data)
    if plot:
        plt.plot(test_data_target, gpr_predicted_data,  '.')
        plt.plot(test_data_target, test_data_target, '.r')
        plt.ylabel('Predicted data')
        plt.xlabel('True data')
        plt.show()

    accuracy0 = metrics.r2_score(test_data_target, gpr_predicted_data)
    mae0 = metrics.mean_absolute_error(test_data_target, gpr_predicted_data)
    mse0 = metrics.mean_squared_error(test_data_target, gpr_predicted_data)

    accuracy_list.append(accuracy0)
    mae_list.append(mae0)
    mse_list.append(mse0)

    print('Cross-Predicted Accuracy for','RBF',':', accuracy0)
    print('Mean absolute error for','RBF',':', mae0)
    print('Mean squared error for','RBF',':', mse0)

    mfinal = mtest
    kfinal = k_rbf
    namefinal = 'RBF'

    for k, name in zip(ks,ks_names):
        #regression and optimisation of model
        mtest = GPy.models.GPRegression(training_data,training_data_target,k)
        mtest.optimize()
        if plot:
            mtest.plot()
            plt.title(name)
            plt.show()

        gpr_predicted_data,_ = mtest.predict(test_data)
        if plot:
            plt.plot(test_data_target, gpr_predicted_data,  '.')
            plt.plot(test_data_target, test_data_target, '.r')
            plt.ylabel('Predicted data')
            plt.xlabel('True data')
            plt.show()

        accuracy = metrics.r2_score(test_data_target, gpr_predicted_data)
        mae = metrics.mean_absolute_error(test_data_target, gpr_predicted_data)
        mse = metrics.mean_squared_error(test_data_target, gpr_predicted_data)

        accuracy_list.append(accuracy)
        mae_list.append(mae)
        mse_list.append(mse)

        print('Cross-Predicted Accuracy for',name,':', accuracy)
        print('Mean absolute error for',name,':', mae)
        print('Mean squared error for',name,':', mse)

        if (accuracy >= accuracy0 and mae < mae0 and mse < mse0):
            accuracy0 = accuracy
            mae0 = mae
            mse0 = mse
            mfinal = mtest
            kfinal = k
            namefinal = name
    results_table = pd.DataFrame(data=[names, accuracy_list, mae_list, mse_list])
    print("Final kernel:",namefinal)
    return mfinal, kfinal, results_table

def return_data_with_no_outliers(dataset1, dataset2,diff, zscore_cutoff):
    #Eliminating outliers using zscore
    zscore = np.abs(stats.zscore(diff))
    mod1 = []
    mod2 = []

    for (d1, d2, z) in zip(dataset1, dataset2, zscore):
         if (z<zscore_cutoff):
            mod1.append(d1)
            mod2.append(d2)

    return np.array(mod1), np.array(mod2)

def return_data_with_no_outliers_2(dataset1, dataset2, zscore_cutoff):
    #Eliminating outliers using zscore
    from scipy import stats
    diff = np.subtract(np.array([i[1] for i in dataset1]), np.array([i[1] for i in dataset2]))
    zscore = np.abs(stats.zscore(diff))
    mod1 = []
    mod2 = []

    for (d1, d2, z) in zip(dataset1, dataset2, zscore):
         if (z<zscore_cutoff):
            mod1.append(d1)
            mod2.append(d2)

    return mod1, mod2

def create_bins_mz_range(dataset, n_bins):
    dataset.sort()
    lowest_value = np.min([i[0] for i in dataset])
    highest_value = np.max([i[0] for i in dataset])
    width = (highest_value - lowest_value)/n_bins
    bins = []
    for i in range(n_bins):
        bins.append((lowest_value+width*i, lowest_value+width*(i+1)))
    return bins

def create_bins_same_width(dataset, n_bins):
    total = len(dataset)
    bin_length = total / n_bins

    return int(bin_length)

def plot_gpr_bins_mz_range(n_bins,dataset1_string, dataset2_string, matches, optimization_restarts, zscore):

    info1, info2 = get_stds_between_datasets(dataset1_string, dataset2_string, matches)
    info1.sort()
    info2.sort()

    info1, info2 = return_data_with_no_outliers_2(info1, info2, zscore)

    bins = create_bins_mz_range(info1, n_bins)

    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    fig, ax = plt.subplots(figsize=(20,10))

    for i in range(n_bins):
        RT1 = []
        RT2 = []
        for inf1, inf2 in zip(info1,info2):
            if inf1[0] >= bins[i][0] and inf1[0] <= bins[i][1]:
                RT1.append(inf1[1])
                RT2.append(inf2[1])
        print("Bin:", i, "\nm/z range:", bins[i])
        kernel = GPy.kern.RBF(input_dim=1, variance = 0.5)
        model = GPy.models.GPRegression(np.array(RT1)[:,None], np.array(RT2)[:,None] - np.array(RT1)[:,None],kernel)
        if optimization_restarts:
            model.optimize_restarts(10)
        else:
            model.optimize()

        print(model.parameters[0])
        model.plot_mean(color = colors[i], ax=ax, label = ('Mean for Bin',i, bins[i]), plot_limits = (2.5, 22.5))
        model.plot_confidence(color = colors[i], ax=ax, label = ('Confidence for Bin',i), plot_limits = (2.5, 22.5))
        model.plot_data(ax=ax)


    plt.xlabel("RT(min)")
    plt.ylabel("RT drift")
    plt.xlim((2.5,22.5))
    plt.show()

def plot_gpr_bins_same_width(n_bins, dataset1_string, dataset2_string, matches, optimization_restarts, zscore):

    info1, info2 = get_stds_between_datasets(dataset1_string, dataset2_string, matches)
    info1.sort()
    info2.sort()

    info1, info2 = return_data_with_no_outliers_2(info1, info2, zscore)

    bins_length = create_bins_same_width(info1, n_bins)

    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9','C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9' ]
    fig, ax = plt.subplots(figsize=(20,10))

    for f in range(n_bins):

        inf1 =  info1[f*bins_length : (f+1)*bins_length]
        inf2 = info2[f*bins_length : (f+1)*bins_length]

        #inf1, inf2 = return_data_with_no_outliers(inf1, inf2, zscore)

        print("Bin:", f, "\nm/z range:",np.min([i[0] for i in inf1]),",", np.max([i[0] for i in inf1]))

        kernel = GPy.kern.RBF(input_dim=1)
        model = GPy.models.GPRegression(np.array([i[1] for i in inf1])[:,None], np.array([i[1] for i in inf2])[:,None] - np.array([i[1] for i in inf1])[:,None],kernel)
        if optimization_restarts:
            model.optimize_restarts(10)
        else:
            model.optimize()

        print(model.parameters[0])
        model.plot_mean(color = colors[f], ax=ax, label = ('Mean for Bin',f,np.min([i[0] for i in inf1]), np.max([i[0] for i in inf1])), plot_limits = (2.5, 22.5))
        model.plot_confidence(color = colors[f], ax=ax, label = ('Confidence for Bin',f), plot_limits = (2.5, 22.5))
        model.plot_data(ax=ax)

        #model.plot()

    plt.xlabel("RT(min)")
    plt.ylabel("RT drift")
    plt.xlim((2.5,22.5))
    plt.show()

def get_non_anomalies(dataset_1_rt, dataset_2_rt):
    from sklearn.ensemble import IsolationForest

    data = np.array([dataset_1_rt, dataset_2_rt]).transpose()
    df = pd.DataFrame(data)
    l = len(dataset_1_rt)

    clf = IsolationForest(behaviour = 'auto', max_samples=86, random_state = 1, contamination= 'auto')
    pred = clf.fit_predict(data)
    df['pred'] = pred
    df_no_anomalies = df[df['pred'] == 1]
    dataset_1_rt_mod = np.array(df_no_anomalies[0])
    dataset_2_rt_mod = np.array(df_no_anomalies[1])
    return dataset_1_rt_mod, dataset_2_rt_mod

def plot_and_get_values_for_regression(reference_dataset,refname, compared_dataset,compname, diff, ml = True, zscore = 1.4):

    plt.plot(np.array(compared_dataset), np.array(reference_dataset) - np.array(compared_dataset), '.')
    plt.xlabel('RT '+refname+' (min)')
    plt.ylabel('RT '+compname+' - RT '+refname+' (min)')

    plt.show()


    if ml:
        compared_dataset_mod, reference_dataset_mod = get_non_anomalies(compared_dataset, reference_dataset)
    else:
        compared_dataset_mod, reference_dataset_mod = return_data_with_no_outliers(compared_dataset, reference_dataset, diff, zscore)

    plt.plot(compared_dataset_mod, reference_dataset_mod - compared_dataset_mod, '.')
    plt.plot(reference_dataset_mod[:,None], reference_dataset_mod[:,None] - reference_dataset_mod[:,None], 'r')
    plt.xlabel('RT '+refname+' (min)')
    plt.ylabel('RT '+compname+' - RT '+refname+' (min)')
    plt.title('Before correction')
    plt.ylim(-4,4)
    return reference_dataset_mod, compared_dataset_mod

def modify_rtdrift_in_mztab(output_dir, output_dir_new, model, name_tag, standards = False):

    if standards == False:
        filelen = len(glob.glob(os.path.join(output_dir,name_tag)))
    else:
        filelen = 6
    print(filelen)

    peak_file_list = glob.glob(os.path.join(output_dir, name_tag))
    peak_file_list.sort()
    peak_file_list_new = glob.glob(os.path.join(output_dir_new, name_tag))
    peak_file_list_new.sort()

    for i in range(filelen):

        with open(peak_file_list[i], 'r') as f:
            with open(peak_file_list_new[i],'w') as g:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('SMH'):
                        tokens = line.split()
                        rt_pos = [tokens.index('retention_time')]
                        rt_pos.append(tokens.index('opt_assay[1]_peak_rt'))
                    if not line.startswith('SML'):
                        g.write(line + '\n')
                    else:
                        tokens = line.split()
                        old_rt = float(tokens[rt_pos[0]])

                        new_rt,_ = model.predict(np.array([[old_rt]]))
                        print(new_rt[0][0])
                        new_rt = new_rt[0][0]
                        print(old_rt)
                        new_rt += old_rt
                        print(new_rt)

                        tokens[rt_pos[0]] = (new_rt)
                        tokens[rt_pos[1]] = (new_rt)
                        new_line = '\t'.join([str(t) for t in tokens])
                        g.write(new_line + '\n')

def modify_rt_in_mztab(output_dir, output_dir_new, model, name_tag = '*.mzTab', standards = False):

    if standards == False:
        filelen = len(glob.glob(os.path.join(output_dir,name_tag)))
    else:
        filelen = 6

    for i in range(filelen):
        with open(glob.glob(os.path.join(output_dir,name_tag))[i], 'r') as f:
            with open(glob.glob(os.path.join(output_dir_new,name_tag))[i],'w') as g:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('SMH'):
                        tokens = line.split()
                        rt_pos = [tokens.index('retention_time')]
                        rt_pos.append(tokens.index('opt_assay[1]_peak_rt'))
                    if not line.startswith('SML'):
                        g.write(line + '\n')
                    else:
                        tokens = line.split()
                        old_rt = float(tokens[rt_pos[0]])
                        print(old_rt)
                        #new_rt = brr.predict(np.array(old_rt).reshape(-1,1))[0]
                        new_rt, = model.predict(np.array(old_rt).reshape(-1,1))[0]
                        print(new_rt[0])

                        tokens[rt_pos[0]] = new_rt[0]
                        tokens[rt_pos[1]] = new_rt[0]
                        new_line = '\t'.join([str(t) for t in tokens])
                        g.write(new_line + '\n')

def modify_rtdrift_in_csv(csv_file, csv_file_new, model):
    with open(csv_file, 'r') as f:
        with open(csv_file_new, 'w') as g:
            for line in f:
                line = line.rstrip()
                if line.startswith('row'):
                    g.write(line+'\n')
                else:
                    tokens = line.split(',')
                    old_rt = float(tokens[2])
                    new_rt,_ = model.predict(np.array([[old_rt]]))
                    new_rt = new_rt[0][0]
                    new_rt += old_rt
                    tokens[2] = new_rt

                    new_line = ','.join([str(t) for t in tokens])
                    g.write(new_line + '\n')

def get_predicted_rt(old_rt, model):
    new_rt,_ = model.predict(np.array([[old_rt]]))
    new_rt = new_rt[0][0]
    new_rt += old_rt
    return new_rt

def modify_rtdrift_in_csv_with_boxes(csv_file, csv_file_new, model):
    with open(csv_file, 'r') as f:
        with open(csv_file_new, 'w') as g:
            for line in f:
                line = line.rstrip()
                if line.startswith('row'):
                    g.write(line+'\n')
                else:
                    tokens = line.split(',')

                    old_rt = float(tokens[2])
                    old_rt_min = float(tokens[4])
                    old_rt_max = float(tokens[5])


                    tokens[2] = get_predicted_rt(old_rt, model)
                    tokens[4] = get_predicted_rt(old_rt_min, model)
                    tokens[5] = get_predicted_rt(old_rt_max, model)



                    new_line = ','.join([str(t) for t in tokens])
                    g.write(new_line + '\n')

def change_rt_in_mgf(mgf_file, model):
    #changed it to work with mass-spec-utils
    for idspec in mgf_file:
        spectrum = mgf_file[idspec]
        old_rt = float(spectrum.metadata['RTINSECONDS'])/60
        new_rt,_ = model.predict(np.array([[old_rt]]))
        new_rt = new_rt[0][0]
        new_rt += old_rt
        new_rt = new_rt*60
        spectrum.metadata['RTINSECONDS'] = new_rt

def modify_rtdrif_in_mzxml(file_path, file_name, model):
    import xml.etree.ElementTree
    file = file_path + file_name
    et = xml.etree.ElementTree
    et.register_namespace('',"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2")
    et.register_namespace('xsi', "http://www.w3.org/2001/XMLSchema-instance")


    tree = et.parse(file)
    root = tree.getroot()

    tag = '{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}'
    for child in root:
        if child.tag == tag+'msRun':
        ##change rt here too. no
            for element in child:

                if element.tag == tag+'scan':


                    old_rt = element.attrib['retentionTime'].split('PT')[1].split('S')[0] ##rt in s

                    old_rt = float(old_rt)/60
                    new_rt,_ = model.predict(np.array([[old_rt]]))
                    new_rt = new_rt[0][0]
                    new_rt += old_rt
                    new_rt = new_rt*60

                    if new_rt < 0: ##some values are negative after gp modification
                        new_rt = old_rt * 60 ##so keep the old values for those (anyway we're not interested in metbolites eluted before first 120 s

                    new_rt = f"{new_rt:.5f}"

                    element.attrib['retentionTime'] = 'PT'+new_rt+'S'

    newfile = file_path + file_name.split('.mzXML')[0] + '_modified.mzXML'

    tree.write(newfile)

def get_good_peaks(file, matches, stds_matches):
    good_peaks = 0
    for match in matches:
        for name in stds_matches:
            if check_peak_is_good(name, matches[match], file, stds_matches):
                good_peaks+=1
    return good_peaks

def check_peak_is_good(name, peak, files, stds_matches):

    for file in files:
        if file in peak:
            peak_id = peak[file][0]
            if file in stds_matches[name]:
                peak_id2 = stds_matches[name][file][0]
                if (peak_id != peak_id2):
                    return False
            else:
                return False
        else:
            return False

    return True

def get_total_stds(file, stds_matches):
    total = 0
    for name in stds_matches:
        if file in stds_matches[name]:
            total+=1

    return total

def get_peaks_for_files(output_dir, files, rt_window, xml_template, stds_matches):

    print(files)
    dict_peaks = {}
    dict_peaks[str(files)] = []



    for rt in rt_window:
        print(rt)
        align(output_dir, xml_template = xml_template, mzmine_command = MZMINE_COMMAND, rt_window = rt)
        aligned_peaks,original_file,f_idx_dict = load_aligned_peaks(output_dir)
        #print(f_idx_dict)
        matches = match_aligned_to_original(aligned_peaks,files,output_dir,f_idx_dict,write_file = True)
        good_peaks = get_good_peaks(files, matches, stds_matches)

        dict_peaks[str(files)].append((len(aligned_peaks),good_peaks))


    return dict_peaks


def extract_peaksets_with_spectra_from_filt_dataset(aligner, dataid):
    peaksets = []
    for peakid in list(dataid):
        new_peakid = peakid - 1
        peaksets.append(aligner.peaksets[new_peakid])

    return peaksets

def extract_scores_pairs_for_dataset(peaksets, dataset1, dataset2):
    from scoring_functions import fast_cosine
    good_scores = {}
    for peakset in peaksets:
        actual_scores = []
        np = peakset.n_peaks
        actual_mz = peakset.mean_mz
        actual_rt = peakset.mean_rt
        for i in range(np-2):
            if (dataset1 in peakset.peaks[i+1].ms2_spectrum.file_name) or (dataset1 in peakset.peaks[i+2].ms2_spectrum.file_name):
                if (dataset2 in peakset.peaks[i+1].ms2_spectrum.file_name) or (dataset2 in peakset.peaks[i+2].ms2_spectrum.file_name):
                    score, used_matches = fast_cosine(peakset.peaks[i+1].ms2_spectrum, peakset.peaks[i+2].ms2_spectrum, 0.2, 2)
                    actual_scores.append(score)
        if len(actual_scores)>0:
            good_scores[(actual_mz, actual_rt)] = actual_scores
    return good_scores

def plot_boxplots(data, metab_id, name_column, ylim, y = False ):
    import seaborn as sns
    plt.figure(figsize=(10,5))

    ax = sns.boxplot(y=data[metab_id], x=name_column, data = data, order = ['controlVL','infectedVL','controlMalaria', 'infectedMalaria', 'controlZika','infectedZika'], palette = 'Pastel2' )
    ax = sns.swarmplot(y=data[metab_id], x=name_column, data = data, order = ['controlVL','infectedVL','controlMalaria', 'infectedMalaria', 'controlZika','infectedZika'], color="black")
    if y:
        plt.ylim(ylim)

    plt.ylabel("Intensity (log2)")
    plt.title(metab_id)
    plt.show()

def plot_boxplots_individual_dataset(data, metab_id, name_column, ylim = (10,30), y = False ):
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5,5))

    ax = sns.boxplot(y=data[metab_id], x=name_column, data = data, order = ['controlHATcontrol','infectedHATS1', 'infectedHATS2'])
    ax = sns.swarmplot(y=data[metab_id], x=name_column, data = data, order = ['controlHATcontrol','infectedHATS1', 'infectedHATS2'], color="black")
    if y:
        plt.ylim(ylim)
    plt.show()

def mergeDict(dict1, dict2):
    ''' Merge dictionaries and keep values of common keys in list'''
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = [value , dict1[key]]

    return dict3

def plot_spectra(score, rt_file1, rt_file2, spectrum_file1, spectrum_file2):
    print('Score:', score, 'RT1:', rt_file1, 'RT2:', rt_file2)
    spectrum_file1.plot()
    spectrum_file2.plot()
    plt.show()

def get_bad_spectral_matches(mgf_file1, mgf_file2, rtdiff, tolerance, minmatch, plot = True):
    from scoring_functions import fast_cosine
    tolerance = tolerance
    minmatch = minmatch
    badscores = []

    for spectrum_id1 in mgf_file1:
        spectrum_file1 = mgf_file1[spectrum_id1]
        mz_file1 = spectrum_file1.precursor_mz
        rt_file1 = spectrum_file1.rt

        if rt_file1 >= 120:

            for spectrum_id2 in mgf_file2:
                spectrum_file2 = mgf_file2[spectrum_id2]
                mz_file2 = spectrum_file2.precursor_mz
                rt_file2 = spectrum_file2.rt

                if rt_file2 >= 120:

                    if mz_file2 >= (mz_file1 - 0.01) and mz_file2 <= (mz_file1 + 0.01):

                        if rt_file2 >= (rt_file1 + rtdiff) or rt_file2 <= (rt_file1 - rtdiff):

                            score, used_matches = fast_cosine(spectrum_file1, spectrum_file2, tolerance, minmatch)

                            if plot:
                                if score > 0.6:
                                    plot_spectra(score, rt_file1, rt_file2, spectrum_file1, spectrum_file2)

                            badscores.append(score)
    return badscores

def remove_small_peaks(mgf_file, percentage):

    for spectrum_id in mgf_file:
        spectrum = mgf_file[spectrum_id]
        spectrum_max_int = max([intensity for mz,intensity in spectrum.peaks])

        threshold = percentage/100 * spectrum_max_int
        #threshold = number

        new_peaks = []


        for mz, intensity in spectrum.peaks:

            if intensity > threshold:
                new_peaks.append((mz, intensity))

        spectrum.peaks = new_peaks
        spectrum.n_peaks = len(spectrum.peaks)
        spectrum.normalised_peaks = sqrt_normalise(spectrum.peaks)

def sqrt_normalise(peaks):
    import math
    temp = []
    total = 0.0
    for mz,intensity in peaks:
        temp.append((mz,math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz,intensity in temp:
         normalised_peaks.append((mz,intensity/norm_facc))
    return normalised_peaks

def extract_peaksets_with_spectra_from_filt_dataset(aligner, dataid):
    peaksets = []
    for peakid in list(dataid):
        new_peakid = peakid - 1
        peaksets.append(aligner.peaksets[new_peakid])
    return peaksets

def extract_peaksets_with_spectra_from_dataset(aligner):
    peaksets = []
    for peakset in aligner.peaksets:
        if peakset.n_peaks > 2:
            peaksets.append(peakset)
    return peaksets

def extract_peaksets_with_spectra(peaksets):
    epeaksets = []
    for peakset in peaksets:
        np = peakset.n_peaks
        if np>2:
            epeaksets.append(peakset)
    return epeaksets

def extract_scores_pairs_for_dataset(peaksets, dataset1, dataset2):
    from scoring_functions import fast_cosine
    good_scores = {}
    for peakset in peaksets:
        actual_scores = []
        np = peakset.n_peaks
        actual_mz = peakset.mean_mz
        actual_rt = peakset.mean_rt
        for i in range(np-2):
            if (dataset1 in peakset.peaks[i+1].ms2_spectrum.file_name) or (dataset1 in peakset.peaks[i+2].ms2_spectrum.file_name):
                if (dataset2 in peakset.peaks[i+1].ms2_spectrum.file_name) or (dataset2 in peakset.peaks[i+2].ms2_spectrum.file_name):
                    score, used_matches = fast_cosine(peakset.peaks[i+1].ms2_spectrum, peakset.peaks[i+2].ms2_spectrum, 0.2, 2)
                    actual_scores.append(score)
        if len(actual_scores)>0:
            good_scores[(actual_mz, actual_rt)] = actual_scores
    return good_scores

def get_bad_spectral_matches_for(mz, rtdiff, mgf_file1, mgf_file2, tolerance, minmatch):
    from scoring_functions import fast_cosine
    tolerance = tolerance
    minmatch = minmatch
    badscores = []
    badsc = {}

    for spectrum_id1 in mgf_file1:
        spectrum_file1 = mgf_file1[spectrum_id1]
        mz_file1 = spectrum_file1.precursor_mz
        rt_file1 = spectrum_file1.rt



        if mz_file1 >= (mz - 0.01) and mz_file1 <= (mz + 0.01):

            if rt_file1 >= 10:

                for spectrum_id2 in mgf_file2:
                    spectrum_file2 = mgf_file2[spectrum_id2]
                    mz_file2 = spectrum_file2.precursor_mz
                    rt_file2 = spectrum_file2.rt

                    if rt_file2 >= 10:

                        if mz_file2 >= (mz_file1 - 0.01) and mz_file2 <= (mz_file1 + 0.01):

                            if rt_file2 >= (rt_file1 + rtdiff) or rt_file2 <= (rt_file1 - rtdiff):

                                score, used_matches = fast_cosine(spectrum_file1, spectrum_file2, tolerance, minmatch)


                                badscores.append(score)

    return badscores

def get_bad_scores_distribution(peaksets, rtdiff, mgf_file1, mgf_file2):
    bs = {}
    for key in peaksets:
        mz = key[0]
        bad_scores = get_bad_spectral_matches_for(mz, rtdiff, mgf_file1, mgf_file2, 0.2, 2)
        print(peaksets[key])
        print('---')
        print(bad_scores)
        print('=====')
        bs[key] = bad_scores
    return bs

def plot_distrib(bs, peaksets):
    for i in peaksets:
        print(i)
        plt.hist(peaksets[i], bins=50)
        plt.hist(bs[i], bins = 50, alpha = 0.3)

        plt.show()

def calculate_overall(bs, peaksets):
    actual_scores = []
    bad_scores = []
    rts = []
    mz = []
    for i in peaksets:

        actual_score = peaksets[i][0]
        mean_bad_score = np.mean(bs[i])
        actual_scores.append(actual_score)
        bad_scores.append(mean_bad_score)
        rts.append(i[1])
        mz.append(i[0])
    return actual_scores, bad_scores , rts , mz

def calculate_counts_bad_scores(bs, peaksets, score_threshold):
    nonexistent = 0
    higher = 0
    lower = 0
    diffh = []
    diffl = []
    for i in peaksets:

        actual_score = peaksets[i][0]
        mean_bad_score = np.mean(bs[i])

        if score_threshold-0.1 <= actual_score <= score_threshold:

            if mean_bad_score > actual_score:
                diff_higher = mean_bad_score - actual_score
                higher +=1
                diffh.append(diff_higher)

            elif mean_bad_score < actual_score:
                diff_lower = actual_score - mean_bad_score
                lower+=1
                diffl.append(diff_lower)

            elif len(bs[i]) == 0:
                nonexistent+=1
    return higher, diffh, lower, diffl, nonexistent

def get_number_of_pairs_with_no_bad_matches_attached(bs):
    z=0
    for i in bs:
        if i >= 0 :
            z+=1
    return z

def print_threshold(bs, peaksets):
    scts = np.arange(0.1,1.1,0.1)
    nhs = []
    dhs = []
    nls = []
    dls = []
    nes = []
    for i in scts:

        nh,dh,nl,dl,ne = calculate_counts_bad_scores(bs, peaksets, i)
        if len(dh)>0:
            dh = np.mean(dh)
        else:
            dh = 'none'
        if len(dl)>0:
            dl = np.mean(dl)
        else:
            dl = 'none'
        nhs.append(nh)
        dhs.append(dh)
        nls.append(nl)
        dls.append(dl)
        nes.append(ne)
    df = pd.DataFrame([nhs, dhs, nls, dls, nes], ['#higher bad score', 'diff','#lower bad score', 'diff', '#no bad scores' ]).transpose()
    df.index = scts
    return df

def get_ids_for_top_percent(dataframe, per):
    idlist = []

    for rowid,row in dataframe.iterrows():

        total_zero = 0
        total = len(row)
        for i in range(len(row)):

            if float(row[i]) == 0.0:
                total_zero += 1

        percentage = total_zero*100/total
        if percentage <= per:
            idlist.append(rowid)
    return np.array(idlist)

def get_ids_for_zero(dataframe):
    idlist = []

    for rowid,row in dataframe.iterrows():

        total_zero = 0
        total = len(row)
        for i in range(len(row)):

            if float(row[i]) == 0.0:
                total_zero += 1

        percentage = total_zero*100/total
        if percentage == 100:
            idlist.append(rowid)
    return np.array(idlist)

def impute_knn(all_samples, condition_name, k_value = 3):
    from fancyimpute import KNN
    filt = all_samples[condition_name]
    filled_knn = KNN(k=k_value).fit_transform(filt)
    df = pd.DataFrame(filled_knn, columns=filt.columns, index = filt.index)
    return df

def print_boxplots_for_pathway(name, limma_topfeatures):
    for n in pathway_dict[name]:
            for i in n:
                pd.set_option('display.max_rows', 1000)
                print(data.loc[i])
                plot_boxplots(limma_topfeatures,i, (15,30), False)
