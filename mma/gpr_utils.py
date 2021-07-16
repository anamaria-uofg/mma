import numpy as np
import pandas as pd
import GPy
from sklearn import metrics
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import train_test_split


def create_bins_mz_range(dataset, n_bins):
    dataset.sort()
    lowest_value = np.min(dataset)-0.001
    highest_value = np.max(dataset)+0.001
    width = (highest_value - lowest_value)/n_bins
    bins = []
    for i in range(n_bins):
        bins.append((lowest_value+width*i, lowest_value+width*(i+1)))
    nb = [i[0] for i in bins]
    nb.append(bins[-1][1])
    return nb

def return_stratcol(rt, mz, ref_rt, ref_mz, nbins):
    df = pd.DataFrame([rt, mz, ref_rt, ref_mz, np.subtract(ref_rt, rt)]).transpose()
    df.columns = ['rt', 'mz', 'ref_rt', 'ref_mz', 'drift']

    #check correlations between columns
    corr_matrix = df.corr()
    colname = corr_matrix['drift'].sort_values(ascending=False).index[1]

    bins = create_bins_mz_range(list(df[colname]), nbins)

    df['rt_cat'] = pd.cut(df[colname],
                              bins=bins,
                              labels=np.arange(1,nbins+1,1))
    return df['rt_cat']


def extract_outliers(old_mzml_rts, new_mzml_rts, stratcol, zscore_cutoff):
    """
    Eliminating outliers using zscore. A default zscore of 2 is used, which keeps at least 95% of the data points (for normally distributed data).
    """
    old_mzml_rts_no_outliers = []
    new_mzml_rts_no_outliers = []
    stratcol_no_outliers = []


    diff = np.subtract(old_mzml_rts, new_mzml_rts)
    zscores = np.abs(stats.zscore(diff))

    for (old_rt, new_rt, strat, zscore) in zip(old_mzml_rts, new_mzml_rts, stratcol, zscores):
         if (zscore<zscore_cutoff):
            old_mzml_rts_no_outliers.append(old_rt)
            new_mzml_rts_no_outliers.append(new_rt)
            stratcol_no_outliers.append(strat)

    return old_mzml_rts_no_outliers, new_mzml_rts_no_outliers, stratcol_no_outliers


def predict_data(training_data, test_data, training_data_target, test_data_target, k, kname, plot):

    mtest = GPy.models.GPRegression(training_data,training_data_target,k)
    mtest.constrain_positive('')
    mtest.optimize_restarts(10)

    if plot:
        mtest.plot()
        plt.title(kname)
        plt.xlabel('RT (min)')
        plt.ylabel('RT drift (min)')
        plt.show()

    gpr_predicted_data,_ = mtest.predict(test_data)

    reference_rt = test_data + test_data_target
    corrected_rt = test_data + gpr_predicted_data
    corrected_drift = test_data + test_data_target - (test_data + gpr_predicted_data)

    if plot:
        plt.plot(reference_rt, corrected_drift,  '.')
        plt.plot(reference_rt, reference_rt - reference_rt, '.r')
        plt.ylabel('New drift')
        plt.xlabel('Reference data')
        plt.show()

    accuracy = metrics.r2_score(reference_rt, corrected_rt)
    mae = metrics.mean_absolute_error(reference_rt, corrected_rt)
    mse = metrics.mean_squared_error(reference_rt, corrected_rt)

    print('Cross-Predicted Accuracy for',kname,':', accuracy)
    print('Mean absolute error for',kname,':', mae)
    print('Mean squared error for',kname,':', mse)

    return accuracy, mae, mse, mtest


def try_gp_regressions(new_rt, drift, plot, stratcol=None):


    X = np.array(new_rt)[:,None]
    Y = np.array(drift)[:,None] #old_rt -new_rt

    #Define kernels (MLP and RBF)
    k_nn = GPy.kern.MLP(1)
    k_rbf = GPy.kern.RBF(input_dim=1, variance=1, lengthscale=1)

    #Define a vector with kernels which might work
    ks = [k_nn, k_nn+k_rbf,k_nn*k_rbf]
    ks_names = ['MLP','MLP+RBF','MLP*RBF']
    names = ['RBF','MLP','MLP+RBF','MLP*RBF']

    #Cross-validation
    accuracy_list = []
    mae_list = []
    mse_list = []

    if stratcol:
        training_data, test_data, training_data_target, test_data_target = train_test_split(X,
                                                                Y, test_size=0.33, random_state=0, stratify = stratcol)
    else:
        training_data, test_data, training_data_target, test_data_target = train_test_split(X,
                                                                Y, test_size=0.33)




    accuracy0, mae0, mse0, mtest = predict_data(training_data, test_data,
                                                training_data_target, test_data_target, k_rbf, 'RBF', plot)
    accuracy_list.append(accuracy0)
    mae_list.append(mae0)
    mse_list.append(mse0)

    mfinal = mtest
    kfinal = k_rbf
    namefinal = names[0]

    for k, name in zip(ks,ks_names):

        accuracy, mae, mse, mtest = predict_data(training_data, test_data,
                                                 training_data_target, test_data_target, k, name, plot)
        accuracy_list.append(accuracy)
        mae_list.append(mae)
        mse_list.append(mse)

        if mae < mae0:
            if accuracy >= accuracy0:
                accuracy0 = accuracy
                mae0 = mae
                mse0 = mse
                mfinal = mtest
                kfinal = k
                namefinal = name

    results_table = pd.DataFrame(data=[names, accuracy_list, mae_list, mse_list])
    print("Final kernel:",namefinal)
    return kfinal


def get_drift_gpr_model(old_mzml_rts, new_mzml_rts,output_dir, zscore, stratcol = None, plot = False):


    #new_mzml_rts_no_outliers, old_mzml_rts_no_outliers, stratcol_no_outliers = extract_outliers(new_mzml_rts, old_mzml_rts, stratcol, zscore)
    print(stratcol_no_outliers)
    print("Choosing the optimal kernel")
    kernel= try_gp_regressions(np.array(new_mzml_rts_no_outliers),
                                               np.subtract(np.array(old_mzml_rts_no_outliers),
                                                           np.array(new_mzml_rts_no_outliers)), stratcol_no_outliers, plot)


    final_model = GPy.models.GPRegression(np.array(new_mzml_rts_no_outliers)[:,None], np.subtract(np.array(old_mzml_rts_no_outliers),np.array(new_mzml_rts_no_outliers))[:,None],kernel)
    #final_model.Gaussian_noise.variance.fix(0.01)
    final_model.optimize_restarts(10)
    #final_model.optimize(optimizer='scg')
    #final_model.log_likelihood()

    if plot:
        fig, ax = plt.subplots()
        nice_fonts = {"text.usetex": True,
            "font.family": "serif",
            "font.serif" : "Times New Roman"}
        plt.rcParams.update(nice_fonts)

        plt.rc('font', family='serif', size = 15)
        plt.rc('xtick', labelsize=15)
        plt.rc('ytick', labelsize=15)

        final_model.plot()

        #final_model.plot_mean(color = '#000000',ax=ax, label = ('Mean'))
        #final_model.plot_confidence(color = '#AFABAB', ax=ax, label = ('Confidence Interval'))
        #final_model.plot_data(ax=ax, marker = '.')
        #plt.xlim(2,19.5)
        plt.xlabel('Reference RT (min)')
        plt.ylabel('Drift RT (min)')

        plt.show()
        fig.savefig(output_dir+'model.png', format = 'png', dpi =400, bbox_inches='tight')
    print(final_model)

    #display(final_model)
    old_mzml_rts = np.array(old_mzml_rts)[:,None]
    new_mzml_rts = np.array(new_mzml_rts)[:,None]

    gpr_predicted_data,_ = final_model.predict(new_mzml_rts) #predict new RTs for the whole data including outliers

    if plot:
        plotdrift(old_mzml_rts, new_mzml_rts, gpr_predicted_data, output_dir)

    accuracy = metrics.r2_score(old_mzml_rts, new_mzml_rts+gpr_predicted_data)
    mae = metrics.mean_absolute_error(old_mzml_rts, new_mzml_rts+gpr_predicted_data)
    mse = metrics.mean_squared_error(old_mzml_rts, new_mzml_rts+gpr_predicted_data)

    return final_model, accuracy, mae, mse


def plotdrift(old_mzml_rts, new_mzml_rts, gpr_predicted_data, output_dir):
    fig, ax = plt.subplots()
    nice_fonts = {"text.usetex": True,
        "font.family": "serif",
        "font.serif" : "Times New Roman"}
    plt.rcParams.update(nice_fonts)

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=19)
    plt.rc('ytick', labelsize=19)

    ax.scatter(old_mzml_rts, old_mzml_rts - new_mzml_rts - gpr_predicted_data, alpha = 0.7, color = '#32E081', label = 'After GPR correction' )
    ax.scatter(old_mzml_rts, old_mzml_rts - new_mzml_rts , alpha = 0.5, color = '#E02A3F', label = 'Before GPR correction')


    ax.legend(loc = 'best', fontsize = 12)
    plt.ylabel('RT drift (min)')
    plt.xlabel('Reference RT (min)')
    plt.axhline(y=0, color='#000000', alpha = 0.8, linestyle='solid')
    plt.axhline(y=0.5, color='#000000', alpha = 0.6, linestyle='dashed')
    plt.axhline(y=-0.5, color='#000000', alpha = 0.6, linestyle='dashed')

    plt.show()
    fig.savefig(output_dir+'plot.png', format = 'png', dpi =400, bbox_inches='tight')


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
                    old_rt_min = float(tokens[5])
                    old_rt_max = float(tokens[6])

                    tokens[2] = get_predicted_rt(old_rt, model)
                    tokens[4] = get_predicted_rt(old_rt_min, model)
                    tokens[5] = get_predicted_rt(old_rt_max, model)

                    new_line = ','.join([str(t) for t in tokens])
                    g.write(new_line + '\n')
