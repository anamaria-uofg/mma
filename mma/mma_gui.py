import streamlit as st
import numpy as np
import pandas as pd
import re
from io import StringIO
import altair as alt
import itertools
import matplotlib.pyplot as plt
import GPy
from gpr_utils import try_gp_regressions
from gpr_utils import extract_outliers



def main():
    #setting the title
    backgroundColor = '#7792E3'
    st.title('GPR RT Drift Correction')
    st.write('This application calculates the retention time drift between two datasets, based on n\
        their standard reference mixture files. Each metabolite\'s m/z and RT are extracted and compared n\
        . The first dataset uploaded is set as the reference dataset. And the RT drift in the other datasets n\
        is plotted. A Gaussian Process regression model is fitted to the data and the final model can be saved.')



    #Sidebar stuff

    st.sidebar.title('Upload SRM files here:')

    datasets_range = [2,3,4,5]

    datasets_number_option = st.sidebar.selectbox(
        'How many datasets would you like to integrate?', datasets_range)

    st.sidebar.write('You selected:', datasets_number_option, 'datasets.')



    datasets = [None]*datasets_number_option
    srms = [None]*datasets_number_option

    for i in range(datasets_number_option):
        datasets[i] = st.sidebar.file_uploader('Upload files containing SRM information for dataset '+str(i+1),
                                                type = ['csv'], key = chr(i+1),
                                                accept_multiple_files = True )

        srms[i] = []




    left_column, right_column = st.sidebar.beta_columns(2)
    pressed = left_column.checkbox('Calculate RT drift')
    if pressed:
        right_column.write("Calculation in progress...")
        for i in range(datasets_number_option):
            if datasets[i] is not None:
                for srm_file in datasets[i]:
                    for line in srm_file:
                        #to skip the comments read only the rows which start with a number
                        string_line = line.decode('utf-8')
                        if re.match(r"^\d+.*$",string_line):
                            name = line.decode().split(',')[2].lower()
                            polarity = line.decode().split(',')[4]
                            mz = float(line.decode().split(',')[6])
                            rt = float(line.decode().split(',')[9])
                            # mzs.append(mz)
                            # rts.append(rt)
                            # names.append(name)

                            srms[i].append((name, mz, rt))
                    srms[i].sort(key = lambda x:int(x[1]))

        res = [(s1[0], s1[2], s1[2] - s2[2]) for s1 in srms[0]
        for s2 in srms[1] if s1[0] == s2[0]]
        #p = list(itertools.chain(*[e for e in zip(srms[len(srms)-1],srms[len(srms)]) if e[0][0] == e[1][0]]))
        #st.write(res)

        df = pd.DataFrame(res, columns = ['Name', 'Reference RT (min)', 'RT drift (min)'])
        c = alt.Chart(df).mark_circle(size =100).encode(x='Reference RT (min)',y=alt.Y('RT drift (min)'),
                                                opacity=alt.value(0.2),
                                                color=alt.value('red'),

                                                tooltip=['Reference RT (min)', 'RT drift (min)', 'Name'])
                    #scale=alt.Scale(domain=[-10, 10])),

        st.altair_chart(c, use_container_width = True)
        right_column.write('')

        drift_mean = np.mean(df['RT drift (min)'])
        if drift_mean <30:
            st.write('The drift is quite minimal. Do you wish to continue with the anaysis?')
        else:
            st.write('Do you wish to continue with the anaysis?')

        left_column, right_column = st.beta_columns(2)
        yes_pressed = left_column.checkbox('Yes. Continue analysis.')
        no_pressed = right_column.checkbox('No.')
        # option = st.selectbox(
        #     'How would you like to be contacted?',
        #         ('Email', 'Home phone', 'Mobile phone'))
        #



        if yes_pressed:

            st.write('Testing GP kernels...')

            # # Add a placeholder
            # latest_iteration = st.empty()
            # bar = st.progress(0)
            # #
            # for i in range(100):
            #
            #     latest_iteration.text(f'Iteration {i+1}')
            #     bar.progress(i + 1)
            #     time.sleep(0.1)
            #     '...and now we\'re done!'

            final_kernel = try_gp_regressions(df['Reference RT (min)']-df['RT drift (min)'], df['RT drift (min)'],
                                                plot=True)




            final_model = GPy.models.GPRegression(np.array(df['Reference RT (min)']-df['RT drift (min)'])[:,None],
            np.array(df['RT drift (min)'])[:,None] ,final_kernel)
            #final_model.Gaussian_noise.variance.fix(0.01)
            final_model.optimize_restarts(10)
            mdf = pd.DataFrame([final_model.param_array], columns = final_model.parameter_names(), index=None)
            st.dataframe(mdf.transpose())

            fig, ax = plt.subplots()
            # nice_fonts = {"text.usetex": True,
            #     "font.family": "serif",
            #     "font.serif" : "Times New Roman"}
            # plt.rcParams.update(nice_fonts)

            plt.rc('font',  size = 15)
            plt.rc('xtick', labelsize=15)
            plt.rc('ytick', labelsize=15)

            #final_model.plot()

            final_model.plot_mean(color = '#000000',ax=ax, label = ('Mean'))
            final_model.plot_confidence(color = '#AFABAB', ax=ax, label = ('Confidence Interval'))
            final_model.plot_data(ax=ax, marker = '.')
            # plt.xlim(2,19.5)
            plt.xlabel('Reference RT (min)')
            plt.ylabel('Drift RT (min)')

            st.pyplot(fig)





    # st.line_chart(chart_data)
    #
    # st.write('Plotting a map')
    # map_data = pd.DataFrame(
    #     np.random.randn(1000, 2) / [50, 50] + [37.76, -122.4],
    #     columns=['lat', 'lon'])
    #
    # st.map(map_data)
    #
    # if st.checkbox('Show dataframe'):
    #     chart_data = pd.DataFrame(
    #        np.random.randn(20, 3),
    #        columns=['a', 'b', 'c'])
    #
    #     st.line_chart(chart_data)
    #
    # option = st.selectbox(
    #     'Which number do you like best?',
    #      df['first column'])
    #
    # 'You selected: ', option
    #
    #
    #
    # left_column, right_column = st.beta_columns(2)
    # pressed = left_column.button('Press me?')
    # if pressed:
    #     right_column.write("Woohoo!")
    #
    # expander = st.beta_expander("FAQ")
    # expander.write("Here you could put in some really, really long explanations...")
    #
    #
    #
    #
    #st.write('Analysis stopped.')

if __name__ == '__main__':
    main()
