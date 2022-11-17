# mma
Meta Metabolomics Analysis
An analysis performed for integrating multiple metabolomics datasets
This analysis makes use of different python packages (GPy, pymzmine) and open source software MzMine2. Additionally several functions (found in metab_utils.py and gpr_utils.py) and classes (peakinfo.py) were created to aid in the analysis of the picked LC-MS peaks.

Peakinfo.py: The information about any given peak was stored in a peakinfo object with the following attributes: id, m/z, RT, p-val, t-val, logFC, mummichog annotation, mummichog pathway, mummichog kegg id, std annotation, std kegg id, spectra, adducts, best ms2 match adduct, ms2 annotation, ms2 kegg id, intensities (from each sample).

The analysis performed in Jupyter notebooks contains 4 steps:
  1) Regression and alignment:
    First the retention time (RT) drift between the datasets was studied in StatsRT.ipynb.
    Next, Gaussian Process regression was used to model the RT drift in GPRegressionRT.ipynb. This makes use of the functions in gpr_utils.py and GPy package.
    The optimal RTWindow parameter used in JoinAligner algorithm from MZmine2 was selected following tests with different times for RT for aligning the internal standards in JoinAlignerRTWindow.ipynb.
    In ModifyingRT.ipynb, the RT drift was modified for D_VL and D_M based on the GP models previously obtained.
    In PickAndAlign.ipynb, the three datasets were aligned through JoinAligner by calling MZmine2 though pymzmine package.
    In PeakSetFiltering.ipynb, the aligned peaksets were filtered based on data sparsity.

  2) Fragmentation and data analysis:
    MS2 data from all three datasets was aligned to the previously aligned peaksets in Step 1.

  3) Annotation and pathway analysis:
    Data was annotated using mummichog and HMDB.

  4) Algorithm evaluation:
    The algorithm was evaluated using MS2 data and individual datasets alignment.
