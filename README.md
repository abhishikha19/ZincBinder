# ZincBinder
Zinc Bnder is SVM based prediction method for predicting zinc metal binding sites using protein sequence information.
Zinc is one the most abundant catalytic cofactor and also an important structural component of a large number of metallo-proteins. Hence prediction of zinc metal binding sites in proteins can be a significant step in annotation of molecular function of a large number of proteins. Majority of existing methods for zinc-binding site predictions are based on a dataset of proteins, which has been compiled nearly a decade ago. Hence there is a need to develop zinc-binding site prediction system using the current updated data to include recently added proteins. Herein, we propose a SVM based method, named as ZincBinder, for prediction of zinc metal-binding site in a protein using sequence profile information. The predictor was trained using five-fold cross validation approach and achieved 85.37% sensitivity with 86.20% specificity during training. Benchmarking on an independent non-redundant dataset, which was not used during training, showed better performance of ZincBinder vis-à-vis existing methods. Executable versions, source code, sample datasets, and usage instructions are available at http://proteininformatics.org/mkumar/znbinder/
