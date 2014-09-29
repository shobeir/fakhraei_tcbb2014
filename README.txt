These files contain implementation of experimental evaluation in the following paper:

  @article{fakharei2014network,
  title={Network-Based Drug-Target Interaction Prediction with Probabilistic Soft Logic},
  author={Fakhraei, Shobeir and Huang, Bert and Raschid, Louiqa and Getoor, Lise},
  journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  year={2014},
  }

To run the files you need to install Probabilistic Soft Logic (PSL) prerequisites first. Please read the instructions here for details:
https://github.com/linqs/psl/wiki

You can then run the first experiment using 'run.sh'. Please change that file to run the other experiments.

PSL models are located here: /src/main/java/edu/umd/cs/psl/fakhraei_tcbb2014

Perlman's method is implemented in matlab and is located here: /src/main/matlab/Perlmans_method

The full dataset in .mat format is located here: /src/main/matlab/Perlmans_method/data

A copy of the paper can be found here:
http://linqs.cs.umd.edu/basilic/web/Publications/2014/fakhraei:tcbb14/
