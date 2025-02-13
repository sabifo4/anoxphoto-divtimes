          seed = -1
       seqfile = /mnt/c/Users/sandr/Sync/00_Collabs/Ed_Rich/anoxphoto-divtimes/01_PAML/01_MCMCtree/dummy_aln/dummy_aln.aln
      treefile = /mnt/c/Users/sandr/Sync/00_Collabs/Ed_Rich/anoxphoto-divtimes/00_data_formatting/01_inp_data/anoxphoto_withArchExclDPANN_calib_MCMCtree.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1    * Number of partitions
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 0    * 0: no data (prior); 1:exact likelihood;
                      * 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates

         model = 3        * models for AAs or codon-translated AAs:
                          *     0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                          *     6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5      * alpha for gamma rates at sites
         ncatG = 4        * No. categories in discrete gamma
    aaRatefile = lg.dat   * Path to the file with the LG matrix

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 2.7   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10    * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = -1       * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
