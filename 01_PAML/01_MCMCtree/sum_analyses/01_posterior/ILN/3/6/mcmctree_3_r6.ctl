          seed = -1
       seqfile = /data/scratch/btx126/spiders_dating/alignments/3/aln_supaln.phy
      treefile = /data/scratch/btx126/spiders_dating/trees/calibrated/3/supaln_calib_MCMCtree.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2 ./in.BV
                      * 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 2

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 7.2   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10    * gammaDir prior for sigma^2     (for clock = 2

         print = 1       * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
