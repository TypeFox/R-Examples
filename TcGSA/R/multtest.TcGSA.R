#'Computing the P-value of the Likelihood Ratios Applying a Multiple Testing
#'Correction
#'
#'This function computes the p-value of the likelihood ratios and apply a
#'multiple testing correction.
#'
#'
#'@param tcgsa 
#'a TcGSA object.
#'
#'@param threshold 
#'
#'the threshold at which the FDR or the FWER should be
#'controlled.
#'
#'@param myproc 
#'a vector of character strings containing the names of the
#'multiple testing procedures for which adjusted p-values are to be computed.
#'This vector should include any of the following: "\code{Bonferroni}",
#'"\code{Holm}", "\code{Hochberg}", "\code{SidakSS}", "\code{SidakSD}",
#'"\code{BH}", "\code{BY}", "\code{ABH}", "\code{TSBH}" or "\code{none}".  
#'"\code{none}" indicates no adjustement for multiple testing.  See
#'\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for details.  Default is
#'"\code{BY}", the Benjamini & Yekutieli (2001) step-up FDR-controlling
#'procedure (general dependency structures).  In order to control the FWER (in
#'case of an analysis that is more a hypothesis confirmation than an
#'exploration of the expression data), we recommand to use "\code{Holm}", the
#'Holm (1979) step-down adjusted p-values for strong control of the FWER.
#'
#'@param nbsimu_pval 
#'the number of observations under the null distribution to
#'be generated in order to compute the p-values. Default is \code{1e+06}
#'
#'@return \code{multtest.TcGSA} returns an dataframe with 5 variables.  The
#'rows correspond to the gene sets under scrutiny.  The 1st column is the
#'likelihood ratios \code{LR}, the 2nd column is the convergence status of the
#'model under the null hypothesis \code{CVG_H0}, the 3rd column is the
#'convergence status of the model under the alternative hypothesis
#'\code{CVG_H1}, the 4th column is the raw p-value of the mixed likelihood
#'ratio test \code{raw_pval}, the 5th column is the adjusted p-value of the
#'mixed likelihood ratio test \code{adj_pval}.
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{TcGSA.LR}}, 
#'\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}},
#'\code{\link{signifLRT.TcGSA}}
#'
#'@importFrom multtest mt.rawp2adjp
#'
#'@importFrom stats rchisq
#'
#'@export multtest.TcGSA
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'                           
#'mtt <- multtest.TcGSA(tcgsa_sim_1grp, threshold = 0.05, 
#'                      myproc = "BY", nbsimu_pval = 1000)
#'mtt
#'
#'
multtest.TcGSA <-
function(tcgsa, threshold=0.05, myproc="BY", nbsimu_pval = 1000000){
  emp <- tcgsa[["fit"]]
  func <- tcgsa[["time_func"]]
  group.var <- tcgsa[["group.var"]]
  separateSubjects <- tcgsa[["separateSubjects"]]
	time_DF <- tcgsa[["time_DF"]]
  
  if(is.null(group.var)){
    if(!separateSubjects){
      if(func=="linear"){
        theodist <- c(stats::rchisq(nbsimu_pval/2,df=1), stats::rchisq(nbsimu_pval/2,df=2))
      }else if(func=="cubic"){
        theodist <- rmixchisq(nbsimu_pval,3,3)
      }else{
        theodist <-rmixchisq(nbsimu_pval,time_DF,time_DF)
      }
    }else{
      if(func=="linear"){
        theodist <- c(stats::rchisq(nbsimu_pval/2,df=0), stats::rchisq(nbsimu_pval/2,df=1) )
      }else if(func=="cubic"){
        theodist <- rmixchisq(nbsimu_pval,0,3)
      }else{
        theodist <-rmixchisq(nbsimu_pval,0,time_DF)
      }
    }
  }else{
    nbgp <- length(levels(group.var))
    if(func=="linear"){
      theodist <- rmixchisq(nbsimu_pval, 1*(nbgp-1), 0)
    }else if(func=="cubic"){
      theodist <- rmixchisq(nbsimu_pval, 3*(nbgp-1), 0)
    }else{
      theodist <-rmixchisq(nbsimu_pval, time_DF*(nbgp-1), 0)
    }
  }
  
  emp$raw_pval <- unlist(lapply(emp$LR, FUN=pval_simu, theo_dist=theodist))
  if(myproc=="none" | length(emp$raw_pval)==1){
  	emp$adj_pval <- emp$raw_pval
  }
  else{
  	adj_pval <- mt.rawp2adjp(emp$raw_pval, proc=c(myproc),alpha=threshold)
  	emp$adj_pval <- adj_pval$adjp[order(adj_pval$index),2]
  }
  return(emp)
}
