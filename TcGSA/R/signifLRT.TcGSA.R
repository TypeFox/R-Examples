#'Identifiying the Significant Gene Sets
#'
#'A function that identifies the significant gene sets in an object of class
#''\code{TcGSA}'.
#'
#'
#'@param tcgsa 
#'a \code{tcgsa} object.
#'
#'@param threshold 
#'the threshold at which the FDR or the FWER should be
#'controlled.
#'
#'@param myproc 
#'a vector of character strings containing the names of the
#'multiple testing procedures for which adjusted p-values are to be computed.
#'This vector should include any of the following: "\code{Bonferroni}",
#'"\code{Holm}", "\code{Hochberg}", "\code{SidakSS}", "\code{SidakSD}",
#'"\code{BH}", "\code{BY}", "\code{ABH}", "\code{TSBH}" or "\code{none}".  
#'"\code{none}" indicates no adjustement for multiple testing. See
#'\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for details.  Default is
#'"\code{BY}", the Benjamini & Yekutieli (2001) step-up FDR-controlling
#'procedure (general dependency structures).  In order to control the FWER(in
#'case of an analysis that is more a hypothesis confirmation than an
#'exploration of the expression data), we recommand to use "\code{Holm}", the
#'Holm (1979) step-down adjusted p-values for strong control of the FWER.
#'
#'@param nbsimu_pval 
#'the number of observations under the null distribution to
#'be generated in order to compute the p-values.  Default is \code{1e+06}.
#'
#'@param write 
#'logical flag enabling the export of the results as a table in a
#'.txt file.  Default is \code{FALSE}.
#'
#'@param txtfilename 
#'a character string with the name of the .txt file in which
#'the results table is to be written, if \code{write} is \code{TRUE}. Default
#'is \code{NULL}.
#'
#'@param directory 
#'if \code{write} is \code{TRUE}, a character string with the
#'directory of the .txt file in which the results table is to be written, if
#'\code{write} is \code{TRUE}. Default is \code{NULL}.
#'
#'@return \code{signifLRT.TcGSA} returns a list. 
#'
#'The fisrt element \code{mixedLRTadjRes} is data frame with \eqn{p} rows (one
#'row for each significant gene set) and the 3 following variables:
#'\itemize{
#'\item GeneSet the significant gene set name from the gmt object.
#'\item AdjPval the adjusted p-value corresponding to the signicant gene
#'set.
#'\item desc the significant gene set description from the gmt object.
#'}
#'
#'The second element \code{multCorProc} passes along the multiple testing 
#'procedure used (from the argument \code{myproc}).
#'
#'The third element \code{threshold} passes along the significance threshold
#' used (from the argument \code{threshold}).
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{multtest.TcGSA}}, \code{\link{TcGSA.LR}}
#'
#'@references Hejblum BP, Skinner J, Thiebaut R, (2015) 
#'Time-Course Gene Set Analysis for Longitudinal Gene Expression Data. 
#'\emph{PLoS Computat Biol} 11(6): e1004310.
#'doi: 10.1371/journal.pcbi.1004310
#'
#'@importFrom gtools mixedorder
#'
#'@importFrom utils write.table
#'
#'@export signifLRT.TcGSA
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'                           
#'sgnifs <- signifLRT.TcGSA(tcgsa_sim_1grp, threshold = 0.05, myproc = "BY",
#'                          nbsimu_pval = 1000, write=FALSE)
#'sgnifs
#'
#'
signifLRT.TcGSA <-
function(tcgsa, threshold=0.05, myproc="BY", nbsimu_pval = 1e+06, write=F, txtfilename=NULL, directory=NULL){  
  gmt  <-  tcgsa[["GeneSets_gmt"]]
  signif <- multtest.TcGSA(tcgsa, threshold, myproc, nbsimu_pval)

  signif_mod <- gmt$geneset.name[which(signif$adj_pval<threshold)]
  if(!is.null(signif_mod)){
    signif_desc <- gmt$geneset.descriptions[which(signif$adj_pval<threshold)]
    if(is.null(signif_desc)){
    	signif_desc <- rep("", length(which(signif$adj_pval<threshold)))
    }
    AdjPval <- signif$adj_pval[which(signif$adj_pval<threshold)]
  }else{
    signif_desc <- NULL
    AdjPval <- NULL
  }

  Res_Linear_Mod_FDR <- cbind.data.frame("GeneSet"=signif_mod, "AdjPval"=AdjPval, "desc"=signif_desc)
  if(dim(Res_Linear_Mod_FDR)[1]>1){
    Res_Linear_Mod_FDR <- Res_Linear_Mod_FDR[mixedorder(gsub('.',"a",as.character(Res_Linear_Mod_FDR$GeneSet), fixed=TRUE)),]
  }
  
 if(write){
   if(!is.null(txtfilename)){
     if (is.null(directory)){
        directory <- getwd()
        cat("Warning: 'directory' argument is empty, output file written in the current working directory")
      }
      utils::write.table(Res_Linear_Mod_FDR, file=paste(directory, txtfilename, sep="/"), row.names=FALSE, sep="\t")
    }else{
      cat("ERROR: could not write the significant results file because the argument 'txtfilename' is empty")
    }
  }
  return(list("mixedLRTadjRes"=Res_Linear_Mod_FDR, "multCorProc"= myproc, "threshold"=threshold))
}
