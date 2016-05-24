#' The \code{compute.p.snp.obs} function performs a Likelihood Ratio Test (LRT) for an interaction term SNP*E (where E is an environmental variable) or for the effect of the SNP.
#' The function return the p-value of the Likelihood Ratio Test performed for all the SNPs in the original data set.
#' The SNPs are considered as continuous variable (coded 0,1,2).
#' @export
#' @title Likelihhod ratio test (or Wald test) for an interaction term (or a single SNP effect)
#' @name compute.p.snp.obs
#' @param data a data frame containing the variables in the model
#' @param model an object of class "formula": a symbolic description of the model to be fitted
#' without the interaction term
#' @param Outcome.model a character string naming the type of outcome considered. It could be "\code{binary}" (by default) or "\code{survival}". 
#' @param indice.snp vector or character indicating the SNPs to be tested 
#' @param var.inter name of the variable which is tested in interaction with the SNPs (SNP:E). By default var.inter=NULL correspond to a test on the SNPs (no interaction)
#' @param class.inter class of the var.inter variable. By default, the variable is considered continuous and a wald test is performed.
#' Use "\code{factor}" to indicate a categorical variable. 
#' @param file.out name of the output file where the result will be saved.
#' @return A data frame with one row containing the pvalue of the likelihood ratio test (or Wald test) for the interaction term (or a SNP effect) for all the SNPs.
#' This data frame is also saved in a txt file (named by the argument \code{file.out}) located in the current directory.
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}
#' @examples data(data.pige)
#' ## Case-control study:
#' data(data.pathway)
#' data(list.gene.snp)
#' res <-data.to.PIGE(data=data.pige,data.pathway=data.pathway,
#' list.gene.snp=list.gene.snp,choice.pathway=c(1,2))
#' formul <- formula(y~factor(cov1)+factor(cov2)+factor(cov3)+factor(cov4)
#' +var_int)
#' p.snp.obs.ex <-  compute.p.snp.obs(data=data.pige,model=formul,
#' indice.snp=res$snp.selected,var.inter="var_int",class.inter=NULL) 
#'
#' ## Survival analysis
#' data(data.surv)
#' data(data.pathway.surv)
#' data(list.gene.snp.surv)
#' res1 <-data.to.PIGE(data=data.surv,data.pathway=data.pathway.surv,
#' list.gene.snp=list.gene.snp.surv,choice.pathway=c(1:7))
#' formul <- formula(Surv(TIME, EVENT) ~ var_int)
#' p.snp.obs.ex <- compute.p.snp.obs(data=data.surv,Outcome.model="surv"
#' ,model=formul,indice.snp=res1$snp.selected,var.inter="var_int"
#' ,class.inter=NULL,file.out="res-obs-surv")


compute.p.snp.obs <- function(data,Outcome.model="binary",model,indice.snp,var.inter,class.inter="NULL",file.out="res-obs"){
  mat <- data[,indice.snp]
  
  
  if(Outcome.model=="binary"){
  if(is.null(var.inter)) {
    result <- apply(mat,MARGIN=2,FUN=LR.cont,formula=model,data=data)#,class.Z=class.inter)    
  }else{
  if(is.null(class.inter)){
    var.inter <- data[,var.inter]
    result <- apply(mat,MARGIN=2,FUN=LR.inter.cont,formula=model,data=data,Z1=var.inter)#,class.Z=class.inter)
  }else{
   var.inter <- factor(data[,var.inter])
  result <- apply(mat,MARGIN=2,FUN=LR.inter.cat,formula=model,data=data,Z1=var.inter)#,class.Z=class.inter)
  } 
  }
  }else{
    if(is.null(var.inter)) {
      result <- apply(mat,MARGIN=2,FUN=LR.cont.surv,formula=model,data=data)#,class.Z=class.inter)    
    }else{
      if(is.null(class.inter)){
        var.inter <- data[,var.inter]
        result <- apply(mat,MARGIN=2,FUN=LR.inter.cont.surv,formula=model,data=data,Z1=var.inter)#,class.Z=class.inter)
      }else{
        var.inter <- factor(data[,var.inter])
        result <- apply(mat,MARGIN=2,FUN=LR.inter.cat.surv,formula=model,data=data,Z1=var.inter)#,class.Z=class.inter)
      } 
    } 
  }

  write.table(as.matrix(t(result)),file=paste(file.out,".txt",sep=""),row.names=FALSE,col.names=TRUE)
  result <- as.data.frame(t(result))
  return(result)
}


