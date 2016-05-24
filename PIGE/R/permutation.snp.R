#' The \code{permutation.snp} function performs on permutation sample a Likelihood Ratio Test (or a Wald test) for an interaction 
#' term SNP*E (where E is an Environment variable) or for the effect of the SNP. This function uses the parallel computing on different CPU of the computer. 
#' This function returns a matrix containing for each permutation and SNP the p-value of the interaction
#' term tested (SNP*E) or the p-value of the SNP effect tested.
#' @export
#' @title Parallel computing of the Likelihhod ratio test (or Wald test) for an interaction term (or a simple SNP effect) on permutation sample
#' @name permutation.snp
#' @param model an object of class "formula" : a symbolic description of the model to be fitted
#' without the interaction term.
#' @param Outcome.model a character string naming the type of outcome considered. It could be "\code{binary}" (by default) or "\code{survival}". 
#' @param data a data frame containing the variables in the model. 
#' @param var.inter name of the variable which is tested in interaction with the SNPs (SNP:E). By default var.inter=NULL correspond to a test on the SNPs (no interaction)
#' @param indice.snp vector or character indicating the SNPs to be tested. 
#' @param class.inter class of the \code{var.inter} variable. By default, the variable is considered as continuous and a Wald test is performed. Use ("factor") to indicate categorical variable.
#' @param method method choice for the permutation. By default "YX" a permutation of the phenotype and the adjusted effect are performed otherwise only the phenotype is permuted.
#' @param nbcpu integer indicating the number of CPU of your computer (-1). By default, the function use only
#' one cpu.
#' @param Npermut number of permutation (1000 by default).
#' @param file.out name of the output file where the result will be saved.
#' @return A matrix containing the p-value, for each permutation (row) and for each SNP (column), of the likelihood ratio test (or the Wald test) for the interaction term or the SNP effect.
#' This matrix is also saved in a txt file (named by the argument \code{file.out}) located in the current directory.
#' 
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}
#' 
#' @examples data(data.pige)
#' data(data.pathway)
#' data(list.gene.snp)
#' res <-data.to.PIGE(data=data.pige,data.pathway=data.pathway,
#' list.gene.snp=list.gene.snp,choice.pathway=c(1,2))
#' formul <- formula(y~factor(cov1)+factor(cov2)+factor(cov3)+factor(cov4)
#' +var_int)
#' debut <- Sys.time() 
#' p.snp.permut.ex <-  permutation.snp(model=formul,data=data.pige,
#' indice.snp=res$snp.selected,var.inter="var_int",class.inter=NULL,nbcpu=3,
#' Npermut=9,file.out="res-permut") 
#' print(Sys.time()-debut)

#' ##Survival example:
#' data(data.surv)
#' data(data.pathway.surv)
#' data(list.gene.snp.surv)
#' res1 <-data.to.PIGE(data=data.surv,data.pathway=data.pathway.surv,
#' list.gene.snp=list.gene.snp.surv,choice.pathway=c(1:7))
#' formul <- formula(Surv(TIME, EVENT) ~ var_int)
#' p.snp.permut.ex <-  permutation.snp(model=formul,Outcome.model="surv"
#' ,data=data.surv,indice.snp=res1$snp.selected,var.inter="var_int",
#' class.inter=NULL,nbcpu=3,Npermut=9,file.out="res-permut-surv")

permutation.snp <- function(model,Outcome.model="binary",data,var.inter=NULL,indice.snp,class.inter=NULL,method="YX",nbcpu=NULL,Npermut=1000,file.out="res-permut"){
  mat.snp <- data[,indice.snp]
  if(is.null(nbcpu)){sfInit(parallel=FALSE, cpus=1)}else{
    sfInit(parallel=TRUE, cpus=nbcpu)
    sfLibrary(package="PIGE",pos=4,character.only=TRUE)
  }
  sfExportAll()
#  print(model)
  
  
  if(is.null(var.inter)) {
    if(method=="YX"){
      result <- sfLapply(1:Npermut,permutation.wrapper.cont.Y.and.X,mat=mat.snp,data=data,model=model,Outcome.model=Outcome.model)
      }else{
        result <- sfLapply(1:Npermut,permutation.wrapper.cont,mat=mat.snp,data=data,model=model,Outcome.model=Outcome.model)
      }
  }else{
  if(is.null(class.inter)){
    var.inter <- data[,var.inter]
    if(method=="YX"){
    result <- sfLapply(1:Npermut,permutation.wrapper.cont.inter.Y.and.X,mat=mat.snp,data=data,model=model,var.inter=var.inter,Outcome.model=Outcome.model)
    }else{
    result <- sfLapply(1:Npermut,permutation.wrapper.cont.inter,mat=mat.snp,data=data,model=model,var.inter=var.inter,Outcome.model=Outcome.model)
    }
    }else{
    var.inter <- factor(data[,var.inter])
    if(method=="YX"){
    result <- sfLapply(1:Npermut,permutation.wrapper.cat.inter.Y.and.X,mat=mat.snp,data=data,model=model,var.inter=var.inter,Outcome.model=Outcome.model)
    }else{
    result <- sfLapply(1:Npermut,permutation.wrapper.cat.inter,mat=mat.snp,data=data,model=model,var.inter=var.inter,Outcome.model=Outcome.model)  
    }
    }}  
    sfStop()
  mat.result <- matrix(unlist(result),ncol=length(indice.snp),nrow=Npermut,byrow=TRUE)
  colnames(mat.result) <- colnames(mat.snp)
  write.table(mat.result,file=paste(file.out,".txt",sep=""),row.names=FALSE,col.names=TRUE)
  return(mat.result)
}


