#' The \code{data.to.list.gene.snp} function generates the list containing for each gene the corresponding SNP ids from a txt file containing the SNP/gene annotation data. This list is required for the functions \code{data.to.PIGE} and \code{ARTP.GE}.
#' @export
#' @title Annotation gene/SNP
#' @name data.to.list.gene.snp
#' @param file txt file containing the SNP/gene annotation (see Details for the format of the txt file).
#' @param header a logical value indicating whether the file contains the names of the variables as its first line. 
#' By default \code{header} is set to TRUE.
#' @param path linking to the directory containing the data (SNP/gene)
#' @details The txt file containing the annotation data for the SNP/gene is a two columns matrix. The first colums is the name of each SNP.
#' The second columns indicates the genes which the corresponding SNP (same row) belongs. It will be noted "Gene1/Gene4/Gene5", if for example a SNP belongs to the genes: Gene1, Gene4, Gene5.
#' @return A list containing for each gene the names of the SNPs belonging to it. This list is required for using the function
#' \link{data.to.PIGE} and the function \link{ARTP.GE}.
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}
#' @examples 
#' ##Example : case-control study data
#' data(data.pige)
#' data(data.pathway)
#' path.in <- paste(system.file("sampleData", package="PIGE"),"/",sep="")
#' file <- "SNP-GENE-annotation.txt"
#' list.gene.snp <- data.to.list.gene.snp(file,path=path.in) 

#' ##Example Survival data
#' data("data.surv")
#' data("data.pathway.surv")
#' path.in <- paste(system.file("sampleData", package="PIGE"),"/",sep="")
#' file="snp.gene.surv.txt"
#' list.gene.snp.surv <- data.to.list.gene.snp(file,path=path.in)      



data.to.list.gene.snp <- function(file,header=TRUE,path=NULL) {
  name.data <- paste(path.expand(path),file,sep="")
  data.snp <- read.table(file=name.data,header=header)
  snp.gene <- function(x) unlist(strsplit(x, split="/"))
  list.snp.gene <- apply(matrix(as.character(data.snp[,2]),nrow=1,ncol=dim(data.snp)[1]),MARGIN=2,FUN=snp.gene)
  names(list.snp.gene) <- as.character(data.snp[,1]) 
  
  nb.snp.gene <- sapply(list.snp.gene,FUN=length)
  names(list.snp.gene) -> name.snp
  name.gene.snp <- rep(name.snp,times=nb.snp.gene)
  gene.snp <- unlist(list.snp.gene)
  names(gene.snp) <- name.gene.snp
  
  list.gene.snp <- NULL
  for (gene in unique(gene.snp)) {
    list.gene.snp <- c(list.gene.snp,list(gene=names(gene.snp)[which(gene.snp%in%gene)]))  
  }
  names(list.gene.snp) <- unique(gene.snp)
  
  return(list.gene.snp)
  
}
