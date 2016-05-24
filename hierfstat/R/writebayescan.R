#' @title Writes a bayescan file
#' 
#' @description write the genotypes in a format suitable for
#' analysis with bayescan
#' 
#' @usage write.bayescan(dat=dat,diploid=TRUE,fn="dat.bsc")
#' 
#' @param dat a genotype data frame
#' @param diploid whether the dataset is diploid or haploid
#' @param fn  file name for output
#' 
#' @return a text file fn is written in the current directory
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @references \href{http://www.genetics.org/content/180/2/977.abstract}{Foll M and OE Gaggiotti (2008) 
#'  Genetics 180: 977-993}
#'
#' 
#' \url{http://cmpg.unibe.ch/software/BayeScan/}
#'
#'@export
#'
#'
##############################################
write.bayescan<-function(dat=dat,diploid=TRUE,fn="dat.bsc"){
  nloc<-dim(dat)[2]-1
  npop<-length(table(dat[,1]))
  alc.dat<-allele.count(dat,diploid)
  nal<-unlist(lapply(alc.dat,function(x) dim(x)[1]))
  nindx<-sapply(alc.dat,function(x) apply(x,2,sum))
  write(paste("[loci]=",nloc,sep=""),fn)
  write("",fn,append=TRUE)
  write(paste("[populations]=",npop,sep=""),fn,append=TRUE)
  write("",fn,append=TRUE)
  for (ip in 1:npop){
    write("",fn,append=TRUE)
    write(paste("[pop]=",ip,sep=""),fn,append=TRUE)
    for (il in 1:nloc){
      tow<-c(il,nindx[ip,il],nal[il],alc.dat[[il]][,ip])
      write(tow,fn,
            append=TRUE,ncolumns=length(tow))    
    }
  }
}
