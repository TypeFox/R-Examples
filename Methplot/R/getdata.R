#' This function read the output from Methpup into R
#'
#' @param filelist the directory where the output files are saved in
#' @param gene the region that you are interested to look at the methylation profile
#' @param n the number of CpG sites in the region that you specified in "gene".

#' @return This function could yield a dataframe saying the read number detected in each methylation pattern in the given region in all samples under "filelist" directory.
#' @export
#' @author Xin Yang \email{xin.yang@@cimr.cam.ac.uk}
#' @examples
#' foxp3<-getdata(system.file("extdata", package="Methplot"), 10, "foxp3")

getdata<-function(filelist,n, gene){
  L<-matrix(rep(c("C", "T"),times=n), ncol=2, byrow=T)
  L<-split(L, rep(1:nrow(L)))
  pattern<-do.call(paste, c(do.call(expand.grid, L), sep=""))
  
  files<-list.files(filelist, full.names=TRUE)
  a<-basename(files)
  a<-do.call("rbind", strsplit(a, split=".", fixed=TRUE))[,1]

  data<-matrix(0,nrow=length(pattern), ncol=length(a), dimnames=list(pattern,a))
  for(i in 1:length(files)) {
    if(file.info(files[i])$size>0){
      x <- read.csv(files[i],sep="\t",quote="",col.names=c("region","pattern"),as.is=TRUE)
      drop <- c(grep("A|G|N|D",x$pattern),which(as.numeric(nchar(x$pattern))!=n))
      if(length(drop)>0) x <- x[-drop,]
      x<-x[grep(gene, x$region),]
      if (nrow(x)>0) {
        x <- data.frame(n=tapply(1:nrow(x),x$pattern,length),
                pattern=tapply(as.character(x$pattern),x$pattern,unique))
        data[rownames(x),i]<-x$n
      }
    }
  }
  return(data)
}

