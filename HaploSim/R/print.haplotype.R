Class <- function(x)
  cat(paste('An object of class',class(x),"\n",sep =" "))

print.haplotype <- function (x,quote = FALSE,digits = 2,
                             na.print = "",zero.print = "0",...)
{
  Class(x)
  cat(paste("hID:", x@hID,"phID0:",x@phID0,"phID1:",x@phID1,"\n",sep = " "))
  if(length(x@snp>0)){
    cat("SNP:\n")
    xx <- x@snp
    class(xx) <- NULL
    print(xx,quote = quote,...)
  }
  if(length(x@qtl)>0){
    nTrait <- max(sapply(x@qtl,function(y)length(y)))
    xx <- matrix(round(unlist(x@qtl),digits = digits),nrow = nTrait)
    colnames(xx) <- names(x@qtl)
    for(r in 1:nrow(xx)){
      cat("QTL for trait",r,"\n",sep = " ")
      print(xx[r,])
    }
  }
  invisible(x)
}

setMethod("print","haplotype",print.haplotype)
setMethod("show","haplotype", function(object)print.haplotype(object))

print.haploList <- function(x)
  {
    Class(x)
    cat(paste("genDist =",x@genDist,"Morgan\n"))
    cat(paste("nDec =",x@nDec,"\n"))
    cat(paste("nChrom =",x@nChrom,"\n"))
    cat(paste('Contains',length(x),'objects of class "haplotype"\n',sep = " "))
    hh <- getAll(x)
    if(!is.null(dim(hh)))
      cat(paste('There are',ncol(getAll(x)),'polymorphic marker loci \n',sep = " "))
    qq <- getAll(x,what = "q")
    if(!is.null(dim(qq)))
      cat(paste('There are',ncol(qq),'polymorphic QTL',
                ifelse(dim(qq)[3] > 1,paste("for",dim(qq)[3],"traits \n"),"\n")))
    invisible(x)
  }
setMethod("print","haploList",print.haploList)
setMethod("show", "haploList", function(object)print.haploList(object))
