`make.geno` <-
function (data, SNPs.sel) 
{
    if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")

    ans<-togeno(data[,SNPs.sel[1]],sep="/",lab=SNPs.sel[1])
    for (i in 2:length(SNPs.sel))
     {
      ans.i<-togeno(data[,SNPs.sel[i]],sep="/",lab=SNPs.sel[i])
      ans<-cbind(ans,ans.i)
     }
    geno<-setupGeno(ans)
    geno
}

