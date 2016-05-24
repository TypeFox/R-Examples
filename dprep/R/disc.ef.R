disc.ef <-
function (data, varcon, k,out=c("symb","num")) 
{
#Fixed by Esgar Acuna, September  2015
#It handles numerical and symbolic discretization
    p <- dim(data)[2]
    n=dim(data)[1]
    f <- p - 1
   ft <- rep(0, f)
#    for (i in 1:length(varcon)) {
 #       ft[varcon[i]] = 1
  #  }
ft[varcon]=1
    for (i in 1:f) {
a=table(data[,i])
ciclo=ceiling(n/k)
if(length(a[a>ciclo])>0){ft[i]=0}
        if (ft[i] > 0) {
            data[, i] <- disc2(as.vector(data[, i]), k,out)
     }
    }
tempo=length(ft[ft>0])
cat("\n","variables to be discretized:",which(ft==1),"\n") 
data
}
