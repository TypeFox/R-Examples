disc.ew <-
function (data, varcon, out=c("symb","num")) 
{
#Fixed by Edgar Acuna, September 2015
#It handles discretization into integer numbers and intervals 
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    p <- dim(data)[2]
    f <- p - 1
    ft <- rep(0, f)
#    for (i in 1:length(varcon)) {
#        ft[varcon[i]] <- 1
#    }
ft[varcon]=1
if(out=="symb"){
    for (i in 1:f) {
        if (ft[i] > 0) {
            grupos <- nclass.scott(data[, i])
            a=min(data[,i])
            b=max(data[,i])
            w=(b-a)/grupos
            cutpoints=seq(a,b,w)
            cutpoints[1]=-Inf
            cutpoints[grupos+1]=Inf
            data[, i] <- as.vector(cut(data[, i], cutpoints))
        }
    }
}
else{
 for (i in 1:f) {
        if (ft[i] > 0) {
            grupos <- nclass.scott(data[, i])
            data[, i] <- as.vector(cut(data[, i], grupos,labels=FALSE))
        }
    }
}
    data
}
