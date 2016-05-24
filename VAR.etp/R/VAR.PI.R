VAR.PI <-
function(fstat,alpha,f){
dims <- dim(fstat)
h <- dims[1]; k <- dims[2]
intv <- matrix(NA,h,2*k)
for(i in 1:h){
    index <- 1:2
    for(j in 1:k){
    tem <- quantile(fstat[i,j,],c(0.5*(1-alpha),(1-0.5*(1-alpha))))
    intv[i,index] <- tem
    index <- index+2}
    }

rownames(intv) <- rownames(f)
tem1 <- paste(colnames(f),"_lower",sep="")
tem2 <- paste(colnames(f),"_upper",sep="")
tem <- character(2*k)
tem[seq(1,(2*k-1),2)] <- tem1;tem[seq(2,(2*k),2)] <- tem2;
colnames(intv) <- tem
return(list(PI=intv,Prob=alpha)) }
