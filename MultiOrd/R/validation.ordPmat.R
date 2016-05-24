validation.ordPmat <-
function(ordPmat){
ncatvec = numeric(0)
w1=ordPmat>1
w2=ordPmat<0
w3=ordPmat==1
w4=is.na(ordPmat)
if(sum(w4)>0){stop("NAs are not permitted in ordPmat")}
if(sum(w1)>0){stop("The probabilities cannot be greater than 1")}
if(sum(w2)>0){stop("The probabilities cannot be less than 0")}
if(sum(w3)>0){stop("The probabilities cannot be equal to 1")}

noutcome = ncol(ordPmat)
for(i in 1:noutcome){
if (sum(ordPmat[,i])<1 | sum(ordPmat[,i])>1) {
stop("\n Invalid probabilities for variable", i,".", "\n Sum of probabilities for each variable must be exactly 1.")
}
ncatvec[i]= length(ordPmat[ordPmat[,i]>0,i])
}
return(list(J=noutcome, K=ncatvec) )
}
