ardec.lm <-
function(x) {

require(stats)

 dat=x-mean(x)
 ndat=length(dat)
 
 p=ar(dat,method="burg")[[1]] # p=order of autoregressive model from AIC (burg method)

 # linear autoregressive model fit
 
 X=t(matrix(dat[rev(rep((1:p),ndat-p)+ rep((0:(ndat-p-1)),rep(p,ndat-p)))],p,ndat-p))   
 y=rev(dat[(p+1):ndat])    
 fit=lm(y~-1+X, x=TRUE)
 
 return(fit)
 }
