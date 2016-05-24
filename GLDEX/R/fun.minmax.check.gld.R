fun.minmax.check.gld <-
function(data,lambdas,param,lessequalmin=1,greaterequalmax=1){

# lessequalmin==1 means fitted GLD min <= min of data
# lessequalmin==0 means fitted GLD min < min of data
# greaterequalmin==1 means fitted GLD min >= max of data
# greaterequalmin==0 means fitted GLD min > max of data

min.d<-min(data)
max.d<-max(data)

result<-switch(param, FKML= , fkml = , freimer = , frm = , FMKL = , fmkl = .C("q_fmkl_gld_minmax_check", 
        as.double(min.d),as.double(max.d),as.integer(lessequalmin),
        as.integer(greaterequalmax),as.double(lambdas[,1,drop=F]), as.double(lambdas[,2,drop=F]), 
        as.double(lambdas[,3,drop=F]), as.double(lambdas[,4,drop=F]), as.integer(nrow(lambdas)),
        as.double(sqrt(.Machine$double.eps)), as.double((rep(0,nrow(lambdas)))))[[11]],
        ramberg = , ram = , RS = , rs = .C("q_rs_gld_minmax_check",  
        as.double(min.d),as.double(max.d),as.integer(lessequalmin),
        as.integer(greaterequalmax),as.double(lambdas[,1,drop=F]), as.double(lambdas[,2,drop=F]), 
        as.double(lambdas[,3,drop=F]), as.double(lambdas[,4,drop=F]), as.integer(nrow(lambdas)),
        as.double(sqrt(.Machine$double.eps)), as.double((rep(0,nrow(lambdas)))))[[11]])
        
return(as.logical(result))

}
