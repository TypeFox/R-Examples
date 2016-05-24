



setGeneric("hbstm",function(Zt,K,newGrid,reglag,seas,spatlags,hyperpriors,initialvalues,nIter,nBurn,fit,plots,posterior,save,control){standardGeneric("hbstm")})

setGeneric("hbstm.fit",function(HBSTM,nIter,nBurn,time,timerem,plots,posterior,save){standardGeneric("hbstm.fit")})

setGeneric("plotRes",function(object,point,ARlags,ARperiod){standardGeneric("plotRes")})

setGeneric("plotFit",function(object,time,values){standardGeneric("plotFit")})

setGeneric("results",function(object,spatTemp,inter,digits,component,plots,file){standardGeneric("results")})

setGeneric("estimation",function(object,digits){standardGeneric("estimation")})

setGeneric("mse",function(object){standardGeneric("mse")})
