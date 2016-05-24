
setMethod("summary", "ur.ers", function(object){
  return(new("sumurca", classname="ur.ers", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=object@model, type=object@type, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))  
})

setMethod("summary", "ca.jo", function(object){
  return(new("sumurca", classname="ca.jo", test.name=object@test.name, testreg=NULL, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=object@model, type=object@type, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=object@lambda, pval=NULL, V=object@V, W=object@W, P=object@P))
})

setMethod("summary", "cajo.test", function(object){
  return(new("sumurca", classname="cajo.test", test.name=object@test.name, testreg=NULL, teststat=object@teststat, cval=NULL, bpoint=NULL, signif=NULL, model=NULL, type=object@type, auxstat=NULL, lag=NULL, H=object@H, A=object@A, lambda=object@lambda, pval=object@pval, V=object@V, W=object@W, P=NULL))
})

setMethod("summary", "ur.kpss", function(object){
  return(new("sumurca", classname="ur.kpss", test.name=object@test.name, testreg=NULL, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=NULL, type=object@type, auxstat=NULL, lag=object@lag, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})

setMethod("summary", "ca.po", function(object){
  return(new("sumurca", classname="ca.po", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=object@model, type=object@type, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})

setMethod("summary", "ur.pp", function(object){
  return(new("sumurca", classname="ur.pp", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=object@model, type=object@type, auxstat=object@auxstat, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})

setMethod("summary", "ur.df", function(object){
  return(new("sumurca", classname="ur.df", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=NULL, model=object@model, type=NULL, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})

setMethod("summary", "ur.sp", function(object){
  return(new("sumurca", classname="ur.sp", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=NULL, signif=object@signif, model=NULL, type=NULL, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})

setMethod("summary", "ur.za", function(object){
  return(new("sumurca", classname="ur.za", test.name=object@test.name, testreg=object@testreg, teststat=object@teststat, cval=object@cval, bpoint=object@bpoint, signif=NULL, model=NULL, type=NULL, auxstat=NULL, lag=NULL, H=NULL, A=NULL, lambda=NULL, pval=NULL, V=NULL, W=NULL, P=NULL))
})
