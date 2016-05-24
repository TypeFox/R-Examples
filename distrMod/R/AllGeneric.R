
### ---------------------------
### "multiple-purpose generics" --- 
###   for  
###      + accessor method  
###      + stats/base function
###      + functional method
### ---------------------------

### intentionally mask confint for additional ... argument P.R. 28-03-06

confint <- function(object, method, ...)
       {dots <- list(...)
        mc <- match.call()[-1]
        nmc <- names(as.list(mc))
        mc0 <- as.list(mc)[!nmc %in% c("object", "parm", "level")]
        argList <- list(object = object)
        parm <- level <- NULL
        if(hasArg(parm)) argList <- c(argList, parm = dots$"parm")
        if(hasArg(level)) argList <- c(argList, level = dots$"level")
        else argList <- c(argList, level = 0.95)   
        if(length(mc0)) argList <- c(argList, mc0)
        do.call(stats::confint, argList)
        }   


## access and replace methods
#if(!isGeneric("param<-")){ 
#    setGeneric("param<-", function(object, value) standardGeneric("param<-"))
#}
if(!isGeneric("distrSymm")){
    setGeneric("distrSymm", function(object) standardGeneric("distrSymm"))
}
if(!isGeneric("props")){
    setGeneric("props", function(object) standardGeneric("props"))
}
if(!isGeneric("props<-")){
    setGeneric("props<-", function(object, value) standardGeneric("props<-"))
}
if(!isGeneric("addProp<-")){
    setGeneric("addProp<-", function(object, value) standardGeneric("addProp<-"))
}
if(!isGeneric("main")){
    setGeneric("main", function(object) standardGeneric("main"))
}
if(!isGeneric("main<-")){
    setGeneric("main<-", function(object, value) standardGeneric("main<-"))
}
if(!isGeneric("nuisance")){
    setGeneric("nuisance", function(object) standardGeneric("nuisance"))
}
if(!isGeneric("nuisance<-")){
    setGeneric("nuisance<-", function(object, value) standardGeneric("nuisance<-"))
}
if(!isGeneric("trafo")){
    setGeneric("trafo", function(object, param, ...) standardGeneric("trafo"))
}
if(!isGeneric("trafo<-")){
    setGeneric("trafo<-", function(object, value) standardGeneric("trafo<-"))
}
#new 20.2.09
if(!isGeneric("trafo.fct")){
    setGeneric("trafo.fct", function(object) standardGeneric("trafo.fct"))
}
if(!isGeneric("modifyParam")){
    setGeneric("modifyParam", function(object) standardGeneric("modifyParam"))
}
if(!isGeneric("fam.call")){
    setGeneric("fam.call", function(object) standardGeneric("fam.call"))
}
if(!isGeneric("dimension")){
    setGeneric("dimension", function(object) standardGeneric("dimension"))
}
if(!isGeneric("L2deriv")){
    setGeneric("L2deriv", function(object, param) standardGeneric("L2deriv"))
}
if(!isGeneric("L2derivSymm")){
    setGeneric("L2derivSymm", function(object) standardGeneric("L2derivSymm"))
}
if(!isGeneric("L2derivDistr")){
    setGeneric("L2derivDistr", function(object) standardGeneric("L2derivDistr"))
}
if(!isGeneric("L2derivDistrSymm")){
    setGeneric("L2derivDistrSymm", function(object) standardGeneric("L2derivDistrSymm"))
}
if(!isGeneric("FisherInfo")){
    setGeneric("FisherInfo", function(object, param) standardGeneric("FisherInfo"))
}
if(!isGeneric("checkL2deriv")){
    setGeneric("checkL2deriv", function(L2Fam, ...) standardGeneric("checkL2deriv"))
}
if(!isGeneric("bound")){ 
    setGeneric("bound", function(object) standardGeneric("bound"))
}
if(!isGeneric("width")){ 
    setGeneric("width", function(object) standardGeneric("width"))
}

if(!isGeneric("sign")){
    setGeneric("sign", function(x) standardGeneric("sign"))
}
if(!isGeneric("nu")){
    setGeneric("nu", function(object) standardGeneric("nu"))
}

if(!isGeneric("sign<-")){
    setGeneric("sign<-", function(object,value) standardGeneric("sign<-"))
}
if(!isGeneric("nu<-")){
    setGeneric("nu<-", function(object,value) standardGeneric("nu<-"))
}

if(!isGeneric("biastype")){
    setGeneric("biastype", function(object) standardGeneric("biastype"))
}

if(!isGeneric("biastype<-")){
    setGeneric("biastype<-", function(object,value) standardGeneric("biastype<-"))
}

if(!isGeneric("solve")){
    setGeneric("solve", function(a,b,...) standardGeneric("solve"))
}

if(!isGeneric("modifyModel")){
    setGeneric("modifyModel", function(model, param, ...) standardGeneric("modifyModel"))
}
if(!isGeneric("existsPIC")){
    setGeneric("existsPIC", function(object,...) standardGeneric("existsPIC"))
}

if(!isGeneric("norm")){
    setGeneric("norm", function(x,...) standardGeneric("norm"))
}

if(!isGeneric("normtype")){
    setGeneric("normtype", function(object) standardGeneric("normtype"))
}

if(!isGeneric("normtype<-")){
    setGeneric("normtype<-", function(object,value) standardGeneric("normtype<-"))
}

if(!isGeneric("QuadForm")){
    setGeneric("QuadForm", function(object) standardGeneric("QuadForm"))
}
if(!isGeneric("QuadForm<-")){
    setGeneric("QuadForm<-", function(object,value) standardGeneric("QuadForm<-"))
}

if(!isGeneric("fct<-")){
    setGeneric("fct<-", function(object,value) standardGeneric("fct<-"))
}
if(!isGeneric("fct")){
    setGeneric("fct", function(object) standardGeneric("fct"))
}
if(!isGeneric("estimate")){
    setGeneric("estimate", function(object) standardGeneric("estimate"))
}
if(!isGeneric("estimate.call")){
    setGeneric("estimate.call", function(object) standardGeneric("estimate.call"))
}
if(!isGeneric("name.estimate")){
    setGeneric("name.estimate", function(object) standardGeneric("name.estimate"))
}
if(!isGeneric("trafo.estimate")){
    setGeneric("trafo.estimate", function(object) standardGeneric("trafo.estimate"))
}
if(!isGeneric("nuisance.estimate")){
    setGeneric("nuisance.estimate", function(object) standardGeneric("nuisance.estimate"))
}
if(!isGeneric("samplesize.estimate")){
    setGeneric("samplesize.estimate", function(object, ...) standardGeneric("samplesize.estimate"))
}
if(!isGeneric("call.estimate")){
    setGeneric("call.estimate", function(object) standardGeneric("call.estimate"))
}
if(!isGeneric("fixed.estimate")){
    setGeneric("fixed.estimate", function(object,... ) standardGeneric("fixed.estimate"))
}
if(!isGeneric("completecases.estimate")){
  	     setGeneric("completecases.estimate", function(object) standardGeneric("completecases.estimate"))
}
if(!isGeneric("Infos")){
    setGeneric("Infos", function(object) standardGeneric("Infos"))
}
if(!isGeneric("Infos<-")){
    setGeneric("Infos<-", function(object, value) standardGeneric("Infos<-"))
}
if(!isGeneric("addInfo<-")){
    setGeneric("addInfo<-", function(object, value) standardGeneric("addInfo<-"))
}
if(!isGeneric("criterion")){
    setGeneric("criterion", function(object) standardGeneric("criterion"))
}
if(!isGeneric("criterion<-")){
    setGeneric("criterion<-", function(object, value) standardGeneric("criterion<-"))
}
if(!isGeneric("completecases")){
  	     setGeneric("completecases", function(object) standardGeneric("completecases"))
}
if(!isGeneric("asvar")){
    setGeneric("asvar", function(object) standardGeneric("asvar"))
}
if(!isGeneric("asvar<-")){
    setGeneric("asvar<-", function(object, value) standardGeneric("asvar<-"))
}
if(!isGeneric("nuisance")){
    setGeneric("nuisance", function(object,... ) standardGeneric("nuisance"))
}
if(!isGeneric("main")){
    setGeneric("main", function(object,... ) standardGeneric("main"))
}
if(!isGeneric("fixed")){
    setGeneric("fixed", function(object,... ) standardGeneric("fixed"))
}
if(!isGeneric("nuisance<-")){
    setGeneric("nuisance<-", function(object, value) standardGeneric("nuisance<-"))
}
if(!isGeneric("main<-")){
    setGeneric("main<-", function(object, value) standardGeneric("main<-"))
}
if(!isGeneric("fixed<-")){
    setGeneric("fixed<-", function(object, value) standardGeneric("fixed<-"))
}

#if(!isGeneric("confint")){
setGeneric("confint", function(object, method, ... ) standardGeneric("confint"))
#}

if(!isGeneric("validParameter")){
    setGeneric("validParameter", function(object, ... ) standardGeneric("validParameter"))
}
if(!isGeneric("untransformed.asvar")){
    setGeneric("untransformed.asvar", 
                function(object) standardGeneric("untransformed.asvar"))
}
if(!isGeneric("untransformed.estimate")){
    setGeneric("untransformed.estimate", 
                function(object) standardGeneric("untransformed.estimate"))
}
if(!isGeneric("startPar")){
    setGeneric("startPar", function(object, ...) standardGeneric("startPar"))
}
if(!isGeneric("makeOKPar")){
    setGeneric("makeOKPar", function(object, ...) standardGeneric("makeOKPar"))
}              
if(!isGeneric("locscalename")){
    setGeneric("locscalename", function(object) standardGeneric("locscalename"))
}
if(!isGeneric("locscalename<-")){
    setGeneric("locscalename<-", function(object, value) standardGeneric("locscalename<-"))
}
if(!isGeneric("LogDeriv")){
    setGeneric("LogDeriv", function(object) standardGeneric("LogDeriv"))
}
if(!isGeneric("LogDeriv<-")){
    setGeneric("LogDeriv<-", function(object, value) standardGeneric("LogDeriv<-"))
}
if(!isGeneric("mleCalc")){
    setGeneric("mleCalc", function(x, PFam, ...) standardGeneric("mleCalc"))
}
if(!isGeneric("mceCalc")){
    setGeneric("mceCalc", function(x, PFam, ...) standardGeneric("mceCalc"))
}
if(!isGeneric("criterion.fct")){
    setGeneric("criterion.fct", function(object) standardGeneric("criterion.fct"))
}
if(!isGeneric("method")){
    setGeneric("method", function(object) standardGeneric("method"))
}
if(!isGeneric("optimwarn")){
    setGeneric("optimwarn", function(object) standardGeneric("optimwarn"))
}
if(!isGeneric("withPosRestr")){
    setGeneric("withPosRestr", function(object) standardGeneric("withPosRestr"))
}
if(!isGeneric("withPosRestr<-")){
    setGeneric("withPosRestr<-", function(object,value) standardGeneric("withPosRestr<-"))
}
if(!isGeneric("scaleshapename")){
    setGeneric("scaleshapename", function(object) standardGeneric("scaleshapename"))
}
if(!isGeneric("scaleshapename<-")){
    setGeneric("scaleshapename<-", function(object, value) standardGeneric("scaleshapename<-"))
}
if(!isGeneric("scalename")){
    setGeneric("scalename", function(object) standardGeneric("scalename"))
}
