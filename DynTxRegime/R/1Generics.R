setGeneric(name = "Base", 
           def = function(object, ...){standardGeneric("Base")})

setGeneric(name= "Call",  
           def = function(object, ...){standardGeneric("Call")})

setGeneric(name = "Classif", 
           def = function(object, ...){standardGeneric("Classif")})

setGeneric(name = "Coef", 
           def = function(object, ...){standardGeneric("Coef")})

setGeneric(name = "DecisionPoint", 
           def = function(object, ...){standardGeneric("DecisionPoint")})

setGeneric(name = "Est", 
           def = function(x, ...){standardGeneric("Est")})

setGeneric(name = "Fit", 
           def = function(object, data, response, ...){standardGeneric("Fit")})

setGeneric(name = "FitObject", 
           def = function(object, ...){standardGeneric("FitObject")})

setGeneric(name = "Fitted",  
           def = function(object, ...){standardGeneric("Fitted")})
           
setGeneric(name = "FittedCont",  
           def = function(object, ...){standardGeneric("FittedCont")})
           
setGeneric(name = "FittedMain", 
           def = function(object, ...){standardGeneric("FittedMain")})
           
setGeneric(name = "FitType",  
           def = function(object, ...){standardGeneric("FitType")})

setGeneric(name = "Genetic", 
           def = function(object, ...){standardGeneric("Genetic")})

setGeneric(name = "IStep", 
           def = function(object, ...){standardGeneric("IStep")})

setGeneric(name = "ModelObject", 
           def = function(object, ...){standardGeneric("ModelObject")})

setGeneric(name = "ModelObjectFit", 
           def = function(object, ...){standardGeneric("ModelObjectFit")})

setGeneric(name = "ModelObjectFitCont", 
           def = function(object, ...){standardGeneric("ModelObjectFitCont")})

setGeneric(name = "ModelObjectFitMain", 
           def = function(object, ...){standardGeneric("ModelObjectFitMain")})

setGeneric(name = "MySummary", 
           def = function(object, ...){standardGeneric("MySummary")})

setGeneric(name = "NVars", 
           def = function(object, ...){standardGeneric("NVars")})

setGeneric(name = "OptTx",
           def = function(x, newdata, ...){standardGeneric("OptTx")})

setGeneric(name = "Outcome", 
           def = function(object, ...){standardGeneric("Outcome")})

setGeneric(name = "Plot", 
           def = function(x, ...){standardGeneric("Plot")})

setGeneric(name = "Predict", 
           def = function(object, newdata, ...){standardGeneric("Predict")})

setGeneric(name = "PredictCont",  
           def = function(object, newdata, ...){standardGeneric("PredictCont")})
           
setGeneric(name = "PredictMain", 
           def = function(object, newdata, ...){standardGeneric("PredictMain")})
           
setGeneric(name = "PredictPropen", 
           def = function(object, newdata, ...){standardGeneric("PredictPropen")})

setGeneric(name = "Print", 
           def = function(x, ...){standardGeneric("Print")})
           
setGeneric(name = "Propen",  
           def = function(object, ...){standardGeneric("Propen")})

setGeneric(name = "PtsSubset",  
           def = function(object, ...){standardGeneric("PtsSubset")})
           
setGeneric(name = "qFuncs", 
           def = function(object, ...){standardGeneric("qFuncs")})

setGeneric(name = "qqPlot", 
           def = function(x, ...){standardGeneric("qqPlot")})

setGeneric(name = "RegFunc", 
           def = function(object, ...){standardGeneric("RegFunc")})

setGeneric(name = "RegimeCoef", 
           def = function(object, ...){standardGeneric("RegimeCoef")})

setGeneric(name = "Regimes", 
           def = function(object, ...){standardGeneric("Regimes")})

setGeneric(name = "Residuals", 
           def = function(object, ...){standardGeneric("Residuals")})

setGeneric(name = "Scale", 
           def = function(object, ...){standardGeneric("Scale")})

setGeneric(name = "Show", 
           def = function(object, ...){standardGeneric("Show")})

setGeneric(name = "StdDev",
           def = function(object, ...){standardGeneric("StdDev")})

setGeneric(name = "Step", 
           def = function(object, ...){standardGeneric("Step")})

setGeneric(name = "Subset",  
           def = function(object, ...){standardGeneric("Subset")})

setGeneric(name = "SubsetRule", 
           def = function(object, ...){standardGeneric("SubsetRule")})

setGeneric(name = "Subsets",  
           def = function(object, ...){standardGeneric("Subsets")})

setGeneric(name = "SuperSet", 
           def = function(object, ...){standardGeneric("SuperSet")})
           
setGeneric(name = "TxName", 
           def = function(object, ...){standardGeneric("TxName")})
           
setGeneric(name = "TxInfo", 
           def = function(object, ...){standardGeneric("TxInfo")})

setGeneric(name = "TxVec", 
           def = function(object, ...){standardGeneric("TxVec")})

setGeneric(name = "VNames", 
           def = function(object, ...){standardGeneric("VNames")})

setGeneric(name = "YTilde",
           def = function(object, ...){standardGeneric("YTilde")})

setGeneric(name = "YTilde<-",
           def = function(object, value){standardGeneric("YTilde<-")})


if(!isGeneric("classif")){
  setGeneric(name = "classif", 
             def = function(object, ...){standardGeneric("classif")})
}

if(!isGeneric("DTRstep")){
  setGeneric(name = "DTRstep", 
             def = function(object, ...){standardGeneric("DTRstep")})
}

if(!isGeneric("estimator")){
  setGeneric(name = "estimator", 
             def = function(x, ...){standardGeneric("estimator")})
}

if(!isGeneric("fittedCont")){
  setGeneric(name = "fittedCont",  
             def = function(object, ...){standardGeneric("fittedCont")})
}
           
if(!isGeneric("fittedMain")){
  setGeneric(name = "fittedMain", 
             def = function(object, ...){standardGeneric("fittedMain")})
}
           
if(!isGeneric("FC")){
  setGeneric(name = "FC",  
             def = function(object, ...){standardGeneric("FC")})
}
           
if(!isGeneric("FM")){
  setGeneric(name = "FM", 
             def = function(object, ...){standardGeneric("FM")})
}
           
if(!isGeneric("genetic")){
  setGeneric(name = "genetic", 
             def = function(object, ...){standardGeneric("genetic")})
}

if(!isGeneric("optTx")){
  setGeneric(name = "optTx", 
             def = function(x, newdata, ...){standardGeneric("optTx")})
}

if(!isGeneric("outcome")){
  setGeneric(name = "outcome", 
             def = function(object, ...){standardGeneric("outcome")})
}

if(!isGeneric("propen")){
  setGeneric(name = "propen", 
             def = function(object, ...){standardGeneric("propen")})
}

if(!isGeneric("regimeCoef")){
  setGeneric(name = "regimeCoef", 
             def = function(object, ...){standardGeneric("regimeCoef")})
}

if(!isGeneric("stdDev")){
  setGeneric(name = "stdDev", 
             def = function(object, ...){standardGeneric("stdDev")})
}

