if(!isGeneric("radius")){ 
    setGeneric("radius", function(object) standardGeneric("radius"))
}
if(!isGeneric("radius<-")){
    setGeneric("radius<-", function(object,value) standardGeneric("radius<-"))
}
if(!isGeneric("center")){
    setGeneric("center", function(object) standardGeneric("center"))
}
if(!isGeneric("center<-")){
    setGeneric("center<-", function(object, value) standardGeneric("center<-"))
}
if(!isGeneric("neighbor")){
    setGeneric("neighbor", function(object) standardGeneric("neighbor"))
}
if(!isGeneric("neighbor<-")){
    setGeneric("neighbor<-", function(object, value) standardGeneric("neighbor<-"))
}
if(!isGeneric("Curve")){
    setGeneric("Curve", function(object) standardGeneric("Curve"))
}
if(!isGeneric("Risks")){
    setGeneric("Risks", function(object) standardGeneric("Risks"))
}
if(!isGeneric("Risks<-")){
    setGeneric("Risks<-", function(object, value) standardGeneric("Risks<-"))
}
if(!isGeneric("addRisk<-")){
    setGeneric("addRisk<-", function(object, value) standardGeneric("addRisk<-"))
}
if(!isGeneric("CallL2Fam")){ 
    setGeneric("CallL2Fam", function(object) standardGeneric("CallL2Fam"))
}
if(!isGeneric("CallL2Fam<-")){ 
    setGeneric("CallL2Fam<-", function(object, value) standardGeneric("CallL2Fam<-"))
}
if(!isGeneric("generateIC")){
    setGeneric("generateIC", function(neighbor, L2Fam, ...) standardGeneric("generateIC"))
}
if(!isGeneric("checkIC")){
    setGeneric("checkIC", function(IC, L2Fam, ...) standardGeneric("checkIC"))
}
if(!isGeneric("evalIC")){
    setGeneric("evalIC", function(IC, x) standardGeneric("evalIC"))
}
if(!isGeneric("makeIC")){
    setGeneric("makeIC", function(IC, L2Fam, ...) standardGeneric("makeIC"))
}
if(!isGeneric("clip")){
    setGeneric("clip", function(x1, ...) standardGeneric("clip"))
}
if(!isGeneric("clip<-")){
    setGeneric("clip<-", function(object, value) standardGeneric("clip<-"))
}
if(!isGeneric("cent")){
    setGeneric("cent", function(object) standardGeneric("cent"))
}
if(!isGeneric("cent<-")){
    setGeneric("cent<-", function(object, value) standardGeneric("cent<-"))
}
if(!isGeneric("stand")){
    setGeneric("stand", function(object) standardGeneric("stand"))
}
if(!isGeneric("stand<-")){
    setGeneric("stand<-", function(object, value) standardGeneric("stand<-"))
}
if(!isGeneric("lowerCase")){
    setGeneric("lowerCase", function(object) standardGeneric("lowerCase"))
}
if(!isGeneric("lowerCase<-")){
    setGeneric("lowerCase<-", function(object, value) standardGeneric("lowerCase<-"))
}
if(!isGeneric("neighborRadius")){
    setGeneric("neighborRadius", function(object) standardGeneric("neighborRadius"))
}
if(!isGeneric("neighborRadius<-")){
    setGeneric("neighborRadius<-", function(object, value) standardGeneric("neighborRadius<-"))
}
if(!isGeneric("clipLo")){
    setGeneric("clipLo", function(object) standardGeneric("clipLo"))
}
if(!isGeneric("clipLo<-")){
    setGeneric("clipLo<-", function(object, value) standardGeneric("clipLo<-"))
}
if(!isGeneric("clipUp")){
    setGeneric("clipUp", function(object) standardGeneric("clipUp"))
}
if(!isGeneric("clipUp<-")){
    setGeneric("clipUp<-", function(object, value) standardGeneric("clipUp<-"))
}
#if(!isGeneric("oneStepEstimator")){
#    setGeneric("oneStepEstimator",
#        function(x, IC, start, ...) standardGeneric("oneStepEstimator"))
#}
#if(!isGeneric("kStepEstimator")){
#    setGeneric("kStepEstimator",
#        function(x, IC, start, ...) standardGeneric("kStepEstimator"))
#}
if(!isGeneric("locMEstimator")){
    setGeneric("locMEstimator", function(x, IC, ...) standardGeneric("locMEstimator"))
}
if(!isGeneric("infoPlot")){
    setGeneric("infoPlot", function(object,...) standardGeneric("infoPlot"))
}
if(!isGeneric("optIC")){
    setGeneric("optIC", function(model, risk, ...) standardGeneric("optIC"))
}


if(!isGeneric("weight")){
    setGeneric("weight",
        function(object, ...) standardGeneric("weight"))
}
if(!isGeneric("weight<-")){
    setGeneric("weight<-",
        function(object, value) standardGeneric("weight<-"))
}
if(!isGeneric("clip<-")){
    setGeneric("clip<-",
        function(object, value, ...) standardGeneric("clip<-"))
}
if(!isGeneric("stand")){
    setGeneric("stand",
        function(object, ...) standardGeneric("stand"))
}
if(!isGeneric("stand<-")){
    setGeneric("stand<-",
        function(object, value, ...) standardGeneric("stand<-"))
}
if(!isGeneric("cent")){
    setGeneric("cent",
        function(object, ...) standardGeneric("cent"))
}
if(!isGeneric("cent<-")){
    setGeneric("cent<-",
        function(object, value, ...) standardGeneric("cent<-"))
}

if(!isGeneric("getweight")){
    setGeneric("getweight",
        function(Weight, neighbor, biastype, ...) standardGeneric("getweight"))
}

if(!isGeneric("minbiasweight")){
    setGeneric("minbiasweight",
        function(Weight, neighbor, biastype, ...) standardGeneric("minbiasweight"))
}
if(!isGeneric("generateIC.fct")){
    setGeneric("generateIC.fct", function(neighbor, L2Fam, ...) standardGeneric("generateIC.fct"))
}
if(!isGeneric("getRiskIC")){
    setGeneric("getRiskIC", 
        function(IC, risk,  neighbor, L2Fam, ...) standardGeneric("getRiskIC"))
}
if(!isGeneric("getBiasIC")){
    setGeneric("getBiasIC", 
        function(IC, neighbor, ...) standardGeneric("getBiasIC"))
}
if(!isGeneric(".evalBiasIC")){
    setGeneric(".evalBiasIC", 
        function(IC, neighbor, biastype, ...) standardGeneric(".evalBiasIC"))
}
if(!isGeneric("comparePlot")){
    setGeneric("comparePlot", function(obj1,obj2,...) standardGeneric("comparePlot"))
}
if(!isGeneric("pIC")){
    setGeneric("pIC", function(object) standardGeneric("pIC"))
}
if(!isGeneric("asbias")){
    setGeneric("asbias", function(object) standardGeneric("asbias"))
}
if(!isGeneric("steps")){
    setGeneric("steps", function(object) standardGeneric("steps"))
}
if(!isGeneric("ksteps")){
    setGeneric("ksteps", function(object,...) standardGeneric("ksteps"))
}
if(!isGeneric("uksteps")){
    setGeneric("uksteps", function(object,...) standardGeneric("uksteps"))
}
if(!isGeneric("start")){
    setGeneric("start", function(x, ...) standardGeneric("start"))
}
if(!isGeneric("startval")){
    setGeneric("startval", function(object) standardGeneric("startval"))
}
if(!isGeneric("ustartval")){
    setGeneric("ustartval", function(object) standardGeneric("ustartval"))
}
if(!isGeneric("ICList")){
    setGeneric("ICList", function(object) standardGeneric("ICList"))
}
if(!isGeneric("pICList")){
    setGeneric("pICList", function(object) standardGeneric("pICList"))
}
if(!isGeneric("Mroot")){
    setGeneric("Mroot", function(object) standardGeneric("Mroot"))
}
if(!isGeneric("modifyIC")){
    setGeneric("modifyIC", function(object) standardGeneric("modifyIC"))
}
if(!isGeneric("cutoff.quantile")){
    setGeneric("cutoff.quantile", function(object) standardGeneric("cutoff.quantile"))
}
if(!isGeneric("cutoff.quantile<-")){
    setGeneric("cutoff.quantile<-", function(object,value)
                standardGeneric("cutoff.quantile<-"))
}
if(!isGeneric("ddPlot")){
    setGeneric("ddPlot", function(data, dist.x, dist.y, cutoff.x, cutoff.y,...)
                                  standardGeneric("ddPlot"))
}
if(!isGeneric("kStepEstimator.start")){
    setGeneric("kStepEstimator.start",
                function(start,...) standardGeneric("kStepEstimator.start"))
}
if(!isGeneric("radius")){
    setGeneric("radius", function(object) standardGeneric("radius"))
}

if(!isGeneric("samplesize<-")){
    setGeneric("samplesize<-",
        function(object, value) standardGeneric("samplesize<-"))
}
if(!isGeneric("getRiskFctBV")){
    setGeneric("getRiskFctBV", function(risk, biastype) standardGeneric("getRiskFctBV"))
}

if(!isGeneric("moveL2Fam2RefParam")){
    setGeneric("moveL2Fam2RefParam", function(L2Fam, ...)
                standardGeneric("moveL2Fam2RefParam"))
}

if(!isGeneric("moveICBackFromRefParam")){
    setGeneric("moveICBackFromRefParam", function(IC, L2Fam, ...)
               standardGeneric("moveICBackFromRefParam"))
}

if(!isGeneric("rescaleFunction")){
    setGeneric("rescaleFunction", function(L2Fam, ...)
               standardGeneric("rescaleFunction"))
}
