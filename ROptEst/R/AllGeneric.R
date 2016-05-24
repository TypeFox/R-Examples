if(!isGeneric("getInfRobIC")){
    setGeneric("getInfRobIC", 
        function(L2deriv, risk, neighbor, ...) standardGeneric("getInfRobIC"))
}
if(!isGeneric("getFixRobIC")){
    setGeneric("getFixRobIC", 
        function(Distr, risk, neighbor, ...) standardGeneric("getFixRobIC"))
}
if(!isGeneric("getAsRisk")){
    setGeneric("getAsRisk", 
        function(risk, L2deriv, neighbor, biastype, ...) standardGeneric("getAsRisk"))
}
if(!isGeneric("getFiRisk")){
    setGeneric("getFiRisk", 
        function(risk, Distr, neighbor, ...) standardGeneric("getFiRisk"))
}
if(!isGeneric("getInfClip")){
    setGeneric("getInfClip", 
        function(clip, L2deriv, risk, neighbor, ...) standardGeneric("getInfClip"))
}
if(!isGeneric("getInfRad")){
    setGeneric("getInfRad",
        function(clip, L2deriv, risk, neighbor, ...) standardGeneric("getInfRad"))
}
if(!isGeneric("getFixClip")){
    setGeneric("getFixClip", 
        function(clip, Distr, risk, neighbor, ...) standardGeneric("getFixClip"))
}
if(!isGeneric("getInfGamma")){
    setGeneric("getInfGamma", 
        function(L2deriv, risk,  neighbor, biastype, ...) standardGeneric("getInfGamma"))
}
if(!isGeneric("getInfCent")){
    setGeneric("getInfCent", 
        function(L2deriv, neighbor, biastype, ...) standardGeneric("getInfCent"))
}
if(!isGeneric("getInfStand")){
    setGeneric("getInfStand", 
        function(L2deriv,  neighbor, biastype, ...) standardGeneric("getInfStand"))
}
if(!isGeneric("getInfV")){
    setGeneric("getInfV", 
        function(L2deriv,  neighbor, biastype, ...) standardGeneric("getInfV"))
}
if(!isGeneric("optIC")){
    setGeneric("optIC", function(model, risk,  ...) standardGeneric("optIC"))
}
if(!isGeneric("optRisk")){
    setGeneric("optRisk", function(model, risk,  ...) standardGeneric("optRisk"))
}
if(!isGeneric("radiusMinimaxIC")){
    setGeneric("radiusMinimaxIC", function(L2Fam, neighbor, risk, ...) 
            standardGeneric("radiusMinimaxIC"))
}
if(!isGeneric("getIneffDiff")){
    setGeneric("getIneffDiff", function(radius, L2Fam, neighbor, risk, ...) 
            standardGeneric("getIneffDiff"))
}
if(!isGeneric("leastFavorableRadius")){
    setGeneric("leastFavorableRadius", function(L2Fam, neighbor, risk, ...) 
            standardGeneric("leastFavorableRadius"))
}
if(!isGeneric("lowerCaseRadius")){
    setGeneric("lowerCaseRadius", function(L2Fam, neighbor, risk,  biastype, ...) 
    standardGeneric("lowerCaseRadius"))
}
if(!isGeneric("minmaxBias")){
    setGeneric("minmaxBias", 
        function(L2deriv, neighbor, biastype, ...) standardGeneric("minmaxBias"))
}
if(!isGeneric("getL1normL2deriv")){
    setGeneric("getL1normL2deriv", 
        function(L2deriv, ...) standardGeneric("getL1normL2deriv"))
}
if(!isGeneric("updateNorm")){
    setGeneric("updateNorm", function(normtype, ...) standardGeneric("updateNorm"))
}
if(!isGeneric("getModifyIC")){
    setGeneric("getModifyIC", function(L2FamIC, neighbor, risk, ...) standardGeneric("getModifyIC"))
}
if(!isGeneric("scaleUpdateIC")){
    setGeneric("scaleUpdateIC", function(neighbor, ...) standardGeneric("scaleUpdateIC"))
}
if(!isGeneric("eff")){
    setGeneric("eff", function(object) standardGeneric("eff"))
}
if(!isGeneric("get.asGRisk.fct")){
    setGeneric("get.asGRisk.fct", function(Risk) standardGeneric("get.asGRisk.fct"))
}
if(!isGeneric("getStartIC")){
    setGeneric("getStartIC", function(model, risk, ...) standardGeneric("getStartIC"))
}
