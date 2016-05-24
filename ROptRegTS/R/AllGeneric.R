if(!isGeneric("ErrorDistr")){
    setGeneric("ErrorDistr", function(object) standardGeneric("ErrorDistr"))
}
if(!isGeneric("ErrorSymm")){
    setGeneric("ErrorSymm", function(object) standardGeneric("ErrorSymm"))
}
if(!isGeneric("RegDistr")){
    setGeneric("RegDistr", function(object) standardGeneric("RegDistr"))
}
if(!isGeneric("RegSymm")){
    setGeneric("RegSymm", function(object) standardGeneric("RegSymm"))
}
if(!isGeneric("Regressor")){
    setGeneric("Regressor", function(object) standardGeneric("Regressor"))
}
if(!isGeneric("ErrorL2deriv")){
    setGeneric("ErrorL2deriv", function(object) standardGeneric("ErrorL2deriv"))
}
if(!isGeneric("ErrorL2derivSymm")){
    setGeneric("ErrorL2derivSymm", function(object) standardGeneric("ErrorL2derivSymm"))
}
if(!isGeneric("ErrorL2derivDistr")){
    setGeneric("ErrorL2derivDistr", function(object) standardGeneric("ErrorL2derivDistr"))
}
if(!isGeneric("ErrorL2derivDistrSymm")){
    setGeneric("ErrorL2derivDistrSymm", function(object) standardGeneric("ErrorL2derivDistrSymm"))
}
if(!isGeneric("radiusCurve")){ 
    setGeneric("radiusCurve", function(object) standardGeneric("radiusCurve"))
}
if(!isGeneric("neighborRadiusCurve")){
    setGeneric("neighborRadiusCurve", function(object) standardGeneric("neighborRadiusCurve"))
}
if(!isGeneric("neighborRadiusCurve<-")){
    setGeneric("neighborRadiusCurve<-", function(object, value) standardGeneric("neighborRadiusCurve<-"))
}
if(!isGeneric("getInfRobRegTypeIC")){
    setGeneric("getInfRobRegTypeIC", 
        function(ErrorL2deriv, Regressor, risk, neighbor, ...) standardGeneric("getInfRobRegTypeIC"))
}
if(!isGeneric("getFixRobRegTypeIC")){
    setGeneric("getFixRobRegTypeIC", 
        function(ErrorDistr, Regressor, risk, neighbor, ...) standardGeneric("getFixRobRegTypeIC"))
}
if(!isGeneric("getAsRiskRegTS")){
    setGeneric("getAsRiskRegTS", 
        function(risk, ErrorL2deriv, Regressor, neighbor, ...) standardGeneric("getAsRiskRegTS"))
}
if(!isGeneric("getFiRiskRegTS")){
    setGeneric("getFiRiskRegTS", 
        function(risk, ErrorDistr, Regressor, neighbor, ...) standardGeneric("getFiRiskRegTS"))
}
if(!isGeneric("getInfClipRegTS")){
    setGeneric("getInfClipRegTS", 
        function(clip, ErrorL2deriv, Regressor, risk, neighbor, ...) standardGeneric("getInfClipRegTS"))
}
if(!isGeneric("getFixClipRegTS")){
    setGeneric("getFixClipRegTS", 
        function(clip, ErrorDistr, Regressor, risk, neighbor, ...) standardGeneric("getFixClipRegTS"))
}
if(!isGeneric("getInfGammaRegTS")){
    setGeneric("getInfGammaRegTS", 
        function(ErrorL2deriv, Regressor, risk, neighbor, ...) standardGeneric("getInfGammaRegTS"))
}
if(!isGeneric("getInfCentRegTS")){
    setGeneric("getInfCentRegTS", 
        function(ErrorL2deriv, Regressor, neighbor, ...) standardGeneric("getInfCentRegTS"))
}
if(!isGeneric("getInfStandRegTS")){
    setGeneric("getInfStandRegTS", 
        function(ErrorL2deriv, Regressor, neighbor, ...) standardGeneric("getInfStandRegTS"))
}
