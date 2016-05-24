
###_ + Plotting
if (!isGeneric("plotTDR")) {
    setGeneric("plotTDR",
               function(x, y, ...) standardGeneric("plotTDR"))
}

if (!isGeneric("plotZOC")) {
    setGeneric("plotZOC",
               function(x, y, ...) standardGeneric("plotZOC"))
}

if (!isGeneric("plotDiveModel")) {
    setGeneric("plotDiveModel",
               function(x, y, ...) standardGeneric("plotDiveModel"))
}

if (!isGeneric("plotBouts")) {
    setGeneric("plotBouts", function(fit, ...) standardGeneric("plotBouts"))
}

###_ + Accessors
if (!isGeneric("getFileName")) {        # File name accessor
    setGeneric("getFileName", function(x) standardGeneric("getFileName"))
}

if (!isGeneric("getTime")) {            # time accessor
    setGeneric("getTime", function(x) standardGeneric("getTime"))
}

if (!isGeneric("getDepth")) {           # Depth accessor
    setGeneric("getDepth", function(x) standardGeneric("getDepth"))
}

if (!isGeneric("getSpeed")) {           # speed accessor
    setGeneric("getSpeed", function(x) standardGeneric("getSpeed"))
}

if (!isGeneric("getDtime")) {           # interval accessor
    setGeneric("getDtime", function(x) standardGeneric("getDtime"))
}

if (!isGeneric("getCCData")) {          # concurrent data accessor
    setGeneric("getCCData", function(x, y) standardGeneric("getCCData"))
}

if (!isGeneric("getTDR")) {             # zoc'ed TDR accessor
    setGeneric("getTDR", function(x) standardGeneric("getTDR"))
}

if (!isGeneric("getGAct")) {            # gross activity accessor
    setGeneric("getGAct", function(x, y) standardGeneric("getGAct"))
}

if (!isGeneric("getDAct")) {            # dive activity accessor
    setGeneric("getDAct", function(x, y) standardGeneric("getDAct"))
}

if (!isGeneric("getDPhaseLab")) {       # dive phase label accessor
    setGeneric("getDPhaseLab",
               function(x, diveNo) standardGeneric("getDPhaseLab"))
}

if (!isGeneric("getDiveModel")) {     # dive model accessor
    setGeneric("getDiveModel",
               function(x, diveNo) standardGeneric("getDiveModel"))
}

if (!isGeneric("getDiveDeriv")) {     # dive derivative accessor
    setGeneric("getDiveDeriv",
               function(x, ...) standardGeneric("getDiveDeriv"))
}

if (!isGeneric("getSpeedCoef")) {       # speed calibration coefs accessor
    setGeneric("getSpeedCoef", function(x) standardGeneric("getSpeedCoef"))
}

###_ + Coercions and Replacements
if (!isGeneric("as.TDRspeed")) {        # coerce to TDRspeed
    setGeneric("as.TDRspeed", function(x) standardGeneric("as.TDRspeed"))
}

if (!isGeneric("depth<-")) {            # depth replacement
    setGeneric("depth<-", function(x, value) standardGeneric("depth<-"))
}

if (!isGeneric("speed<-")) {            # speed replacement
    setGeneric("speed<-", function(x, value) standardGeneric("speed<-"))
}

if (!isGeneric("ccData<-")) {           # concurrent data replacement
    setGeneric("ccData<-", function(x, value) standardGeneric("ccData<-"))
}

###_ + Generators and Summaries
if (!isGeneric("extractDive")) {        # extract a dive
    setGeneric("extractDive",
               function(obj, diveNo, id) standardGeneric("extractDive"))
}

if (!isGeneric("timeBudget")) {
    setGeneric("timeBudget",
               function(obj, ignoreZ) standardGeneric("timeBudget"))
}

###_ + bec2, and bec3
if (!isGeneric("bec2")) {
    setGeneric("bec2", function(fit) standardGeneric("bec2"))
}

if (!isGeneric("bec3")) {
    setGeneric("bec3", function(fit) standardGeneric("bec3"))
}


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
