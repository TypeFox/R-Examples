########################################################################
### intentional maskings --- purpose: add formal arg "..." to generic...
########################################################################

############################################################################
### ------------------------------------------------
### In this comment substitute 'xxx' by 
###  'var', 'median', 'IQR', 'mad', respectively
### ------------------------------------------------
### We intentionally mask function 'xxx' from stats in order to add a formal
### argument "...". 
### functionality of 'stats::xxx' is completely retained, however
### for help to the original 'stats::xxx' function write
###       'help("xxx", package="stats")'                             
### for code to the original 'stats::xxx' function write
###       'stats::xxx'
############################################################################

## included from version 1.7 of distrEx on

### intentionally mask functionals for additional ... argument P.R. 28-03-06

var <- function(x , ...)
       {dots <- list(...)
        if(hasArg(y)) y <- dots$"y"
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        if(!hasArg(use)) 
             use <- ifelse (na.rm, "complete.obs","all.obs")
        else use <- dots$"use"
        if(hasArg(y))       
           stats::var(x = x, y = y, na.rm = na.rm, use)
        else
           stats::var(x = x, y = NULL, na.rm = na.rm, use)
        }   

## sd already masked in NormalDistribution.R in package "distr"

median <- function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        stats::median(x = x, na.rm = na.rm)}

IQR <- function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        stats::IQR(x = x, na.rm = na.rm)}

mad <- function(x , ...)
       {dots <- list(...)
        na.rm     <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        low       <-  ifelse(hasArg(low), dots$"low", FALSE)
        high      <-  ifelse(hasArg(high), dots$"high", FALSE)
        center    <-  ifelse(hasArg(center), dots$"center", median(x))
        constant  <-  ifelse(hasArg(constant), dots$"constant", 1.4826)
        stats::mad(x = x, center = center, constant = constant , na.rm = na.rm, 
                   low = low, high = high)}

### --------- registration as generics ------------------

if(!isGeneric("var")){ 
   setGeneric("var", function(x, ...) standardGeneric("var"))
}

##sd already registered as generic in package "distr"

if(!isGeneric("median")){ 
   setGeneric("median", function(x, ...) standardGeneric("median"))
}

if(!isGeneric("IQR")){ 
   setGeneric("IQR", function(x, ...) standardGeneric("IQR"))
}

if(!isGeneric("mad")){ 
   setGeneric("mad", function(x, ...) standardGeneric("mad"))
}


if(!isGeneric("skewness")){ 
   setGeneric("skewness", function(x, ...) standardGeneric("skewness"))
}

if(!isGeneric("kurtosis")){ 
   setGeneric("kurtosis", function(x, ...) standardGeneric("kurtosis"))
}

############################################################################
# Access methods
############################################################################

if(!isGeneric("cond")){
   setGeneric("cond", function(object) standardGeneric("cond"))
}

if(!isGeneric("Range")){ 
   setGeneric("Range", function(object, ...) standardGeneric("Range"))
}



############################################################################
# generics to  "usual"  methods
############################################################################

if(!isGeneric("ContaminationSize")){ 
   setGeneric("ContaminationSize", 
               function(e1, e2, ...) standardGeneric("ContaminationSize"))
}
if(!isGeneric("TotalVarDist")){
   setGeneric("TotalVarDist", 
               function(e1, e2, ...) standardGeneric("TotalVarDist"))
}
if(!isGeneric("KolmogorovDist")){
   setGeneric("KolmogorovDist", 
               function(e1, e2, ...) standardGeneric("KolmogorovDist"))
}
if(!isGeneric("HellingerDist")){
   setGeneric("HellingerDist", 
               function(e1, e2, ...) standardGeneric("HellingerDist"))
}
if(!isGeneric("CvMDist")){
   setGeneric("CvMDist", 
               function(e1, e2, ...) standardGeneric("CvMDist"))
}

if(!isGeneric("ConvexContamination")){ 
   setGeneric("ConvexContamination", 
               function(e1, e2, size) standardGeneric("ConvexContamination"))
}

if(!isGeneric("AsymTotalVarDist")){
   setGeneric("AsymTotalVarDist", 
               function(e1, e2, ...) standardGeneric("AsymTotalVarDist"))
}
if(!isGeneric("OAsymTotalVarDist")){
   setGeneric("OAsymTotalVarDist", 
               function(e1, e2, ...) standardGeneric("OAsymTotalVarDist"))
}

if(!isGeneric("E")){ 
   setGeneric("E", function(object, fun, cond, ...) standardGeneric("E"))
}

if(!isGeneric("m1df")){
   setGeneric("m1df", function(object, upper, ...) standardGeneric("m1df"))
}
if(!isGeneric("m2df")){
   setGeneric("m2df", function(object, upper, ...) standardGeneric("m2df"))
}

#if(!isGeneric("illustrateCLT")){
#   setGeneric("illustrateCLT", function(Distr, ...) 
#      standardGeneric("illustrateCLT"))
#}
#if(!isGeneric("plotCLT")){
#  setGeneric("plotCLT", function(Tn, ...) standardGeneric("plotCLT"))
#}

