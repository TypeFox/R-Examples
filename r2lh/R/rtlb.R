######################################################################
###                           Dispatching                          ###
###                               S4                               ###
######################################################################

### Ces deux classes sont définies dans "functions.R" parce qu'univ en a aussi besoin.
#setClass(Class="continuous",contains=c("numeric"))
#setClass(Class="discrete",contains=c("numeric"))

setGeneric(name="r2lBiv",def=function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){standardGeneric("r2lBiv")})

setMethod(f="r2lBiv",signature=c("logical","logical"),def=r2lBivLogicalLogical)
setMethod(f="r2lBiv",signature=c("logical","factor"),def=r2lBivLogicalFactor)
setMethod(f="r2lBiv",signature=c("logical","ordered"),def=r2lBivLogicalOrdered)
setMethod(f="r2lBiv",signature=c("logical","discrete"),def=r2lBivLogicalDiscrete)
setMethod(f="r2lBiv",signature=c("logical","continuous"),def=r2lBivLogicalContinuous)

setMethod(f="r2lBiv",signature=c("factor","logical"),def=r2lBivFactorLogical)
setMethod(f="r2lBiv",signature=c("factor","factor"),def=r2lBivFactorFactor)
setMethod(f="r2lBiv",signature=c("factor","ordered"),def=r2lBivFactorOrdered)
setMethod(f="r2lBiv",signature=c("factor","discrete"),def=r2lBivFactorDiscrete)
setMethod(f="r2lBiv",signature=c("factor","continuous"),def=r2lBivFactorContinuous)

setMethod(f="r2lBiv",signature=c("ordered","logical"),def=r2lBivOrderedLogical)
setMethod(f="r2lBiv",signature=c("ordered","factor"),def=r2lBivOrderedFactor)
setMethod(f="r2lBiv",signature=c("ordered","ordered"),def=r2lBivOrderedOrdered)
setMethod(f="r2lBiv",signature=c("ordered","discrete"),def=r2lBivOrderedDiscrete)
setMethod(f="r2lBiv",signature=c("ordered","continuous"),def=r2lBivOrderedContinuous)

setMethod(f="r2lBiv",signature=c("discrete","logical"),def=r2lBivDiscreteLogical)
setMethod(f="r2lBiv",signature=c("discrete","factor"),def=r2lBivDiscreteFactor)
setMethod(f="r2lBiv",signature=c("discrete","ordered"),def=r2lBivDiscreteOrdered)
setMethod(f="r2lBiv",signature=c("discrete","discrete"),def=r2lBivDiscreteDiscrete)
setMethod(f="r2lBiv",signature=c("discrete","continuous"),def=r2lBivDiscreteContinuous)

setMethod(f="r2lBiv",signature=c("continuous","logical"),def=r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("continuous","factor"),def=r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("continuous","ordered"),def=r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("continuous","discrete"),def=r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("continuous","continuous"),def=r2lBivContinuousContinuous)


#setGeneric(name="r2lBiv",def=function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){standardGeneric("r2lBiv")})
#setMethod(f="r2lBiv",signature=c("ANY","ANY"),def=
#          function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){
#              cat("Y=",class(y)," X=",class(x)," Display=",displayStyle)
#          }
#          )




################################
### Function rtlb, highest level.
### It opens a file, creates the directories, cuts the formula in two variables then calls r2lBiv

rtlb <- r2lb <- r2latexBiv <- function(formula,fileOut="bivAnalysis.tex",textBefore="",textAfter="",graphDir="graphBiv",graphName="V",type="png",displayStyle=7,limDiscreteY=10,limDiscreteX=10){
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }

    if (graphDir!="") {dir.create(graphDir,showWarnings=FALSE)}else{}

    y <- eval(formula[[2]],envir=parent.frame(n=1))
    nameY <- deparse(formula[[2]])
    x <- eval(formula[[3]],envir=parent.frame(n=1))
    if(class(x)[1]=="data.frame"){nameX <- colnames(x)}else{nameX <- deparse(formula[[3]])}

    rtlBiv(y,x,nameY,nameX,textBefore=textBefore,textAfter=textAfter,graphDir=graphDir,graphName=graphName,type=type,out="latex",displayStyle=displayStyle,limDiscreteY=limDiscreteY,limDiscreteX=limDiscreteX)
    return(invisible())
}

rthb <- r2hb <- r2htmlBiv <- function(formula,fileOut="bivAnalysis.html",textBefore="",textAfter="",graphDir="graphBiv",graphName="V",type="png",displayStyle=7,limDiscreteY=10,limDiscreteX=10){
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }
    if (graphDir!="") {dir.create(graphDir,showWarnings=FALSE)}else{}

    y <- eval(formula[[2]],envir=parent.frame(n=1))
    nameY <- deparse(formula[[2]])

    x <- eval(formula[[3]],envir=parent.frame(n=1))
    if(class(x)[1]=="data.frame"){nameX <- colnames(x)}else{nameX <- deparse(formula[[3]])}

    rtlBiv(y,x,nameY,nameX,textBefore=textBefore,textAfter=textAfter,graphDir=graphDir,graphName=graphName,type=type,out="html",displayStyle=displayStyle,limDiscreteY=limDiscreteY,limDiscreteX=limDiscreteX)
    return(invisible())
}

rtlBiv <- function(y,x,nameY,nameX,textBefore="",textAfter="",graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle=7,limDiscreteY=10,limDiscreteX=10){
    if(class(x)[1]=="data.frame"){
        cat(r2lComment("r2lb.data.frame",out=out))
        nbVar <- length(x)
        if (length(textBefore)==1) {textBefore <- rep(textBefore,time=nbVar)}else{}
        if (length(textAfter)==1) {textAfter <- rep(textAfter,time=nbVar)}else{}
        if (length(graphName)==1) {graphName <- paste(graphName,1:nbVar,sep="")}else{}
        if (length(limDiscreteX)==1) {limDiscreteX <- rep(limDiscreteX,time=nbVar)}else{}
        if (length(displayStyle)==1) {displayStyle <- rep(displayStyle,time=nbVar)}else{}

        if (class(displayStyle)=="list") {
            if(is.null(displayStyle$lim)){
                displayAux <- as.character(rep(7,time=nbVar))
            }else{
                displayAux <- as.character(rep(displayStyle$lim,time=nbVar))
            }

            displayAux[names(x)%in%displayStyle$wide] <- "wide"
            displayAux[names(x)%in%displayStyle$long] <- "long"
            displayStyle <- displayAux
        }
        for (i in 1:nbVar) {
            rtlBiv(y,x[,i],nameY,nameX[i],textBefore=textBefore[i],textAfter=textAfter[i],graphDir=graphDir,graphName=graphName[i],
                   type=type,out=out,displayStyle=displayStyle[i],limDiscreteY=limDiscreteY,limDiscreteX=limDiscreteX[i])
        }
    }else{
        y <- changeClass(y,limDiscreteY)
        x <- changeClass(x,limDiscreteX)
        if(!displayStyle%in%c("wide","long")){
            if(length(table(y)) < as.integer(displayStyle)-2 & length(table(x)) < as.integer(displayStyle)){
                displayStyle <- "wide"
            }else{
                displayStyle <- "long"
            }
        }
        nameY <- paste(unlist(strsplit(nameY,"$",fixed=TRUE)),collapse="\\$")
        nameX <- paste(unlist(strsplit(nameX,"$",fixed=TRUE)),collapse="\\$")
        tabTitle <- paste(nameY,ifelse(out=="latex"," $\\sim$ ","~"),nameX,sep="")

        cat(textBefore,"\n")
        r2lBiv(y=y,x=x,tabTitle=tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out,displayStyle=displayStyle)
        cat(textAfter,"\n")
    }
    return(invisible())
}



######################################################################
###                            Discrete                            ###
######################################################################

##############################
### r2lBiv.discrete does call r2lBivDiscrete but changes the order of the first two arguments
### So r2lBivDiscrete does dispatch according to the type of vX
