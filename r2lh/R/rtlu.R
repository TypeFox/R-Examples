################################
### Basic class
###  - logical and character are like factor
###  - integer and numeric are changed in discrete or continuous according to the number of modalities

r2lUnivFactor <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
    cat(r2lComment("rtlu.factor",out=out))
    cat(r2lUnivBeginStruct(x,tabTitle,nbColumn=2,tabSpec="|c|c|",out=out))
    if(length(levels(x))>2){
        x <- ordered(x,level=names(sort(table(x),decreasing=TRUE)))
    }else{}
    cat(r2lBuildRow(x=c(r2lBold("Frequency",out),r2lBold("Histogram",out)),hline=FALSE,out=out))
    rowTxt <- r2lBuildRow(x=c(
        r2lUnivFrequency(x, out=out),
        r2lGraphBarplot(x, graphDir=graphDir, graphName=graphName, type=type, out=out)
    ),out=out)
    cat(rowTxt)
    cat(r2lEndStruct(out=out))
    return(invisible())
}





#r2lUnivLogical <- function(x,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
#    r2lUnivFactor(x=x,textBefore=textBefore,graphDir=graphDir,graphName=graphName,type=type,out=out)
#}


#r2lUnivCharacter <- function(x,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
#    r2lUnivFactor(x=x,textBefore=textBefore,graphDir=graphDir,graphName=graphName,type=type,out=out)
#}


r2lUnivOrdered <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
    cat(r2lComment("rtlu.ordered",out=out))
    cat(r2lUnivBeginStruct(x,tabTitle,nbColumn=3,tabSpec="|c|c|c|",out=out))
    cat(r2lBuildRow(x=c(r2lBold("Frequency",out),r2lBold("Summary",out),r2lBold("Histogram",out)),hline=FALSE,out=out))
    rowTxt <- r2lBuildRow(x=c(
        r2lUnivFrequency(x, out=out),
        r2lUnivSummary(x, out=out),
        r2lGraphBarplot(x, graphDir=graphDir, graphName=graphName, type=type, out=out)
    ), out=out)
    cat(rowTxt)
    cat(r2lEndStruct(out=out))
    return(invisible())
}


r2lUnivDiscrete <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
    cat(r2lComment("rtlu.discrete",out=out))
    cat(r2lUnivBeginStruct(x,tabTitle,nbColumn=4,tabSpec="|c|c|cc|",out=out))
    cat(r2lBuildRow(x=c(r2lBold("Frequency",out),r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Histogram",out)),hline=FALSE,out=out))
    rowTxt <- r2lBuildRow(x=c(
        r2lUnivFrequency(x, out=out),
        r2lUnivSummary(x, out=out),
        r2lGraphBoxplot(x, graphDir=graphDir, graphName=graphName, type=type, out=out),
        r2lGraphBarplot(x, graphDir=graphDir, graphName=graphName, type=type, out=out)
    ), out=out)
    cat(rowTxt)
    cat(r2lEndStruct(out=out))
    return(invisible())
}


r2lUnivContinuous <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
    cat(r2lComment("rtlu.continuous",out=out))
    cat(r2lUnivBeginStruct(x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))
    cat(r2lBuildRow(x=c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Histogram",out)),hline=FALSE,out=out))
    rowTxt <- r2lBuildRow(x=c(
        r2lUnivSummary(x, out=out),
        r2lGraphBoxplot(x, graphDir=graphDir, graphName=graphName, type=type, out=out),
        r2lGraphHist(x, graphDir=graphDir, graphName=graphName, type=type, out=out)
    ), out=out)
    cat(rowTxt)
    cat(r2lEndStruct(out=out))
    return(invisible())
}


## r2lUnivIntegerNumeric <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
##     if (length(table(x))<=limDiscrete && identical(as.numeric(x),unclass(round(x)))) {
##         r2lUnivDiscrete(x=x,varName=varName,textBefore=textBefore,graphDir=graphDir,graphName=graphName,type=type,out=out)
##     }else{
##         r2lUnivContinuous(x=x,varName=varName,textBefore=textBefore,graphDir=graphDir,graphName=graphName,type=type,out=out)
##     }
## }


## r2lUnivDataFrame <- function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10) {
##     nbVar <- length(x)
##     if (length(textBefore)==1) {textBefore <- rep(textBefore,time=nbVar)}
##     if (length(graphName)==1) {graphName <- paste(graphName,1:nbVar,sep="")}
##     if (length(limDiscrete)==1) {limDiscrete <- rep(limDiscrete,time=nbVar)}
## 	cat(r2lComment("r2lu.data.frame",out=out))
##     for (i in 1:nbVar) {
##         r2lUniv(x[,i],varName=varName[i],textBefore=textBefore[i],graphDir=graphDir,graphName=graphName[i],type=type,out=out,limDiscrete=limDiscrete[i])
##     }
## }

setGeneric(name="r2lUniv",def=function(x,tabTitle,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10){standardGeneric("r2lUniv")})

setMethod(f="r2lUniv",signature="logical",def=r2lUnivFactor)
setMethod(f="r2lUniv",signature="factor",def=r2lUnivFactor)
#setMethod(f="r2lUniv",signature="character",def=r2lUnivFactor)
setMethod(f="r2lUniv",signature="ordered",def=r2lUnivOrdered)
#setMethod(f="r2lUniv",signature="integer",def=r2lUnivIntegerNumeric)
#setMethod(f="r2lUniv",signature="numeric",def=r2lUnivIntegerNumeric)
setMethod(f="r2lUniv",signature="discrete",def=r2lUnivDiscrete)
setMethod(f="r2lUniv",signature="continuous",def=r2lUnivContinuous)
#setMethod(f="r2lUniv",signature="data.frame",def=r2lUnivDataFrame)



##########################
#   Exported functions   #
##########################

r2lu <- r2latexUniv <- rtlu <- function(x,fileOut="univAnalysis.tex",textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",limDiscrete=10) {
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }
    if (graphDir!="") {
        dir.create(graphDir,showWarnings=FALSE)
    }

    if(class(x)[1]=="data.frame"){nameX <- colnames(x)}else{nameX <- deparse(substitute(x))}
    rtlUniv(x,nameX=nameX,textBefore=textBefore,textAfter=textAfter,graphDir=graphDir,graphName=graphName,type=type,out="latex",limDiscrete=limDiscrete)
    return(invisible())
}

r2hu <- r2htmlUniv <- rthu <- function(x,fileOut="",textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",limDiscrete=10) {
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }
    if (graphDir!="") {
        dir.create(graphDir,showWarnings=FALSE)
    }

    if(class(x)[1]=="data.frame"){nameX <- colnames(x)}else{nameX <- deparse(substitute(x))}
    rtlUniv(x=x,nameX=nameX,textBefore=textBefore,textAfter=textAfter,graphDir=graphDir,graphName=graphName,type=type,out="html",limDiscrete=limDiscrete)
    return(invisible())
}

rtlUniv <- function(x,nameX,textBefore="",textAfter="",graphDir="graphUniv",graphName="V",type="png",out="latex",limDiscrete=10){
    if(class(x)[1]=="data.frame"){
        cat(r2lComment("r2lu.data.frame",out=out))
        nbVar <- length(x)
        if (length(textBefore)==1) {textBefore <- rep(textBefore,time=nbVar)}else{}
        if (length(textAfter)==1) {textAfter <- rep(textAfter,time=nbVar)}else{}
        if (length(graphName)==1) {graphName <- paste(graphName,1:nbVar,sep="")}else{}
        if (length(limDiscrete)==1) {limDiscrete <- rep(limDiscrete,time=nbVar)}else{}

        for (i in 1:nbVar) {
            rtlUniv(x[,i],nameX[i],textBefore=textBefore[i],textAfter=textAfter[i],graphDir=graphDir,graphName=graphName[i],
                   type=type,out=out,limDiscrete=limDiscrete[i])
        }
    }else{
#        cat("ZZZ",x,"ZZZ",limDiscrete)
        x <- changeClass(x,limDiscrete)
 #       cat("WWWWWWWWWWW")
        if(out=="latex"){
            tabTitle <- paste(unlist(strsplit(nameX,"$",fixed=TRUE)),collapse="\\$")
        }else{
            tabTitle <- nameX
        }

        cat(textBefore,"\n")
        r2lUniv(x,tabTitle=tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
        cat(textAfter,"\n")
    }
    return(invisible())
}


r2lMainFile <- function(fileOut="main.tex",text="\\input{univ.tex}",sweave=FALSE) {
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }
    cat("
\\documentclass[a4paper, 10pt]{article}
\\usepackage{graphicx}
")
    if(sweave){cat("\\usepackage{Sweave}")}else{}
    cat("
\\title{R to LaTeX}

\\begin{document}
", text, "
\\end{document}
")
    return(invisible())
}

r2hMainFile <- function(fileOut="main.html",text="<OBJECT data='univ.html' type='text/html'></OBJECT>",sweave=FALSE) {
    if (fileOut!="") {
        on.exit(sink())
        sink(fileOut)
    }
    cat("
<HTML>
<HEAD>
<TITLE>R to HTML</TITLE>
</HEAD>
<BODY>
", text, "
</BODY>
</HTML>
")
    return(invisible())
}
