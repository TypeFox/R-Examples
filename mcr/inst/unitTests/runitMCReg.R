# Automatically generate Runit test functions from XML test cases
# 
# Author: Fabian Model
###############################################################################


cat("\n\nMCReg method comparison test cases\n\n")

library(XML)

dir.xml <- "./TestCaseCollection/"

##
## Define generic test function for one MC algorithm
##
genericMethCompTest <- function(xml.InputData,xml.Algo.Parameter,xml.Algo.RefResults,equ.prec=0.00001) {
	## InpuData
	data.x <- as.numeric(strsplit(xmlValue(getNodeSet(xml.InputData,"X")[[1]]),",",fixed=T)[[1]])
	data.y <- as.numeric(strsplit(xmlValue(getNodeSet(xml.InputData,"Y")[[1]]),",",fixed=T)[[1]])

    # If a local equivalence precision is specified use this value instead of the value of 'equ.prec'.
    # This had to be added because the approximative Passing-Bablok algorithm (PaBaLarge) differs more
    # from referencer esults than the exact algorithm (PaBa). PaBaLarge is tested against the exact
    # implementation using 34 testcases (see file "runit.PaBaLarge.R").

    if("LocalEquivalencePrecision" %in% names(xmlChildren(xml.Algo.Parameter)))
    {
        cat("\n\nUse following equ.prec value: ")
        equ.prec <- as.numeric(xmlValue(getNodeSet(xml.Algo.Parameter, "LocalEquivalencePrecision")[[1]]))
        cat(equ.prec, "\n\n")
    }
       
	## Parameter
	method.reg <- xmlValue(getNodeSet(xml.Algo.Parameter,"RegMethod")[[1]])
	method.ci <- xmlValue(getNodeSet(xml.Algo.Parameter,"CIMethod")[[1]])
   
	if(length(getNodeSet(xml.Algo.Parameter,"ErrorRatio")) > 0) 
    {
		error.ratio <- as.numeric(xmlValue(getNodeSet(xml.Algo.Parameter,"ErrorRatio")[[1]]))
	} 
    else error.ratio <- 1
    
	if(length(getNodeSet(xml.Algo.Parameter,"Bias")) > 0) 
    {
		bias.points <- as.numeric(strsplit(xmlValue(getNodeSet(xml.Algo.Parameter,"Bias")[[1]]),",",fixed=T)[[1]])
	} 
    else bias.points <- NULL
	
	## Run regression
	result <- mcreg(data.x,data.y,error.ratio=error.ratio,alpha=0.05,
					method.reg=method.reg,method.ci=method.ci,method.bootstrap.ci="Student")

	## Compare results and reference
	## Intercept
	checkEqualsNumeric(getCoefficients(result)["Intercept","EST"],
			           as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"Intercept")[[1]])),
					   msg="Check Intercept EST",tolerance=equ.prec)
	checkEqualsNumeric(getCoefficients(result)["Intercept","LCI"],
					   as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"InterceptL")[[1]])),
					   msg="Check Intercept LCI",tolerance=equ.prec)
	checkEqualsNumeric(getCoefficients(result)["Intercept","UCI"],
					   as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"InterceptU")[[1]])),
					   msg="Check Intercept UCI",tolerance=equ.prec)
	## Slope
	checkEqualsNumeric(getCoefficients(result)["Slope","EST"],
			as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"Slope")[[1]])),
			msg="Check Slope EST",tolerance=equ.prec)
	checkEqualsNumeric(getCoefficients(result)["Slope","LCI"],
			as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"SlopeL")[[1]])),
			msg="Check Slope LCI",tolerance=equ.prec)
	checkEqualsNumeric(getCoefficients(result)["Slope","UCI"],
			as.numeric(xmlValue(getNodeSet(xml.Algo.RefResults,"SlopeU")[[1]])),
			msg="Check Slope UCI",tolerance=equ.prec)
    
	## Bias at decision points
    if(length(getNodeSet(xml.Algo.Parameter,"Bias"))>0) {
	    dPoints <- as.numeric(strsplit(xmlValue(getNodeSet(xml.Algo.Parameter,"Bias")[[1]]),",",fixed=T)[[1]])		   
	    if(length(dPoints)>0) {
		    bias.result <- calcBias(result,x.levels=dPoints,alpha=0.05)
		    for(i in 1:length(dPoints)) {
			    checkEqualsNumeric(bias.result[i,"Bias"],
				    	as.numeric(strsplit(xmlValue(getNodeSet(xml.Algo.RefResults,"Bias")[[1]]),",",fixed=T)[[1]])[i],
					    msg=paste("Check Bias EST at decision point",dPoints[i]),tolerance=equ.prec)
			    checkEqualsNumeric(bias.result[i,"LCI"],
				    	as.numeric(strsplit(xmlValue(getNodeSet(xml.Algo.RefResults,"BiasL")[[1]]),",",fixed=T)[[1]])[i],
					    msg=paste("Check Bias LCI at decision point",dPoints[i]),tolerance=equ.prec)
			    checkEqualsNumeric(bias.result[i,"UCI"],
				    	as.numeric(strsplit(xmlValue(getNodeSet(xml.Algo.RefResults,"BiasU")[[1]]),",",fixed=T)[[1]])[i],
					    msg=paste("Check Bias UCI at decision point",dPoints[i]),tolerance=equ.prec)	
		    }
	    }			   
    }
}

##
## Helper function to set up closure of generic test function called with specific parameters
##
getTestFunction <- function(xml.InputData,xml.Algo.Parameter,xml.Algo.RefResults,equ.prec) {
	## Store input parameters to capture them in local environment
	## otherwise reference to global environment would be attached to test functions
	inData <- xml.InputData
	algoPara <- xml.Algo.Parameter
	algoRes <- xml.Algo.RefResults
    local.equ.prec <- equ.prec
	## Return test function that uses current local parameters
	tf <- function(){genericMethCompTest(inData,algoPara,algoRes,local.equ.prec)}
	return(tf)
}

##
## Setup dynamic test functions for all xml test files and all specified algorithms
##

testFiles <- list.files(dir.xml,pattern="^MC_TestCase.*\\.xml$",full.names=TRUE)

for (testFile in testFiles) {
	xmltc <- xmlInternalTreeParse(file=testFile)

	tc.name <- xmlValue(getNodeSet(xmltc,"/RDx_TestCase_MethodComparison/TestDescription/Name")[[1]])
	xml.InputData <- getNodeSet(xmltc,"/RDx_TestCase_MethodComparison/InputData")[[1]]
	xml.Algo <- getNodeSet(xmltc,"/RDx_TestCase_MethodComparison/TestAlgorithms/Algorithm")

    if(length(getNodeSet(xmltc,"/RDx_TestCase_MethodComparison/TestAlgorithms/DefaultEquivalencePrecision"))>0) {
        equ.prec <- as.numeric(xmlValue(getNodeSet(xmltc,"/RDx_TestCase_MethodComparison/TestAlgorithms/DefaultEquivalencePrecision")[[1]]))
    } else equ.prec <- 0.00001
    
	## For each algorithm
	for(an in 1:length(xml.Algo)) {
		## Get algo specs
		algo.name <- xmlValue(getNodeSet(xml.Algo[[an]],"Name")[[1]])
		xml.Algo.Parameter <- getNodeSet(xml.Algo[[an]],"Parameter")[[1]]
		xml.Algo.RefResults <- getNodeSet(xml.Algo[[an]],"ReferenceResults")[[1]]
		## Function name from test case and algorithm name
		fname <- paste("test.MCRegXML",tc.name,algo.name,sep=".")
		## Put test function into environment
		assign(fname,getTestFunction(xml.InputData,xml.Algo.Parameter,xml.Algo.RefResults,equ.prec))
	}
}

