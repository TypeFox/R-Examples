cat("\n\nMCResultMethods.R, MCResultJackknifeMethods.R, MCResultAnalyticalMethods.R,\nMCResultBCaMethods.R and MCResultResamplingMethods.R method comparison test cases\n\n")

## Check if method 'getCoefficients'of class 'MCResult' returns correct values computed for dummy-data.
## Assumes a (2 x 4) matrix with rownames "Intercept" and "Slope" and colnames "EST", "SE", "LCI", "UCI" for the test data

test.getCoefficients.call <- function() 
{    
	## Some test data
	data.x <- 1:10
	data.y <- c(1,3,2,4,6,5,7,9,8,10)
	
	obj.lr.num.jk <- mcreg(data.x, data.y, method.reg="LinReg", method.ci="jackknife")
    obj.lr.num.an <- mcreg(data.x, data.y, method.reg="LinReg", method.ci="analytical")

    
    coef.mat.jk <- matrix(  c(0.20000000, 0.96363636, 0.55758112, 0.08763732, -1.08578438, 0.76154434, 1.48578438, 1.16572839),
                            ncol=4, dimnames=list(c("Intercept", "Slope"), c("EST", "SE", "LCI", "UCI")))
    coef.mat.an <- matrix(  c(0.2, 0.9636364, 0.5862051, 0.0944755, -1.1517913, 0.7457755, 1.5517913, 1.1814973),
                            ncol=4, dimnames=list(c("Intercept", "Slope"), c("EST", "SE", "LCI", "UCI")))
	getCoefficients(obj.lr.num.jk)
	checkEquals(getCoefficients(obj.lr.num.jk), coef.mat.jk, tolerance=10e-7)
    checkEquals(getCoefficients(obj.lr.num.an), coef.mat.an, tolerance=10e-7)
}

## Check if method 'getData' of class 'MCResult' returns the correct data use for 
## Assumes a (10 x 2) matrix with rownames S1, ..., S10 and colnames "X" and "Y" for the test data.

test.getData.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    mat.xy <- matrix(c(data.x, data.y), ncol=2, dimnames=list(paste("S", 1:10, sep=""), c("X", "Y")))
 
    obj.lr.num <- mcreg(data.x, data.y, method.reg="LinReg")

    checkEquals(getData(obj.lr.num), mat.xy)    
}

## Check if method 'getErrorRatio' of class 'MCResult' returns the correct error ratio 
## Assumes error ratio=0.8 for deming regression


test.getErrorRatio.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    lmbd <- 0.8
    
    obj.dr.num <- mcreg(data.x, data.y, method.reg="Deming", error.ratio=0.8)
    obj.lr.num <- mcreg(data.x, data.y, method.reg="LinReg")
    

    checkEquals(getErrorRatio(obj.dr.num), lmbd)
    checkEquals(getErrorRatio(obj.lr.num), 1)    
}

test.getWeights.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    w <- 1/data.x^2
    w1 <- rep(1,length(data.x))
    names(w) <-  paste("S", 1:10, sep="")
    names(w1) <- paste("S", 1:10, sep="")
    
    obj.wlr.num <- mcreg(data.x, data.y, method.reg="WLinReg")
    obj.lr.num <- mcreg(data.x, data.y, method.reg="LinReg")
    

    checkEquals(getWeights(obj.wlr.num), w)
    checkEquals(getWeights(obj.lr.num), w1)    
}



## Check if method 'getResiduals' of class 'MCResult' yields the correct result for dummy-data.

test.getResiduals.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10)
       
    obj.lr.num <- mcreg(data.x, data.y, method.reg="Deming", error.ratio=0.8)

    ## currently no reference data available
}



## Check if method 'getFitted' of class 'MCResult' yields consistent results.

test.getFitted.call <- function()
{

    data(creatinine)
    crea <- creatinine[complete.cases(creatinine),]
    
    obj.lr.num <- mcreg(crea[,1],crea[,2], method.reg="Deming", error.ratio=0.8)
    
    d_hat1 <- crea-getResiduals(obj.lr.num)[,c("x","y")]
    d_hat2 <- getFitted(obj.lr.num)
    
    checkEquals(d_hat1[,1],d_hat2[,1],tolerance = 10^-7)
    checkEquals(d_hat1[,2],d_hat2[,2],tolerance = 10^-7)
}


## Check if method 'getRegmethod' of class 'MCResult' yields the correct names of the regression.

test.getRegmethod.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10)
    
    obj.LR <- mcreg(data.x, data.y, method.reg = "LinReg")
    obj.WLR <- mcreg(data.x, data.y, method.reg = "WLinReg")
    obj.PaBa <- mcreg(data.x, data.y, method.reg = "PaBa")
    obj.Dem <- mcreg(data.x, data.y, method.reg = "Deming")
    obj.WDem <- mcreg(data.x, data.y, method.reg = "WDeming")
    
    checkEquals( getRegmethod(obj.LR), "LinReg")
    checkEquals( getRegmethod(obj.WLR), "WLinReg")
    checkEquals( getRegmethod(obj.PaBa), "PaBa")
    checkEquals( getRegmethod(obj.Dem), "Deming")
    checkEquals( getRegmethod(obj.WDem), "WDeming")
}

## Check if method 'calcCUSUM' of class 'MCResult' yields the correct result for dummy-data.
## Assumes a list with specific elements.

test.calcCUSUM.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    
    obj.jk <- mcreg(data.x, data.y, method.reg = "LinReg")
    obj.res <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap")
    obj.bca <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap", method.bootstrap.ci="BCa")
    obj.ana <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="analytical")
    
    list.obj <- list(nPos=5, nNeg=5, cusum=c(-1, 0,-1, -2, -1, -2, -1, 0, -1, 0), max.cusum=2)
    #names(list.obj$cusum) <- paste("S", 1:10, sep="")
    
    checkEquals(calcCUSUM(obj.jk), list.obj)
}

## Check if method 'calcBias' of class 'MCResult' yields the correct result for dummy-data.
## Assumes a matrix, see below.

test.calcBias.call <- function()
{
    ## Some test data
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    obj.jk <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="jackknife")
    obj.res <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap")
    obj.bca <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap", method.bootstrap.ci="BCa")
    obj.ana <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="analytical")
    
    mat.abs <- matrix(  c(1, 2, 3, 0.16363636, 0.12727273, 0.09090909, 0.48383232, 0.41551331, 0.35576609, -0.95208296, -0.83090269, -0.72948899, 1.27935569, 1.08544814, 0.91130718),
                    ncol=5, dimnames=list(c("X1", "X2", "X3"), c("Level", "Bias", "SE", "LCI", "UCI")))

    mat.rel1 <- matrix(  c(1, 2, 3, 0.16363636, 0.12727273, 0.09090909, 0.48383232, 0.41551331, 0.35576609, -0.95208296, -0.83090269, -0.72948899, 1.27935569, 1.08544814, 0.91130718),
                    ncol=5, dimnames=list(c("X1", "X2", "X3"), c("Level", "Prop.bias(%)", "SE", "LCI", "UCI")))
    mat.rel1[,c("Prop.bias(%)", "SE", "LCI", "UCI")] <- 100*mat.rel1[,c("Prop.bias(%)", "SE", "LCI", "UCI")]/1:3               

     mat.rel2 <- matrix(  c(1, 2, 3, 0.16363636, 0.12727273, 0.09090909, 0.48383232, 0.41551331, 0.35576609, -0.95208296, -0.83090269, -0.72948899, 1.27935569, 1.08544814, 0.91130718),
                    ncol=5, dimnames=list(c("X1", "X2", "X3"), c("Level", "Prop.bias", "SE", "LCI", "UCI")))
   
     mat.rel2[,c("Prop.bias", "SE", "LCI", "UCI")] <- mat.rel2[,c("Prop.bias", "SE", "LCI", "UCI")]/1:3               

            
    checkEquals( calcBias(obj.jk, x.levels=1:3), mat.abs)
    checkEquals( calcBias(obj.jk, x.levels=1:3, type="proportional"), mat.rel1)
    checkEquals( calcBias(obj.jk, x.levels=1:3, type="proportional", percent=FALSE), mat.rel2)
    
    
    checkException(calcBias(obj.jk, x.levels=1:3, type="bla"))
    checkException(calcBias(obj.jk, x.levels=1:3, percent=as.numeric(NA)))
    checkException(calcBias(obj.jk, x.levels=1:3, percent=numeric(0)))
    checkException(calcBias(obj.jk, x.levels=1:3, percent="non-numeric"))
    
    
    checkException(calcBias(obj.jk, x.levels=as.numeric(NA)))
    checkException(calcBias(obj.jk, x.levels=numeric(0)))
    checkException(calcBias(obj.jk, x.levels="non-numeric"))
    
    checkException(calcBias(obj.jk, alpha=as.numeric(NA)))
    checkException(calcBias(obj.jk, alpha=numeric(0)))
    checkException(calcBias(obj.jk, alpha="non-numeric"))
    checkException(calcBias(obj.jk, alpha=-1))
    checkException(calcBias(obj.jk, alpha=2))
}

## Check whether incorrect specification of parameters throws an error when
## method 'calcResponse' is used for objects with classes having common superclass
# 'MCResult'.

test.calcResponse.call <- function()
{
    data.x <- 1:10
    data.y <- c(1,3,2,4,6,5,7,9,8,10) 
    obj.jk <- mcreg(data.x, data.y, method.reg = "LinReg")
    obj.res <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap")
    obj.bca <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="bootstrap", method.bootstrap.ci="BCa")
    obj.ana <- mcreg(data.x, data.y, method.reg = "LinReg", method.ci="analytical")
    
    checkException( calcResponse(obj.jk, x.levels=as.numeric(NA)) )
    checkException( calcResponse(obj.res, x.levels=as.numeric(NA)) )
    checkException( calcResponse(obj.bca, x.levels=as.numeric(NA)) )
    checkException( calcResponse(obj.ana, x.levels=as.numeric(NA)) )
    
    checkException( calcResponse(obj.jk, x.levels=numeric(0)) )
    checkException( calcResponse(obj.res, x.levels=numeric(0)) )
    checkException( calcResponse(obj.bca, x.levels=numeric(0)) )
    checkException( calcResponse(obj.ana, x.levels=numeric(0)) )
    
    checkException( calcResponse(obj.jk, x.levels="non-numeric")  )
    checkException( calcResponse(obj.res, x.levels="non-numeric") )
    checkException( calcResponse(obj.bca, x.levels="non-numeric") )
    checkException( calcResponse(obj.ana, x.levels="non-numeric") )
    
    checkException( calcResponse(obj.jk, alpha="non-numeric")  )
    checkException( calcResponse(obj.res, alpha="non-numeric") )
    checkException( calcResponse(obj.bca, alpha="non-numeric") )
    checkException( calcResponse(obj.ana, alphas="non-numeric") )
    
    checkException( calcResponse(obj.jk, alpha=2)  )
    checkException( calcResponse(obj.res, alpha=2) )
    checkException( calcResponse(obj.bca, alpha=2) )
    checkException( calcResponse(obj.ana, alpha=2) )
    
    checkException( calcResponse(obj.jk, alpha=-2)  )
    checkException( calcResponse(obj.res, alpha=-2) )
    checkException( calcResponse(obj.bca, alpha=-2) )
    checkException( calcResponse(obj.ana, alpha=-2) )
    
    checkException( calcResponse(obj.jk, alpha=as.numeric(NA)) )
    checkException( calcResponse(obj.res, alpha=as.numeric(NA)) )
    checkException( calcResponse(obj.bca, alpha=as.numeric(NA)) )
    checkException( calcResponse(obj.ana, alpha=as.numeric(NA)) )
}
