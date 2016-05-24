cat("\n\nmcDeming.R method comparison test cases\n\n")

test.mcDeming.call <- function() 
{
    
    data(creatinine)
    crea <- creatinine[complete.cases(creatinine),]
    
    checkException( mcr:::mc.deming(crea[,1], crea[,2], error.ratio=numeric(0)))
    checkException( mcr:::mc.deming(crea[,1], crea[,2], error.ratio=0))
    checkException( mcr:::mc.deming(crea[,1], crea[,2], error.ratio=as.numeric(NA)))
    checkException( mcr:::mc.deming(crea[,1], crea[,2], error.ratio="1"))
    
    ## checking for correct results
   
    res <- list(b0=-0.05891341, b1=1.054539, se.b0=0.04604315, se.b1=0.03534361, xw=1.221111, weight=rep(1,length(crea[,1])))
    
    checkEquals( mcr:::mc.deming(crea[,1], crea[,2], error.ratio=1), res, tolerance=10e-7 )
}