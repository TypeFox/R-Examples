cat("\n\nmcWDeming.R method comparison test cases\n\n")

## check implementation of weighted deming-regression

test.mc.analytical.mc.call <- function() 
{
    data(creatinine)
    crea <- creatinine[complete.cases(creatinine),]
    
    res <- mcr:::mc.wdemingConstCV(crea[,1], crea[,2], error.ratio=1)
    
    checkEquals( res$b1, 1.111956340733 )
    checkEquals( res$b0, -0.1254944949597)
}