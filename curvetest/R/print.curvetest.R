print.curvetest <-
structure(function(x,...){
    tests<-x
    ss<-function(x, d=3) signif(x,digits=d)
    mycat<-function(...) cat(paste(..., sep=""))
    mycat("\n \t \t ======================\n")
    mycat("\t\t Curve Test  Procedures \n")
    mycat("\t \t ====================== ")
    mycat("\n The p-value to test H0:", if (!"equal.var"%in%names(tests)) 
        "f(x)=0 is \t"
    else "f1(x)=f2(x) is ", ss(tests$p), ".\n")
    mycat("\n With test statistics equals\t\t", ss(tests$Statistic), ",", 
     "\n Estimated degree of freedom is \t", 
        ss(tests$eDF), ".\n")
    if ("equal.var"%in%names(tests)){
     if(tests$equal.var) {
        mycat("\n Equal variances assumed. \n")
        mycat(" Estimated common sigma^2 is \t \t", ss(tests$sigma.square), 
            ".\n")
        mycat(" ====================\n")
     } else if (!tests$equal.var) {
        mycat("\n Unequal variances assumed.\n")
        mycat(" Estimated sigma^2 are\t \t \t", ss(tests$esigma1^2), "  \n")
        mycat(" and \t\t\t\t \t", ss(tests$esigma2^2), ".\n")
        mycat(" =========================\n")
    }
    }
}, modifiers = "public")
