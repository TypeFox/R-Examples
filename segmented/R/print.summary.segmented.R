`print.summary.segmented` <-
function(x, short = x$short, var.diff = x$var.diff, 
    digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...){
    cat("\n\t***Regression Model with Segmented Relationship(s)***\n\n")
    cat( "Call: \n" )
    print( x$call )
    cat("\nEstimated Break-Point(s):\n ")
    print(round(x$psi[,-1],3)) #era "signif(,4)"
#    cat("\nt value for the gap-variable(s) V: ",x$gap[,3],"\n")
#if(any(abs(x$gap[,3])>1.96)) cat("    Warning:", sum(abs(x$gap[,3])>1.96),"gap coefficient(s) significant at 0.05 level\n")
    if(short){ 
    cat("\nDifference-in-slopes parameter(s):\n")
    #print(x$Ttable[(nrow(x$Ttable)-nrow(x$psi)+1):nrow(x$Ttable),])}
    nome<-rownames(x$psi)
    #nome<-as.character(parse("",text=nome))
    #aa<-grep("U",rownames(x$Ttable))
    #bb<-unlist(sapply(nome,function(xx){grep(xx,rownames(x$Ttable))},simplify=FALSE,USE.NAMES=FALSE))
    #cc<-intersect(aa,bb) #indices of diff-slope parameters
    nomiU<-rownames(x$gap)
    #idU<-match(nomiU,rownames(x$Ttable))
    print(x$Ttable[nomiU,])
      } else {cat("\nMeaningful coefficients of the linear terms:\n")
        if(is.null(dim(x$Ttable))){
        print(x$Ttable)
        #printCoefmat(matrix(x$Ttable,nrow=1,ncol=4,dimnames=list(" ",names(x$Ttable))),has.Pvalue=FALSE)
        } else {
        printCoefmat(x$Ttable, digits = digits, signif.stars = signif.stars,na.print = "NA", ...)
        }
        
        }
if("summary.lm"%in%class(x)){ #for lm
    if(var.diff){
    for(i in 1:length(x$sigma.new)){
    cat("\nResidual standard error ",i,":", format(signif(x$sigma.new[i], 
        digits)), "on", x$df.new[i], "degrees of freedom")}
    cat("\n")    
    } else {
    cat("\nResidual standard error:", format(signif(x$sigma, 
        digits)), "on", x$df[2], "degrees of freedom\n")}
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
        cat(",  Adjusted R-squared:", formatC(x$adj.r.squared, 
            digits = digits), "\n")}
        }
if("summary.glm"%in%class(x)){ #for glm
    cat("(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format.default(c("Null", 
            "Residual"), width = 8, flag = ""), "deviance:"), 
            format(unlist(x[c("null.deviance", "deviance")]), 
                digits = max(5, digits + 1)), " on", format(unlist(x[c("df.null", 
                "df.residual")])), " degrees of freedom\n"), 
            1, paste, collapse = " "), "AIC: ", format(x$aic, 
            digits = max(4, digits + 1)), "\n", sep = "")
        }
if("summary.Arima"%in%class(x)){#for Arima 
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS") 
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
            ",  log likelihood = ", format(round(x$loglik, 2)), 
            ",  aic = ", format(round(x$aic, 2)), "\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits), 
        ",  part log likelihood = ", format(round(x$loglik, 2)), 
        "\n", sep = "")
    }
invisible(x) 
cat("\nConvergence attained in",x$it,"iterations with relative change",x$epsilon,"\n")
}

