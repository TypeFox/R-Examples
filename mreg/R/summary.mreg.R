"summary.mreg" <-
function (object, digits = max(3, getOption("digits") - 3), symbolic.cor = object$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    #cat("\nResiduals:\n")
    #print(quantile(object$residuals,na.rm=TRUE))
cat("\nCoefficients:\n")
coefs<- matrix( 0, nrow=length(object$coefficients), ncol=4)
coefs[,1]<-object$coefficients
coefs[,2]<-object$se
coefs[,3]<-object$coefficients/object$se
coefs[,4]<-2*pnorm(abs( coefs[,3]), lower.tail=FALSE)
rownames(coefs)<-names(object$coefficients)
colnames(coefs)<-c("Estimate","S.E.","Z-value","Pr(>|Z|)")
if(!is.null(object$nuisance)){
  nuis <- matrix( 0, nrow=length(object$nuisance$estimate), ncol=4)
  nuis[,1] <- object$nuisance$estimate
  nuis[,2] <- object$nuisance$se
  nuis[,3] <- object$nuisance$estimate/ object$nuisance$se
  nuis[,4] <- 2*pnorm(abs( nuis[,3]), lower.tail=FALSE)
  colnames(coefs)<-c("Estimate","S.E.","Z-value","Pr(>|Z|)")
  if( object$density.name=="negbin"){
    rownames(nuis) <- c("log.disp")
  }
  if( object$density.name=="negbin.ncar"){
    rownames(nuis) <- c("log.disp","logis.int", "logis.x")
  }

  coefs <- rbind(coefs, nuis)
}

    
    
  
 
	
printCoefmat(coefs,digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
cat("\n\nDeviance:  ", object$deviance,"\n\n")


}

