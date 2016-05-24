"PermTest"<-function(obj,B=1000,...){
UseMethod("PermTest")
}

PermTest.lm<-function(obj,B=1000,...){

    # Giraudoux 24.10.2004 permutation test for linear models
    # obj = a lm object
    # B = permutation number (default = 1000)
    # values: a list including a data.frame of the p.values of each independent variable F
    # and the call
    an<-anova(obj)
    Fobs<-an[[4]]
    n<-rep(0,length(an[[4]])-1)
    
    for (i in 1:B){
        lm1<-update(obj,sample(.)~.)
        an<-anova(lm1)
        Frdm<-an[[4]]
            for (j in 1:(length(an[[4]])-1)){
                if (Frdm[j]>=Fobs[j]) n[j]<-n[j]+1
            }
        }
    n<-n/B
    names(n)<-row.names(an)[1:length(an[[4]])-1]
    n<-as.data.frame(n)
    names(n)<-"p.value"
    output<-list(resultats=n, B=B, call=match.call())
    class(output)<-c("PermTest","list")
    output
    
}


"PermTest.lme" <- function(obj, B = 1000,...){
# Giraudoux & Lancelot 18.4.2004 permutation test for linear mixed effect models,
# modified 24.10.04
# obj = a lme object
# B = permutation number (default = 1000)
# value: a data.frame of the p.values of each independent variable F
# values: a list including a data.frame of the p.values of each independent variable F
# and the call

  an <- anova(obj, type = "marginal")
  Fobs <- an[[3]]
  n <- rep(0, length(an[[3]]))
  for(i in seq(B)){
    fm <- update(obj, sample(.) ~ . )
    an <- anova(fm, type = "marginal")
    Fperm <- an[[3]]
    for(j in seq(length(an[[3]]))){
      if(Fperm[j] >= Fobs[j])
        n[j] <- n[j] + 1
      }
    }
  n <- n / B
  names(n) <- row.names(an)
  n <- as.data.frame(n)
  names(n) <- "p.value"
    output<-list(resultats=n, B=B,call=match.call())
    class(output)<-c("PermTest","list")
    output
  }


"PermTest.glm"<-function(obj,B=1000,...){

    # Giraudoux 24.10.2004 permutation test for general linear models
    # obj = a glm object of poisson or binomial family
    # B = permutation number (default = 1000)
    # value: a data.frame of the p.values of each independent variable F
    # values: a list including a data.frame of the p.values of each independent variable F
    # and the call
    if (!any(family(obj)$family=="poisson",family(obj)$family=="binomial",family(obj)$family=="quasipoisson"))
        stop("method not implemented for ",family(obj)$family," family")
       
    an<-anova(obj)
    Fobs<-an$Deviance
    n<-rep(0,length(an$Deviance)-1)
    resp<-eval(terms(formula(obj))[[2]])
    
    for (i in 1:B){
        if (is.matrix(resp)) glm1<-update(obj,permcont(.)~.)
        else glm1<-update(obj,sample(.)~.)
        an<-anova(glm1)
        Frdm<-an$Deviance
            for (j in 2:(length(an$Deviance))){
                if (Frdm[j]>=Fobs[j]) n[j-1]<-n[j-1]+1
            }
        }
    n<-n/B
    names(n)<-row.names(an)[2:length(an$Deviance)]
    n<-as.data.frame(n)
    names(n)<-"p.value"
    output<-list(resultats=n, B=B,call=match.call())
    class(output)<-c("PermTest","list")
    output
}



"print.PermTest"<-function (x, ...) {
# Giraudoux 24.10.04
# print method for PermTest objects
    if (!inherits(x, "PermTest")) 
        stop("Object must be of class 'PermTest'")
    cat("\nMonte-Carlo test\n\n")
    cat("Call: \n")
    print(x$call)
    cat("\nBased on", x$B, "replicates\n")
    cat("Simulated p-value:\n")
    print(x$resultats)
}
