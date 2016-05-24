anova.CNVassoc<-function(object,...){

  if (length(list(object, ...)) <= 1L) 
     stop("anova requires two nested models")
  objects<-list(object,...)
 
  responses <- as.character(lapply(objects, function(x) deparse(x$terms[[2L]])))

  sameresp<-responses==responses[1L]

  if (!all(sameresp)) {
     stop("models might have the same response ")
    }

  ns <- sapply(objects, function(x) length(x$y))
  if (any(ns != ns[1L])) 
     stop("models were not all fitted to the same size of dataset")

  mod1<-objects[[1]]
  mod2<-objects[[2]]
  chi<-2*(logLik(mod1)[1]-logLik(mod2)[1])
  chi<-as.double(abs(chi))
  df<-logLik(mod1)[2]-logLik(mod2)[2]
  df<-as.integer(abs(df))
  pvalue<-pchisq(abs(chi),abs(df),lower.tail=FALSE)

  out<-c("chi"=chi,"df"=df,"pvalue"=pvalue)
  
  out<-list(chi.test=out,mod1=mod1,mod2=mod2)
  
  class(out)<-"anova.CNVassoc"
  
  out

}
  
