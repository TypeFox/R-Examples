# '
A_clere<-function(y=y,x=x,g=NULL,analysis="aic"){
   requireNamespace("clere")
   x=1*as.matrix(x)
   if(is.null(g)){
      g=min(5,ncol(x))#ncol(x)
   }
   model <- clere::fitClere(y = as.numeric(y), x = x, g = g, plotit = FALSE, analysis=analysis)
   A=c(model@intercept,rowMeans(model@Bw))
   return(A)
}