estimateModel  <- function(x,graph,nresp=1)
{   
#  library(earth)

  input <- cbind(x$xpop,x$ypop)

  fmla <- as.formula(paste(paste("data.matrix(input[,",ncol(x$xpop)+1,":",ncol(x$xpop)+ncol(x$ypop),"])"),"~", paste(colnames(x$xpop),collapse= "+")))

  mars.mod <- earth(fmla,degree=2,trace=0,minspan=1,data=input,keepxy=TRUE)

  if(graph=="yes" & ncol(mars.mod$bx)>1) plotmo(mars.mod,nresponse=nresp)

  forecast <- predict(mars.mod,newdata=x$xspace)
  predicted <- data.frame(forecast)
  rownames(predicted) <- rownames(x$xspace)
  colnames(predicted) <- colnames(x$yspace)

  return(predicted)
}

