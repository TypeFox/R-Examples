`m4plPersonParameters` <-
function(x,s=1/1.702,b=0,c=0,d=1,m=0,model="T",prior="uniform", more=FALSE){
 if (model != "T") {
  if (more == FALSE) res <- data.frame(t(apply(x,1,m4plEstimate,     b=b, s=s, c=c, d=d, m=m, model=model, prior=prior)))
  if (more == TRUE)  res <- apply(x,1,             m4plEstimateMore, b=b, s=s, c=c, d=d, m=m, model=model, prior=prior)
  }
 if (model == "T") {
  if (more == FALSE) {
   res <- t(data.frame(t(apply(x,1,m4plEstimate, b=b, s=s, c=c, d=d, m=m, model=model, prior=prior))))
   rownames(res) <- (1:length(res)); colnames(res) <- "T"
   res <- data.frame(res)
   }
  if (more == TRUE)  res <- apply(x,1,m4plEstimateMore, b=b, s=s, c=c, d=d, m=m, model=model,    prior=prior)
  }
 # Class definitions to be implemented later so that summary functions
 # will be replaced by print.m4pl() and summary.m4pl()
 # if (more == FALSE) class(res) <- "m4pl"
 # if (more == TRUE)  class(res) <- "m4plMore"
 if (more == TRUE) class(res) <- c("m4plMore", "list")
 return(res)
 }

 #model <- "C"
 #test  <- m4plPersonParameters(x=X, b=b, s=1/a, c=c, d=d, m=0, model=model,prior="uniform", more=TRUE)

#`is.m4plMore` <-
#function(x) {
# inherits(x, "m4plMore")
# }
