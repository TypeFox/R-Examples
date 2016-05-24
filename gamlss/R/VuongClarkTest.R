# Created from the function of Giampierro Marra
VC.test <- function (obj1, obj2, sig.lev = 0.05) 
{
  # check for gamlss models 
  if (!is.gamlss(obj1))  stop(paste("Object 1 is not a gamlss object", "\n", ""))
  if (!is.gamlss(obj2))  stop(paste("Object 2 is not a gamlss object", "\n", ""))
  # get the likelihood function
  l1 <- logLik(obj1)[1]
  l2 <- logLik(obj2)[1]
  if (l1 == l2) 
    stop("The two competing models have identical log-likelihoods!")
  # get degrees of freedom for the fit?
  p1 <- obj1$df.fit 
  p2 <- obj2$df.fit
  # no of obsernarions
  n <- obj1$N
 # I need to calculate the likelihood increments here  
 # us preict to get the fitted parameters
    pobj1 <- predictAll(obj1)
    pobj2 <- predictAll(obj2)
 # the names of te families
   fname1 <- obj1$family[1]
   fname2 <- obj2$family[1]
   nopar1 <- eval(gamlss.family(fname1))$nopar
   nopar2 <- eval(gamlss.family(fname2))$nopar
#  param1 <- eval(gamlss.family(fname1))$parameter
#  param2 <- eval(gamlss.family(fname2))$parameter
 pdfName1 <- paste("d", fname1, sep="") 
 pdfName2 <- paste("d", fname2, sep="") 
binomial1 <- fname1%in%.gamlss.bi.list
binomial2 <- fname2%in%.gamlss.bi.list
  model1 <- deparse(substitute(obj1))
  model2 <- deparse(substitute(obj2))  
if (binomial1)
{
  linc1 <- switch( nopar1, 
                   eval(call( pdfName1, pobj1$y, bd=obj1$bd, ,mu=pobj1$mu, log=TRUE)),
                   eval(call( pdfName1, pobj1$y, bd=obj1$bd, mu=pobj1$mu, sigma=pobj1$sigma, log=TRUE )),
                   eval(call( pdfName1, pobj1$y, bd=obj1$bd,mu=pobj1$mu, sigma=pobj1$sigma, nu=pobj1$nu, log=TRUE )),
                   eval(call( pdfName1, pobj1$y, bd=obj1$bd,mu=pobj1$mu, sigma=pobj1$sigma,  nu=pobj1$nu, tau=pobj1$tau,log=TRUE))
  ) 
} else
{  
linc1 <- switch( nopar1, 
  eval(call( pdfName1, pobj1$y, mu=pobj1$mu, log=TRUE)),
  eval(call( pdfName1, pobj1$y, mu=pobj1$mu, sigma=pobj1$sigma, log=TRUE )),
  eval(call( pdfName1, pobj1$y, mu=pobj1$mu, sigma=pobj1$sigma, nu=pobj1$nu, log=TRUE )),
  eval(call( pdfName1, pobj1$y, mu=pobj1$mu, sigma=pobj1$sigma, nu=pobj1$nu, tau=pobj1$tau, log=TRUE))
        )
}
if (binomial2)
{
  linc2 <- switch( nopar1, 
                   eval(call( pdfName2, pobj2$y, bd=obj2$bd, ,mu=pobj2$mu, log=TRUE)),
                   eval(call( pdfName2, pobj2$y, bd=obj2$bd, mu=pobj2$mu, sigma=pobj2$sigma, log=TRUE )),
                   eval(call( pdfName2, pobj2$y, bd=obj2$bd,mu=pobj2$mu, sigma=pobj2$sigma, nu=pobj2$nu, log=TRUE )),
                   eval(call( pdfName2, pobj2$y, bd=obj2$bd,mu=pobj2$mu, sigma=pobj2$sigma,  nu=pobj2$nu, tau=pobj2$tau,log=TRUE))
  )   
} else
{
linc2 <- switch( nopar2, 
                 eval(call( pdfName2, pobj2$y, mu=pobj2$mu, log=TRUE)),
                 eval(call( pdfName2, pobj2$y, mu=pobj2$mu, sigma=pobj2$sigma, log=TRUE )),
                 eval(call( pdfName2, pobj2$y, mu=pobj2$mu, sigma=pobj2$sigma, nu=pobj2$nu, log=TRUE )),
                 eval(call( pdfName2, pobj2$y, mu=pobj2$mu, sigma=pobj2$sigma, nu=pobj2$nu, tau=pobj2$tau, log=TRUE))
)
}
   li12 <- linc1 - linc2 # deviance increments
      w <- sqrt(var(li12) * (n - 1))
     vt <- (l1 - l2 - (p1 - p2)/2 * log(n))/w # the Vuong est
  li12b <- li12 - (p1 - p2)/(2 * n) * log(n)
      b <- sum(li12b > 0)                     # the Clarke test
# conditiong on the P values Vuong test
  if (abs(vt) <= qnorm(1 - sig.lev/2)) # case 1
cat(" Vuong's test:", round(vt,3), "it is not possible to discriminate between models:", model1, "and", model2, "\n")
  if (vt > qnorm(1 - sig.lev/2))       # case 2
cat(" Vuong's test:",  round(vt,3), "model", model1, "is preferred over", model2, "\n")    
  if (vt < -qnorm(1 - sig.lev/2))      # case 3
    cat(" Vuong's test:",  round(vt,3), "model", model2, "is preferred over", model1, "\n")     
#    dvt <- paste(paste(paste(" Vuong's test:", model2),  "is preferred over"), model1)
#"Clarke's test: it is not possible to discriminate between the two models."
 db <- paste(paste(paste("Clarke's test: it is not possible to discriminate between models:", model1),"and"), model2) 
  if (b >= n/2) 
    {
    pvalue <- 2 * (1 - pbinom(b - 1, n, 0.5))
    if (pvalue <= sig.lev) 
      cat("Clarke's test:",b,"p-value=",  round(pvalue,4), model1, "is preferred over", model2, "\n")
     # db <- paste(paste(paste("Clarke's test:", model1),  "is preferred over"), model2)
    else 
     cat("Clarke's test:",b,"p-value=", round(pvalue,4), "it is not possible to discriminate between models:", model1,"and", model2, "\n")  
   }
  if (b < n/2) 
    {
    pvalue <- 2 * (pbinom(b, n, 0.5))
    if (pvalue <= sig.lev) 
    cat("Clarke's test:",b,"p-value=", round(pvalue,4), model2, "is preferred over", model1, "\n")
    else 
      cat("Clarke's test:",b,"p-value=", round(pvalue,4), "it is not possible to discriminate between models:", model1,"and", model2, "\n")    
    #  db <-  paste(paste(paste("Clarke's test:", model2),  "is preferred over"), model1)
      #"Clarke's test: Model 2 is preferred over Model 1."
   }
  #cat("\n", dvt, "\n", sep = "")
  #cat(db, "\n\n", sep = "")
}
