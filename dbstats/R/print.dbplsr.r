

 ##############################
 #### print.dbplsr function ####
 ##############################

 ## description:
 ##
 ##   print generic method. Show the most relevant attributes of a dblm
 ##   object in a pretty format.
 ##     - the call
 ##     - the using method
 ##     - metric
 ##     - Ordinary cross-validation (if method="OCV")
 ##     - Generalized cross-validation (if method="GCV")
 ##     - Aikaike information criterium (if method="AIC")
 ##     - Bayesian information criterium (if method="BIC")


print.dbplsr<-function(x,...){

  # the call
  cat("\ncall:   ")
  x$call[[1]]<-as.name("dbplsr")
  print(x$call)

  # the using method
    cat(gettextf("\nnumber of components: %i",x$ncomp),"\n")

  # metric
  if(attr(x,"way")=="Z")
    cat(gettextf("metric: %s",attr(x,"metric")),"\n")
  
  if(x$method!="ncomp"){
    cat(gettextf("\noptimal number of components using method %s: ",x$method," "))
    cat(x$ncomp.opt)
  }
    
  # print the appropriate statistic according to the using method
   switch(x$method,
	  "OCV"= cat(gettextf("\noptimal Ordinary cross-validation : %s ",
              format(x$ocv[x$ncomp.opt],scientific=TRUE)),"\n"),
    "GCV"=cat(gettextf("\noptimal Generalized cross-validation : %s ",
              format(x$gcv[x$ncomp.opt],scientific=TRUE)),"\n"),
		"AIC"=cat(gettextf("\noptimal Akaike Information Criterion : %s ",
              format(x$aic[x$ncomp.opt],scientific=TRUE)),"\n"),
		"BIC"= cat(gettextf("\noptimal Bayesian Information Criterion : %s ",
              format(x$bic[x$ncomp.opt],scientific=TRUE)),"\n")
  )

 cat("\n")
 return (invisible())
}