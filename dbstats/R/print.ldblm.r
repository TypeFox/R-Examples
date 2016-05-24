

    ##############################
    #### print.ldblm function ####
    ##############################

 ## description:
 ## 
 ##   print generic method. Show the most relevant attributes of a ldblm 
 ##   object in a pretty format. 
 ##     - the using method
 ##     - the kind of kernel to compute the weights
 ##     - optimal bandwidth
 ##     - Ordinary cross-validation (if method="OCV")
 ##     - Generalized cross-validation (if method="GCV")
 ##     - Aikaike information criterium (if method="AIC")
 ##     - Bayesian information criterium (if method="BIC")
 ##


print.ldblm<-function(x,...){

  # stop if the object is not a dblm object.
  if (!inherits(x, "ldblm")&&!inherits(x, "ldbglm")) 
    stop("use only with \"ldblm\" or \"ldbglm\" objects")                              
 
  #the kind of kernel to compute the weights
  kindkernel <- switch(attr(x,"kind.of.kernel"),
     "(1) Epanechnikov",
	   "(2) Biweight",
		 "(3) Triweight",
		 "(4) Normal",
		 "(5) Triangular",
		 "(6) Uniform")
		 
  # print the call
  cat("\ncall:   ")
   if (class(x)[1] == "ldblm")
    x$call[[1]]<-as.name("ldblm")
   else
    x$call[[1]]<-as.name("ldbglm")
  
  print(x$call)		
  
  # print the using method
  cat(gettextf("\nmethod.h= %s, \t kind of kernel= %s",
          attr(x,"method.h"),kindkernel),"\n")  
  
  # metric
  if(attr(x,"way")=="Z"){
   cat(gettextf("metric1: %s",attr(x,"metric1")),"\n")
   cat(gettextf("metric2: %s",attr(x,"metric2")),"\n")
  }
  
  # print the used bandwidth 
  if (attr(x,"method.h")!="user.h")
   cat(gettextf("optimal bandwidth h : %f",x$h.opt),"\n") # print h.opt  
  else 
   cat(gettextf("user bandwidth h : %f",x$h.opt),"\n") # print h.opt  
 
  # print the appropriate statistic according to the using method 
  if (!is.null(attr(x,"OCV_opt")))
   cat(gettextf("Ordinary cross-validation estimate of the prediction error : %s ",
            format(attr(x,"OCV_opt"),scientific=TRUE)),"\n")
  if (!is.null(attr(x,"GCV_opt")))
   cat(gettextf("Generalized cross-validation estimate of the prediction error : %s ",
            format(attr(x,"GCV_opt"),scientific=TRUE)),"\n")  
  if (!is.null(attr(x,"AIC_opt")))
   cat(gettextf("Akaike Information Criterion : %f ",attr(x,"AIC_opt")),"\n")
  if (!is.null(attr(x,"BIC_opt")))
   cat(gettextf("Bayesian Information Criterion : %f ",attr(x,"BIC_opt")),"\n")  
  cat("\n")
  
  return(invisible())
}