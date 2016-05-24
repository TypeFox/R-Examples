

 #############################
 #### print.dblm function ####
 #############################

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


print.dblm<-function(x,...){
 
  # the call 
  cat("\ncall:   ")
  x$call[[1]]<-as.name("dblm")
  
  print(x$call)	

  # the using method

   
  if (attr(x,"method")=="eff.rank") 
    cat(gettextf("\nmethod: %s = %i",attr(x,"method"),attr(x,"ini_eff.rank")),"\n") 
  else if(attr(x,"method")=="rel.gvar")
   cat(gettextf("\nmethod: %s = %s",attr(x,"method"),round(attr(x,"ini_rel.gvar"),digits=10)),"\n") 
  else 
   cat(gettextf("\nmethod: %s",attr(x,"method")),",\t")  
   
  if (attr(x,"method")!="eff.rank"&&attr(x,"method")!="rel.gvar"){
    if (attr(x,"full.search"))
      cat("search: full \n")
    else 
      cat("search: optimize \n")
  }
 
  # metric
  if(attr(x,"way")=="Z")
    cat(gettextf("metric: %s",attr(x,"metric")),"\n")
 
  # effective rank used
  if (attr(x,"method")=="eff.rank"||attr(x,"method")=="rel.gvar")
    cat(gettextf("Used effective rank = %i",x$eff.rank),"\n")
  else    
    cat(gettextf("Optimal effective rank = %i",x$eff.rank),"\n")
  
  # Relative geometric variability
  cat(gettextf("Relative geometric variability = %f",x$rel.gvar),"\n")
  
  # print the appropriate statistic according to the using method 
  if (attr(x,"method")=="OCV")
   cat(gettextf("Ordinary cross-validation estimate of the prediction error : %s ",
          format(x$ocv,scientific=TRUE)),"\n")
  if (attr(x,"method")=="GCV")
   cat(gettextf("Generalized cross-validation estimate of the prediction error : %s ",
          format(x$gcv,scientific=TRUE)),"\n")  
  if (attr(x,"method")=="AIC")
   cat(gettextf("Akaike Information Criterion : %f ",x$aic),"\n")
  if (attr(x,"method")=="BIC")
   cat(gettextf("Bayesian Information Criterion : %f ",x$bic),"\n")
   
 cat("\n")
 return (invisible())
}