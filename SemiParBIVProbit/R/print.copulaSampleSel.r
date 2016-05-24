print.copulaSampleSel <- function(x, ...){

  if(x$BivD=="FGM")  {cop <- "FGM"                ;lind <- "atanh"} 
  if(x$BivD=="AMH")  {cop <- "AMH"                ;lind <- "atanh"} 
  if(x$BivD=="N")    {cop <- "Gaussian"           ;lind <- "atanh"}
  if(x$BivD=="F")    {cop <- "Frank"              ;lind <- "identity"}       
  if(x$BivD=="C0")   {cop <- "Clayton"            ;lind <- "log"}   
  if(x$BivD=="C90")  {cop <- "90\u00B0 Clayton"   ;lind <- "log(- \u00B7)"}                 
  if(x$BivD=="C180") {cop <- "180\u00B0 Clayton"  ;lind <- "log"}                    
  if(x$BivD=="C270") {cop <- "270\u00B0 Clayton"  ;lind <- "log(- \u00B7)"}    
  if(x$BivD=="J0")   {cop <- "Joe"                ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="J90")  {cop <- "90\u00B0 Joe"       ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="J180") {cop <- "180\u00B0 Joe"      ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="J270") {cop <- "270\u00B0 Joe"      ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="G0")   {cop <- "Gumbel"             ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="G90")  {cop <- "90\u00B0 Gumbel"    ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="G180") {cop <- "180\u00B0 Gumbel"   ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="G270") {cop <- "270\u00B0 Gumbel"   ;lind <- "log(- \u00B7 - 1)"}   
  
  
     bin.link <- x$bl
  
     as.p <- x$theta.a
  
     main.t <- "\nCOPULA:  "     

    if(x$margins[1]=="probit")                                  m1l <- "probit"
    if(x$margins[1]=="logit")                                   m1l <- "logit"
    if(x$margins[1]=="cloglog")                                 m1l <- "cloglog"
    if(x$margins[1]=="cauchit")                                 m1l <- "cauchit"  

 
    if(x$margins[2] %in% c("N","GU","rGU","LO","GAi") )               m2l <- "identity"
    if(x$margins[2] %in% c("LN","WEI","iG","GA","DAGUM","SM","FISK") ) m2l <- "log"   
    if(x$margins[2] %in% c("BE") )                              m2l <- "qlogis"   
    
      cat(main.t,cop) 
    
      if(x$margins[1] %in% bin.link) cat("\nMARGIN 1: Bernoulli")  
      
      if(x$margins[2]=="N")      cat("\nMARGIN 2: Gaussian")  
      if(x$margins[2]=="GU")     cat("\nMARGIN 2: Gumbel")    
      if(x$margins[2]=="rGU")    cat("\nMARGIN 2: reverse Gumbel")  
      if(x$margins[2]=="LO")     cat("\nMARGIN 2: logistic")   
      if(x$margins[2]=="LN")     cat("\nMARGIN 2: log-normal") 
      if(x$margins[2]=="WEI")    cat("\nMARGIN 2: Weibull") 
      if(x$margins[2]=="iG")     cat("\nMARGIN 2: inverse Gaussian") 
      if(x$margins[2]%in%c("GA","GAi"))     cat("\nMARGIN 2: gamma")   
      if(x$margins[2]=="BE")     cat("\nMARGIN 2: beta")    
      if(x$margins[2]=="DAGUM")  cat("\nMARGIN 2: Dagum")
      if(x$margins[2]=="SM")     cat("\nMARGIN 2: Singh-Maddala") 
      if(x$margins[2]=="FISK")     cat("\nMARGIN 2: Fisk") 
      
    
    
    
  

  cat("\n\nEQUATION 1")
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])

  cat("\nEQUATION 2")
  
  
  
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  
  
  
  if(!is.null(x$X3) && is.null(x$X4) ){
  
  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  
  }
  
  
  
  
  
  
  if(!is.null(x$X3) && !is.null(x$X4) &&  is.null(x$X5)){
  

  cat("\nEQUATION 3")
  if(x$margins[2] != "BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n") 
  cat("Formula: ");  print(x$formula[[3]]) 
  

  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: ");  print(x$formula[[4]]) 
  
  } 
  
  
  
  
  
  
  if(!is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)){
  
  cat("\nEQUATION 3")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  
  
  cat("\nEQUATION 4")
  if(x$margins[2] %in% c("DAGUM","SM")) cat("\nLink function for nu:","log","\n")  
  cat("Formula: "); print(x$formula[[4]]) 
    

  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
  
  }  
  
  
  
  
  cont2par <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")  
  cont3par <- c("DAGUM","SM")  
  
  
  cat("\n")
            
  if(x$margins[2] %in% cont2par ) cat("n = ",x$n,"  n.sel = ", x$n.sel,"\nsigma2 = ",x$sigma2.a, "\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$margins[2] %in% cont3par ) cat("n = ",x$n,"  n.sel = ", x$n.sel,"\nsigma2 = ",x$sigma2.a, "  nu = ",x$nu.a, "\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

invisible(x)

}

