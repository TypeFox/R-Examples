print.summary.copulaReg <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

   cont2par <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")  
   cont3par <- c("DAGUM","SM") 
 
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
     
     
   main.t <- "\nCOPULA:  "     
   cp <- "  theta = "; as.p <- x$theta.a
   ct <- "  tau = "; kt.p <- x$tau.a
   s1 <- "sigma2.1 = "; s1.p <- x$sigma21.a
   s2 <- "sigma2.2 = "; s2.p <- x$sigma22.a 
   n1 <- "nu.1 = "; n1.p <- x$nu1.a
   n2 <- "nu.2 = "; n2.p <- x$nu2.a 
   
 
   if(x$margins[1] %in% c("N","GU","rGU","LO","GAi") )  m1l <- "identity"
   if(x$margins[2] %in% c("N","GU","rGU","LO","GAi") )  m2l <- "identity"
   if(x$margins[1] %in% c("LN","WEI","iG","GA","DAGUM","SM","FISK") ) m1l <- "log" 
   if(x$margins[2] %in% c("LN","WEI","iG","GA","DAGUM","SM","FISK") ) m2l <- "log" 
   if(x$margins[1] %in% c("BE") )                              m1l <- "qlogis" 
   if(x$margins[2] %in% c("BE") )                              m2l <- "qlogis"    
   

  
  cat(main.t,cop) 
  
  if(x$margins[1]=="N")      cat("\nMARGIN 1: Gaussian")  
  if(x$margins[1]=="GU")     cat("\nMARGIN 1: Gumbel")    
  if(x$margins[1]=="rGU")    cat("\nMARGIN 1: reverse Gumbel")  
  if(x$margins[1]=="LO")     cat("\nMARGIN 1: logistic")   
  if(x$margins[1]=="LN")     cat("\nMARGIN 1: log-normal") 
  if(x$margins[1]=="WEI")    cat("\nMARGIN 1: Weibull") 
  if(x$margins[1]=="iG")     cat("\nMARGIN 1: inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))     cat("\nMARGIN 1: gamma")  
  if(x$margins[1]=="BE")     cat("\nMARGIN 1: beta")    
  if(x$margins[1]=="DAGUM")  cat("\nMARGIN 1: Dagum")  
  if(x$margins[1]=="SM")     cat("\nMARGIN 1: Singh-Maddala") 
  if(x$margins[1]=="FISK")     cat("\nMARGIN 1: Fisk") 
  
  

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
  if(x$margins[2]=="FISK")   cat("\nMARGIN 2: Fisk") 
  
  
  
  
  cat("\n\nEQUATION 1")    
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]]) # fg <- x$formula1; fg[[2]] <- as.symbol("mu.1"); print(fg) 
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    


  

   cat("\nEQUATION 2") 
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]]) # fg <- x$formula2; fg[[2]] <- as.symbol("mu.2"); print(fg)   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    




  if( x$X3.null == FALSE ){
  
  
  
     if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par ){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma2.1:","log","\n") else cat("\nLink function for sigma2.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
    cat("\n")
  
      if(x$l.sp3!=0){
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
      }
    
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2.2:","log","\n") else cat("\nLink function for sigma2.2:","qlogis","\n")  
  cat("Formula: "); print(x$formula[[4]]) 
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp4!=0){
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp5!=0){
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  


}









     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par ){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma2.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2.2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
  cat("\nEQUATION 6")
  cat("\nLink function for nu.2:","log","\n") 
  cat("Formula: "); print(x$formula[[6]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 7")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp7!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}












     if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par ){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma2.1:","log","\n") else cat("\nLink function for sigma2.1:","qlogis","\n")
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2.2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
 
     
  cat("\nEQUATION 5")
  cat("\nLink function for nu.2:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}








     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par ){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma2.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Smooth components' approximate significance:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2.2:","log","\n") else cat("\nLink function for sigma2.2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
 
  
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Smooth components' approximate significance:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}







}






  CIrs    <- colMeans(x$CItheta, na.rm = TRUE)
  CIkt    <- colMeans(x$CIkt, na.rm = TRUE)
  CIsig21 <- colMeans(x$CIsig21, na.rm = TRUE)
  CIsig22 <- colMeans(x$CIsig22, na.rm = TRUE)
  
  
  if(x$margins[1] %in% cont3par)  CInu1 <- colMeans(x$CInu1, na.rm = TRUE)
  if(x$margins[2] %in% cont3par)  CInu2 <- colMeans(x$CInu2, na.rm = TRUE)  

  nodi <- 3
  
  
  if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par ) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par ) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                     "  ",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")",   
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par ) cat(  s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")", 
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
                          
  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par ) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
       
invisible(x)
                
}





















