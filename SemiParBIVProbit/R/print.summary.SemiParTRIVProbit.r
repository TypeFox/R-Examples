print.summary.SemiParTRIVProbit <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

nodi <- 3
cop  <- "Gaussian"
lind <- "atanh"
  
  as.p12 <- x$theta12.a
  as.p13 <- x$theta13.a
  as.p23 <- x$theta23.a
  
      main.t <- "\nCOPULA:  "     

      m1l <- m2l <- m3l <- "probit"
      
      cat(main.t,cop) 
      cat("\nMARGIN 1: Bernoulli")  
      cat("\nMARGIN 2: Bernoulli")
      cat("\nMARGIN 3: Bernoulli")       
      
  
  cat("\n\nEQUATION 1")  
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]]) 
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
  cat("Formula: "); print(x$formula[[2]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
  cat("\nEQUATION 3")
  cat("\nLink function for mu.3:",m3l,"\n")
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
    
    
    
    
  


  #CI12 <- colMeans(x$CI12s, na.rm = TRUE)
  #CI13 <- colMeans(x$CI13s, na.rm = TRUE)
  #CI23 <- colMeans(x$CI23s, na.rm = TRUE)




  CI12 <- x$CI12s
  CI13 <- x$CI13s
  CI23 <- x$CI23s




  
  
  
cat("\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),
                 "\ntheta12 = ",format(as.p12,digits=nodi),"(",format(CI12[1],digits=nodi),",",format(CI12[2],digits=nodi),")",
                 "\ntheta13 = ",format(as.p13,digits=nodi),"(",format(CI13[1],digits=nodi),",",format(CI13[2],digits=nodi),")",
                 "\ntheta23 = ",format(as.p23,digits=nodi),"(",format(CI23[1],digits=nodi),",",format(CI23[2],digits=nodi),")",                 
                 "\n\n", sep="")  

       
invisible(x)
                
}


