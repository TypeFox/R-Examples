print.SemiParTRIVProbit <- function(x, ...){

  cop <- "Gaussian"; lind <- "atanh"
  
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

  cat("\nEQUATION 2")
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  
  cat("\nEQUATION 3")
  cat("\nLink function for mu.3:",m3l,"\n")
  cat("Formula: "); print(x$formula[[3]])  
  
  cat("\n")
            
  cat("n = ",x$n,"\ntheta12 = ", format(as.p12, digits=3),"  theta13 = ", format(as.p13, digits=3),"  theta23 = ", format(as.p23, digits=3),"\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
 
invisible(x)

}

