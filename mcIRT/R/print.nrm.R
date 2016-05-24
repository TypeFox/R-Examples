print.nrm <-
function(x, ...)
{

  RESnrm <- x
  rm(x)
  
  # 2 not (-2) because it is fnscale=-1 in optim  
  min2logL <- 2*RESnrm$last_mstep$value 
  
  
  # Parameters <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  
  names1a <- lapply(1:length(RESnrm$reshOBJ$aDD),function(x)
  {
    actu <- RESnrm$reshOBJ$aDD[[x]]
    paste("Item",x, "|categ",1:actu$anz_cat,sep="")  
  })
  
  
  form1a <- mapply(function(eachG,eachSE)
  {
    forEg <- do.call("rbind",lapply(eachG,function(x)matrix(x,ncol=2)))  
    rownames(forEg) <- unlist(names1a)
    colnames(forEg) <- c("zeta","lambda")
    
    forSE <- do.call("rbind",lapply(eachSE,function(x)matrix(x,ncol=2))) 
    colnames(forSE) <- c("SE|zeta","SE|lambda")
    
    allto <- cbind(forEg,forSE)[,c(1,3,2,4)]
    
    allto
  },eachG=RESnrm$ZLpar,eachSE=attr(RESnrm$SE,"listform"),SIMPLIFY=FALSE)
  
  
  whereM1 <- which(colSums(RESnrm$reshOBJ$Qmat) > 1)
  whereM1

  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  
  firstpart <- matrix(c(min2logL,as.integer(RESnrm$n_steps),ncol(RESnrm$reshOBJ$Qmat)))
  rownames(firstpart) <- c("-2logLikelihood:","Number of EM-cycles:","Number of estimated parameters: ")
  colnames(firstpart) <- ""
  
  ### number of parameters
  if(!RESnrm$ctrl$nonpar)
  {
    ### number of parameters
    nme  <- length(RESnrm$erg_distr$mean_est) - 1
    #RESnrm$ctrl$sigmaest
    nva  <- RESnrm$ctrl$sigmaest * (length(RESnrm$erg_distr$sig_est) -1)
    npar <- ncol(RESnrm$reshOBJ$Qmat) + nme + nva  - length(RESnrm$ctrl$Clist)
    
  } else 
      {
        pardist <- length(RESnrm$QUAD$A$nodes)*length(RESnrm$QUAD) - 3 - (length(RESnrm$QUAD) - 1)  
        # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
        npar <- ncol(RESnrm$reshOBJ$Qmat) + pardist - length(RESnrm$ctrl$Clist)
      }
  
  ######### OUTPUT:
  

  cat("\n Call:",deparse(RESnrm$call),"\n\n")
  
  cat("\n Global Informations")
  cat("\n -------------------------------------------------------------------- \n")
  
  cat("\n-2logLikelihood: ",min2logL,"\n")
  cat("Number of EM-cycles: ",RESnrm$n_steps,"\n")
  cat("Number of estimated parameters: ", npar ,"\n")
  cat("Estimation design:\n")
  print(RESnrm$reshOBJ$design)
  

  cat("\n\n Category Parameter estimates and SE")
  cat("\n -------------------------------------------------------------------- \n")
  print(form1a)
  
  
  
}
