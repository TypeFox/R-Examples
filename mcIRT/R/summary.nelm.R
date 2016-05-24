summary.nelm <-
function(object, RETURN=FALSE, ...)
{
#browser()
  RESnlm <- object
  # 2 not (-2) because it is fnscale=-1 in optim  
  min2logL <- 2*RESnlm$last_mstep$value 
  
  
  # Parameters <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  

  
  # alpha and beta parameters
  albeP <- mapply(function(AL,BE,eachSE)
            {
            albepar <- cbind(AL,BE)
            
            extra <- sapply(eachSE,function(xyz)
                        {
                        woal <- grep("alpha",names(xyz))
                       xyz[woal]
                        })
                      
            extrb <- sapply(eachSE,function(xyz)
            {
              wobe <- grep("beta",names(xyz))
              xyz[wobe]
            })
            
            cbind(albepar,extra,extrb)[,c(1,3,2,4)]
            },AL=RESnlm$ZLpar$alpha,BE=RESnlm$ZLpar$beta,eachSE=attr(RESnlm$SE,"listform"),SIMPLIFY=FALSE)

  
  albePm   <- do.call("rbind",albeP)
  innerlev <- rep(paste("Item",1:length(RESnlm$ZLpar$alpha[[1]]),sep=""),nlevels(RESnlm$reshOBJ$gr))
  backdran <- rep(levels(RESnlm$reshOBJ$gr),each=length(RESnlm$ZLpar$alpha[[1]]))
  rownames(albePm) <- paste(innerlev,"|",backdran,sep="")
  colnames(albePm) <- c("alpha","SE|alpha","beta","SE|beta")
  
  
  


catnrgroup <- lapply(levels(RESnlm$reshOBJ$gr),function(nixi)
{
  
  catallIT <- as.vector(mapply(function(x,y)
  {
    #paste0("Item",y, "|categ",1: (x$anz_cat-1))
    paste0("Item",y, "|categ",gsub(".*(\\d{1,}).*","\\1",x$categ[-1]))
  },x=RESnlm$reshOBJ$aDD,y=1:length(RESnlm$reshOBJ$aDD)))
  
  catallIT 
  
})



  
  form1a <- mapply(function(eachG,eachSE,cnrg)
  {
    forEg <- do.call("rbind",lapply(eachG,function(x)matrix(unlist(x),ncol=2)))  
    rownames(forEg) <- unlist(cnrg)
    colnames(forEg) <- c("zeta","lambda")
    
    # extract
    seextract <- lapply(eachSE,function(x){
      woalbe <- grep("(alpha|beta)",names(x))
      matrix(x[-woalbe],ncol=2)
      })
    
    forSE <- do.call("rbind",seextract) 
    colnames(forSE) <- c("SE|zeta","SE|lambda")
    
    allto <- cbind(forEg,forSE)[,c(1,3,2,4)]
    
    allto
  },eachG=RESnlm$ZLpar$nrmpar,eachSE=attr(RESnlm$SE,"listform"), cnrg=catnrgroup,SIMPLIFY=FALSE)
  

  # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  
  Mest <- RESnlm$erg_distr$mean_est
  Sest <- RESnlm$erg_distr$sig_est
  meansig <- rbind(Mest,Sest)
  
  rownames(meansig) <- c("mean","sigma^2")
  colnames(meansig) <- paste("group|",levels(RESnlm$reshOBJ$gr),sep="")
  
    SEmat <- RESnlm$erg_distr$errmat
  rownames(SEmat) <- c("SE|mean","SE|sigma^2")
  colnames(SEmat) <- paste("group|",levels(RESnlm$reshOBJ$gr),sep="")
  

  
  ### number of parameters <--- NEW
  nme  <- length(RESnlm$erg_distr$mean_est) - 1
  nva  <- RESnlm$ctrl$sigmaest *(length(RESnlm$erg_distr$sig_est) -1)
  npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva - length(RESnlm$ctrl$Clist)
  
  
  if(!RESnlm$ctrl$nonpar)
  {
    
    ### number of parameters <--- NEW
    nme  <- length(RESnlm$erg_distr$mean_est) - 1
    nva  <- RESnlm$ctrl$sigmaest *(length(RESnlm$erg_distr$sig_est) -1)
    npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva - length(RESnlm$ctrl$Clist)
    
  } else 
  {
    pardist <- length(RESnlm$QUAD[[1]]$nodes)*length(RESnlm$QUAD) - 3 - (length(RESnlm$QUAD) - 1)  
    # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
    npar <- ncol(RESnlm$reshOBJ$Qmat) + pardist - length(RESnlm$ctrl$Clist)
  }
  
  
  nonparametric <- ifelse(RESnlm$ctrl$nonpar,"nonparametric","parametric")
  
  ### number of parameters
#   nme  <- length(RESnlm$erg_distr$mean_est) - 1
#   nva  <- length(RESnlm$erg_distr$sig_est) - 1
#   npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva

  firstpart <- matrix(c(min2logL,as.integer(RESnlm$n_steps),npar))
  rownames(firstpart) <- c("-2logLikelihood:","Number of EM-cycles:","Number of estimated parameters: ")
  colnames(firstpart) <- ""
  
  
  ######### OUTPUT:
  
  cat("\n Call:",deparse(RESnlm$call),"\n- job started @",attr(RESnlm$call,"date"),"\n\n\n")
  
  cat("\n Global Informations")
  cat("\n -------------------------------------------------------------------- \n")
  
  print(firstpart)
  
  cat(">>",attr(RESnlm$call,"convergence"),"<<")
  
  cat("\n\n Parameter estimates for latent distributions")
  cat("\n -------------------------------------------------------------------- \n")
  cat("\n Point estimators:\n")
  print(meansig)
  cat("\n Standard Errors:\n")
  print(SEmat)
  
  cat("\nPrior:", nonparametric , "&", attr(RESnlm$QUAD,"wherefrom"))  
  
  cat("\n\n Parameter estimates for the 2-PL part")
  cat("\n -------------------------------------------------------------------- \n")
  print(albePm)
  
  
  
  cat("\n\n Category Parameter estimates and SE")
  cat("\n -------------------------------------------------------------------- \n")
  print(form1a)
  
if(RETURN)  
{
return(list(firstpart=firstpart,meansig=meansig,SEmat=SEmat,albePm=albePm,form1a=form1a))
  
  
}

}
