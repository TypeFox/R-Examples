print.nelm <-
function(x, ...)
{
  RESnlm <- x
  rm(x)
  #
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
  
  
#   catnrgroup <- lapply(RESnlm$reshOBJ$recm, function(grou)
#   {
#     
#     catnr <- mapply(function(x,numbers) 
#     {
#       cn <- colnames(x)[-1]
#       cnn <- gsub(".*_(\\d+)","\\1",cn, perl=TRUE)
#       paste("Item",numbers, "|categ",cnn,sep="") 
#       
#     },x=grou, numbers=1:length(RESnlm$reshOBJ$aDD),SIMPLIFY=FALSE)
#     catnr  
#   })
  
#   catnrgroup <- lapply(levels(RESnlm$reshOBJ$gr),function(nixi)
#   {
#     
#     catallIT <- as.vector(mapply(function(x,y)
#     {
#       paste0("Item",y, "|categ",1: (x$anz_cat-1))
#     },x=RESnlm$reshOBJ$aDD,y=1:length(RESnlm$reshOBJ$aDD)))
#     
#     catallIT 
#     
#   })
  
  
catnrgroup <- lapply(levels(RESnlm$reshOBJ$gr),function(nixi)
{
  
  catallIT <- as.vector(mapply(function(x,y)
  {
    #paste0("Item",y, "|categ",1: (x$anz_cat-1))
    paste0("Item",y, "|categ",gsub(".*(\\d{1,}).*","\\1",x$categ[-1]))
  },x=RESnlm$reshOBJ$aDD,y=1:length(RESnlm$reshOBJ$aDD)))
  
  catallIT 
  
})



  #   
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
  
  firstpart <- matrix(c(min2logL,as.integer(RESnlm$n_steps),ncol(RESnlm$reshOBJ$Qmat)))
  rownames(firstpart) <- c("-2logLikelihood:","Number of EM-cycles:","Number of estimated parameters: ")
  colnames(firstpart) <- ""
  
  ### number of parameters
#   nme  <- length(RESnlm$erg_distr$mean_est) - 1
#   nva  <- length(RESnlm$erg_distr$sig_est) - 1
#   npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva
  
  ### number of parameters <--- NEW
  if(!RESnlm$ctrl$nonpar)
  {
    
    ### number of parameters <--- NEW
    nme  <- length(RESnlm$erg_distr$mean_est) - 1
    nva  <- RESnlm$ctrl$sigmaest *(length(RESnlm$erg_distr$sig_est) -1)
    npar <- ncol(RESnlm$reshOBJ$Qmat) + nme + nva - length(RESnlm$ctrl$Clist)
    
  } else 
      {
        pardist <- length(RESnlm$QUAD$A$nodes)*length(RESnlm$QUAD) - 3 - (length(RESnlm$QUAD) - 1)  
        # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
        npar <- ncol(RESnlm$reshOBJ$Qmat) + pardist - length(RESnlm$ctrl$Clist)
      }
  
  
  ######### OUTPUT:
  
  cat("\n Call:",deparse(RESnlm$call),"\n- job started @",attr(RESnlm$call,"date"),"\n\n\n")
  
  cat("\n Global Informations")
  cat("\n -------------------------------------------------------------------- \n")
  
  cat("\n-2logLikelihood:\t \t \t",min2logL,"\n")
  cat("Number of EM-cycles: \t",RESnlm$n_steps,"\n")
  cat("Number of estimated parameters:", npar ,"\n")
     cat("Estimation design:\n")
      print(RESnlm$reshOBJ$design)


  cat("\n\n Parameter estimates for the 2-PL part")
  cat("\n -------------------------------------------------------------------- \n")
  print(albePm)
  
  
  
  cat("\n\n Category Parameter estimates and SE")
  cat("\n -------------------------------------------------------------------- \n")
  print(form1a)
  
  
}
