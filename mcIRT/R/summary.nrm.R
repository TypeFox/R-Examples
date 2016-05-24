summary.nrm <-
function(object, RETURN=FALSE, ...)
{

  RESnrm <- object
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


Mest <- RESnrm$erg_distr$mean_est
Sest <- RESnrm$erg_distr$sig_est
meansig <- rbind(Mest,Sest)

rownames(meansig) <- c("mean","sigma^2")
colnames(meansig) <- paste("group|",levels(RESnrm$reshOBJ$gr),sep="")

SEmat <- RESnrm$erg_distr$errmat
rownames(SEmat) <- c("SE|mean","SE|sigma^2")
colnames(SEmat) <- paste("group|",levels(RESnrm$reshOBJ$gr),sep="")


  
if(!RESnrm$ctrl$nonpar)
{
  
### number of parameters
nme  <- length(RESnrm$erg_distr$mean_est) - 1
#RESnrm$ctrl$sigmaest
nva  <- RESnrm$ctrl$sigmaest * (length(RESnrm$erg_distr$sig_est) -1)
npar <- ncol(RESnrm$reshOBJ$Qmat) + nme + nva  - length(RESnrm$ctrl$Clist)
  
} else 
    {
      pardist <- length(RESnrm$QUAD[[1]]$nodes)*length(RESnrm$QUAD) - 3 - (length(RESnrm$QUAD) - 1)  
      # anzahl der bins - 3 für die erste gruppe und anzahl - 1 für die restlichen gruppen + itpar - constants
      npar <- ncol(RESnrm$reshOBJ$Qmat) + pardist - length(RESnrm$ctrl$Clist)
    }

nonparametric <- ifelse(RESnrm$ctrl$nonpar,"nonparametric","parametric")
  
firstpart <- matrix(c(min2logL,as.integer(RESnrm$n_steps),npar))
rownames(firstpart) <- c("-2logLikelihood:","Number of EM-cycles:","Number of estimated parameters:")
colnames(firstpart) <- ""
  
  
######### OUTPUT:

cat("\n Call:",deparse(RESnrm$call),"\n- job started @",attr(RESnrm$call,"date"),"\n\n\n")

cat("\n Global Informations")
cat("\n -------------------------------------------------------------------- \n")


print(firstpart)

cat(">>",attr(RESnrm$call,"convergence"),"<<")

cat("\n\n Parameter estimates for latent distributions")
cat("\n -------------------------------------------------------------------- \n")
cat("\n Point estimators:\n")
print(meansig)
cat("\n Standard Errors:\n")
print(SEmat)
  
cat("\nPrior:", nonparametric , "&", attr(RESnrm$QUAD,"wherefrom"))  
  
cat("\n\n Category Parameter estimates and SE")
cat("\n -------------------------------------------------------------------- \n")
print(form1a)


if(RETURN)  
{
  return(list(firstpart=firstpart,meansig=meansig,SEmat=SEmat,form1a=form1a))
 
}

}
