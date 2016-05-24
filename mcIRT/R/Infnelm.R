Infnelm <- function(ESTlist, fromto=c(-5,5), gran=200)
{

  ALPHAS  <- ESTlist$ZLpar$alpha   
  BETAS   <- ESTlist$ZLpar$beta    
  NRMS    <- ESTlist$ZLpar$nrmpar  
  

  thetas <- seq(fromto[1], fromto[2], length.out=gran)
 
  
  catinfG <- mapply(function(levs, alphaG, betaG, nrmG)
  { # loops all groups
    

    
    catinfI <- mapply(function(alI, betI, nrmI)
    { # loops all items
      
      pitem <- c(nrmI$zetas,nrmI$lambdas)
      Km  <- matrix(c(rep(1,length(thetas)), thetas),ncol=2)
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      # 2PL
      abpar <- c(betI,alI)
      P2pl <- twoplpart(Km=Km, abpar=abpar)
      Q2pl <- (1-P2pl)
      
      TwoPLInf <- alI^2 * P2pl * P2pl * Q2pl # sihe Psychometrika Artikel S.462
      #####
      #### Normally the item information function is:
      # alpha^2 * P * Q
      # but for the Nested logit model it is alpha^2 * P * P * Q --> because you are modeling the distractors as well.
      # A similar thing is true for the nrm - part. as we know, the nrm is modeled in case the correct answer was not found --> so it depends on 1-P !
      # --> this can be seen in line labels with: # ***
      

      LAMs <- pitem[(length(pitem)/2 + 1):length(pitem)]
      
      W_g <- sapply(1:nrow(ZQstern),function(zei) # geht die nodes durch
      {
        Zqrow <- ZQstern[zei,]
        z     <- thetas[zei]
        Pqrep <- matrix(-Zqrow, length(Zqrow), length(Zqrow)) 
        diag(Pqrep) <- 1-Zqrow 
        Pdi <- diag(Zqrow)
        Wi  <- Pqrep %*% Pdi
        as.vector(LAMs %*% Wi %*% LAMs) * Q2pl[zei] * ZQstern[zei,] # ***
        
      }) 
      

      
      
      cbind(TwoPLInf,t(W_g)) # Category Informations incl 2PL part
      
    },alI = alphaG, betI = betaG, nrmI = nrmG,  SIMPLIFY = F) 
    
    
    catinfI
    
  },levs=levels(ESTlist$reshOBJ$gr), alphaG=ALPHAS ,betaG = BETAS, nrmG = NRMS ,SIMPLIFY = FALSE)
  
# this approach does not work with unequal number of categories 
#   TIFall <- lapply(catinfG,function(GRs)
#               {
#               apply(simplify2array(GRs, higher=TRUE),1,sum)
#               })

  TIFall <- lapply(catinfG,function(GRs)
    {
      rowSums(do.call("cbind",GRs))
    })
  
  # category informations - for different thetas for each group
  class(catinfG) <- "infnlm"
  
  
  return(list(catinfG=catinfG, thetas=thetas, TestInfGROUPS=TIFall))
  
}  


