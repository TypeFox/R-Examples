Infnrm <- function(ESTlist, fromto=c(-5,5), gran=200)
{
  # ESTlist --> internal ESTlist created by nrm. must contain the estimated values which stem from the EM procedure as well as the centered category parameters
  

  
  #relstv <- relist(opp,SKEL)
  relstv  <- ESTlist$ZLpar
  thetas <- seq(fromto[1], fromto[2], length.out=gran)
  

  
  catinfG <- mapply(function(levs,stvl)
  { # loops all groups
    
    # pitem = stvl[[1]]
    catinfI <- mapply(function(pitem)
    { # loops all items
      
      Km  <- matrix(c(rep(1,length(thetas)), thetas),ncol=2)
      LAM <- matrix(pitem,nrow=2,byrow=T)
      
      Z <- Km %*% LAM
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      

      LAMs <- pitem[(length(pitem)/2 + 1):length(pitem)]
      
      W_g <- sapply(1:nrow(ZQstern),function(zei) # geht die nodes durch
      {
        Zqrow <- ZQstern[zei,]
        z     <- thetas[zei]
        Pqrep <- matrix(-Zqrow, length(Zqrow), length(Zqrow)) 
        diag(Pqrep) <- 1-Zqrow 
        Pdi <- diag(Zqrow)
        Wi  <- Pqrep %*% Pdi
        as.vector(LAMs %*% Wi %*% LAMs) * ZQstern[zei,]
        
      }) 
      
     t(W_g) # Category Information
      
    },pitem=stvl,  SIMPLIFY = F) ### was T before
    
    catinfI
  },levs=levels(ESTlist$reshOBJ$gr), stvl=relstv ,SIMPLIFY = FALSE)
  
#browser()
# GRs <- catinfG[[1]]
#   TIFall <- lapply(catinfG,function(GRs)
#   {
#     apply(simplify2array(GRs, higher=TRUE),1,sum)
#   })
  
    TIFall <- lapply(catinfG,function(GRs)
      {
        rowSums(do.call("cbind",GRs))
      })
      
  # category informations - for different thetas for each group
  class(catinfG) <- "infnrm"
  
  
  return(list(catinfG=catinfG, thetas=thetas, TestInfGROUPS=TIFall))
  
}  



