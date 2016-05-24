fitcml <- function(mt_ind, nrlist, x_mt, rtot, W, ngroups, gind, x_mtlist, NAstruc, g_NA, st.err, etaStart, gby){

  #cml function for call in nlm
  cml <- function(eta){

    beta <- as.vector(W %*% eta)
    #FIXME!!! gby??
    beta.list <- split(beta,gind)  #gind index for treatment groups
    beta.list1 <- beta.list
   
    #beta and NAstructure (over Groups): 1st line parameter values, 2nd line which item NA 
    betaNA <- mapply(function(x,y) {rbind(x,y)},beta.list1,NAstruc,SIMPLIFY=FALSE)  
   
    #likelihood term based on gamma functions for each Group x NAgroup combination
    Lg <- lapply(betaNA, function(betaNAmat) {      
      beta.vec <- betaNAmat[1,]                #get parameter vector beta
   
      #gamma functions for each NAgroup within Groups 
      Lg.NA <- apply(matrix(betaNAmat[-1,],ncol=length(beta.vec)),1, function(NAvec) {
        
        #list of virtual item-category parameters per item
        beta_list <- as.list(split(beta.vec[NAvec==1],mt_ind[1:(length(beta.vec[NAvec==1]))]))       
        parlist <- lapply(beta_list,exp)                                #initial epsilon as list
   
               #------------------gamma functions----------------------
               g_iter <- NULL                                                  #computation of the gamma functions
               K <- length(parlist)
               for (t in 1:(K-1)) {                                            #building up J1,...,Jt,...,Js
   
                 if (t==1) {                                                   #first iteration step
                   gterm <- c(1,parlist[[t]])                                  #0th element included
                 }else
                 {
                  gterm <- g_iter                                   #gamma previous iteration with 0th el
                  g_iter <- NULL
                 }
   
                 parvek <- c(1,parlist[[t+1]])                      #eps vector in current iteration with 0th el
                 h <- length(parvek)                                #dimensions for matrix
                 mt <- length(gterm)
                 rtot1 <- h+mt-1                                    #number of possible raw scores (0 included)
   
                 gtermvek <- rep(c(gterm,rep(0,h)),h)                          #building up matrix for gamma term
                 gtermvek <- gtermvek[-((length(gtermvek)-h+1):length(gtermvek))]      #eliminating last h 0's
                 gmat <- matrix(gtermvek,nrow=rtot1,ncol=h)
                 emat <- matrix(rep(parvek,rep(rtot1,h)),ncol=h,nrow=rtot1)    #building up matrix for eps term
                 gmat_new <- gmat*emat                                                 #merge matrices
                 g_iter <- rowSums(gmat_new)                     #gamma functions in current iteration are rowsums
               }
              #----------------- end gamma functions ------------------
   
              Lg.NA <- as.vector(g_iter[2:(rtot+1)])     #final gamma vector stored in gamma (without gamma0)
              return(Lg.NA)
              })
    })
    #----------------- compute likelihood components -----------------------
    L1 <- sum(mapply(function(x,z) {
                      x[!is.na(z)]%*%na.exclude(z)
                      },nrlist,lapply(Lg,log)))        #sum up L1-terms (group-wise)
   
    L2 <- sum(mapply("%*%",x_mtlist,beta.list1))        #sum up L2-terms (group-wise)
    L1-L2                                               #final likelihood value
  }
  #----------------- end likelihood -----------------------

  eta <- etaStart                                     #starting values for eta parameters

  if(!exists("fitctrl")) fitctrl <- "nlm"   # if fitctrl undefined, set it to "nlm"   ### MjM 2014-01-27

  if(fitctrl == "nlm"){
    suppressWarnings(
      fit <- nlm(cml, eta, hessian = st.err, iterlim = 5000L)   # NLM optimizer
    )
  } else if(fitctrl == "optim"){
    suppressWarnings( 
      fit <- optim(eta, cml, method = "BFGS", hessian = TRUE, control = list(maxit = 5000L))
    )
    fit$counts<-fit$counts[1]
    names(fit)<-c("estimate","minimum","iterations","code","message","hessian")
  } else stop("optimizer misspecified in fitctrl\n")

  fit

}
