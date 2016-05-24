#' @title Simulation routine for CUBE models without covariates
#' @aliases cubeforsim
#' @description Fit CUBE models without covariates to given ordinal data. It is useful for simulation experiments 
#' since it performs the same steps as the CUBE function, but with no printed output.
#' @usage cubeforsim(m, ordinal, starting, maxiter = 500, toler = 1e-06, expinform=FALSE)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param starting Vector of initial parameters estimates to start the optimization algorithm
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' (defaul: maxiter = 500)
#' @param toler Fixed error tolerance for final estimates (default: toler = 1e-6)
#' @param expinform Logical: if TRUE, the function returns the expected variance covariance matrix
#'  (default is FALSE)
#' @return An object of the class "CUBE", with null output for $BIC since the routine is only for
#'  simulation purposes
#' @export cubeforsim
#' @import stats
#' @keywords htest
#' @seealso \code{\link{CUBE}}, \code{\link{loglikCUBE}}
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,37])
#' starting<-rep(0.1,3) 
#' simul<-cubeforsim(m,ordinal,starting)
#' param<-simul$estimates    # Estimated parameters vector (pai,csi,phi)
#' ####################
#' data(univer)
#' m<-7
#' ordinal<-univer[,11]
#' starting<-rep(0.1,3)
#' simul<-cubeforsim(m,ordinal,starting)
#' param<-simul$estimates   # Estimated parameters vector (pai,csi,phi) 


cubeforsim <-
function(m,ordinal,starting,maxiter=500,toler=1e-6,expinform=FALSE){
  freq<-tabulate(ordinal,nbins=m); n<-sum(freq); 
  aver<-mean(ordinal)
  #(0)# initial estimates
  pai<-starting[1]; csi<-starting[2]; phi<-starting[3];
  #(0)# log-lik
  loglik<-loglikcube(m,freq,pai,csi,phi)
  # ********************************************************************
  # *************   E-M algorithm for CUBE     *************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    
    #(1)# betar
    bb<-betar(m,csi,phi)
    aa<-(1-pai)/(m*pai*bb)
    
    #(2)# taunor
    tauno<-1/(1+aa)
    
    #(3)# pai(k+1)
    pai<-sum(freq*tauno)/n  # updated pai estimate
    
    paravecjj<-c(csi,phi)
    
    #(4)# Q(k+1)
    dati<-cbind(tauno,freq)
    ################ EFFECUBE is Q(csi,phi) ###########################
    
    #(5)# (csi(k+1),phi(k+1))
    paravec<-paravecjj
    ################################## maximize w.r.t. paravec  ########
    
    ### prova simulazione ###
    #optimestim=optim(paravec,effecube,dati=dati,m=m,method = "L-BFGS-B",lower=c(0.01,0.0001),upper=c(0.99,0.5)) #prima era c(0.99,0.3))
    optimestim<-optim(paravec,effecube,dati=dati,m=m,method = "L-BFGS-B",lower=c(0.01,0.01),upper=c(0.99,0.5)) #prima era c(0.99,0.3))
    ################################################################         
    
    #(6)# theta(k+1)
    paravecjj<-optimestim$par   # updated paravec estimates
    csi<-paravecjj[1];   phi<-paravecjj[2];
    
    ##########################################
    if(pai<0.001){pai<-0.001; nniter<-maxiter-1}
    #       if(csi<0.001){csi=0.001; nniter=maxiter-1}
    #       if(phi<0.001){phi=0.001; nniter=maxiter-1}
    if(pai>0.999){pai<-0.99}         ### to avoid division by 0 !!!
    #       if(csi>0.999){csi=0.99}         ### to avoid division by 0 !!!
    ###################################### print(c(nniter,pai,csi,phi))
    
    # if(phi>1000){phi=1000; nniter=maxiter-1}    ### to avoid not-convergence !!!
    
    #(7)# elle(theta(k+1))
    liknew<-loglikcube(m,freq,pai,csi,phi)
    
    #(8)# test
    testll<-abs(liknew-likold)            # OPTIONAL printing: print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi,phi));
    if(testll<=toler) break else {loglik<-liknew} # OPTIONAL printing: print(loglik);
    nniter<-nniter+1
  }
  loglik<-liknew
  ###### End of E-M algorithm for CUBE ***********************************************
  
  # ********************************************************
  # Compute ML var-cov matrix and print result for CUBE
  # ********************************************************
  if (expinform==FALSE){
    varmat<-varcovcubeobs(m,pai,csi,phi,freq)
  } else {
    varmat<-varcovcubeexp(m,pai,csi,phi,n)
  }
  
  #####################################################################
  # Assignments as global variables: assign('name',value,pos=1)
  #####################################################################
  #   assign('pai',pai,pos=1)
  #   assign('csi',csi,pos=1)
  #   assign('phi',phi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  #   assign('nniter',nniter,pos=1)
  stime<-c(pai,csi,phi)
  results<-list('estimates'=round(stime,digits=5), 'loglik'=loglik, 'niter'=nniter, 'varmat'=varmat)
  
  ####################################
}
