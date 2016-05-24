EventRateVec <- function(risk.t0, risk.t1, f.d, rho = NULL, event=NULL, trt=NULL){
  #calculate event rate for v percentile of delta 
  # calculate this for all study designs, which we can tell based on the length of rho 

   N = length(f.d)
   order_fd <- order(f.d)
   
  if(all(rho==-9999) | is.null(rho)){

    cs_riskT0 <- cumsum(risk.t0[order_fd])/seq_along(f.d)
    cs_riskT1 <- rev(cumsum(rev(risk.t1[order_fd]))/seq_along(f.d)) #'reverse' cumsum
    
    wght <- sort(f.d)/100
    
    #plot(sort(f.d), cs_riskT0, ylim = c(0,1))
    #points(sort(f.d), cs_riskT1, col = "red" )
    
    #plot(f.d, risk.t0, ylim = c(0,1))
    #points(f.d, risk.t1, col = "red" )
  
    out <- (cs_riskT0*wght+ cs_riskT1*(1-wght))
    
    #plot(f.d, out[rank(f.d)])

  }else if(all(rho[5:7]==-9999)){ #case cohort
    ## need to use cutoffs and calculate sums across different time points 
    ## because I am taking averages here and need to weight by Pr(D==1) in the cohort
    
    ## recall rho[3]  is Pr(D==1) in cohort
    #sample weights
    cc.adjust = numeric(N)
    cc.adjust[event==1] = rho[3]/(mean(event==1)*N/rho[1])
    cc.adjust[event==0] = (1- rho[3])/(mean(event==0)*N/rho[1])

    sorted_risk.t0 <- risk.t0[order_fd]
    sorted_risk.t1 <- risk.t1[order_fd]
    sorted_weights = cc.adjust[order_fd]

    cs_riskT0 <- cumsum(sorted_risk.t0*(sorted_weights))/cumsum(sorted_weights)
    cs_riskT1 <- rev(cumsum(rev(sorted_risk.t1*sorted_weights))/cumsum(rev(sorted_weights))) #'reverse' cumsum

    wght <- sort(f.d)/100
    
    out <- (cs_riskT0*wght+ cs_riskT1*(1-wght))
    

    
  }else{ #stratified case cohort
    n.cc <- length(event)
    Pr.S <- n.cc/rho[1] # Pr(S=1)
    
    Pr.D1.trt1 <- rho[5]
    Pr.D0.trt1 <- rho[4]
    Pr.D1.trt0 <- rho[3]
    Pr.D0.trt0 <- rho[2]
    
    #sampling weights for scc design
    scc.adjust <- numeric(N)
    
    scc.adjust[event==1 & trt==1] <- Pr.D1.trt1/(mean(event==1 & trt==1)*Pr.S)
    scc.adjust[event==1 & trt==0] <- Pr.D1.trt0/(mean(event==1 & trt==0)*Pr.S)              
    scc.adjust[event==0 & trt==1] <- Pr.D0.trt1/(mean(event==0 & trt==1)*Pr.S)
    scc.adjust[event==0 & trt==0] <- Pr.D0.trt0/(mean(event==0 & trt==0)*Pr.S)
    
    sorted_risk.t0 <- risk.t0[order_fd]
    sorted_risk.t1 <- risk.t1[order_fd]
    sorted_weights = scc.adjust[order_fd]
    
    cs_riskT0 <- cumsum(sorted_risk.t0*(sorted_weights))/cumsum(sorted_weights)
    cs_riskT1 <- rev(cumsum(rev(sorted_risk.t1*sorted_weights))/cumsum(rev(sorted_weights))) #'reverse' cumsum
    
    wght <- sort(f.d)/100
    
    out <- (cs_riskT0*wght+ cs_riskT1*(1-wght))

    
  }
   #have to resort by the original data.
    out[rank(f.d)]
}