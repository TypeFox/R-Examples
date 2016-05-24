get.summary.measures.stratified.case.control <-
function(data, rho, d=0){  
  
  event <- data$event
  trt <- data$trt
  trt.effect <- data$trt.effect

  if(is.null(data[["marker.pos"]])){
    neg <- data$marker.neg
    pos <- 1 - neg
  }else{
    pos <- data$marker.pos
    neg <- 1 - pos
  }
  
   
    # get probabilites for all strata based on rho1 and rho2
    Pr.D1.trt1 <- rho[5]
    Pr.D0.trt1 <- rho[4]
    Pr.D1.trt0 <- rho[3]
    Pr.D0.trt0 <- rho[2]

    Pr.D1.givT0 <- rho[3]/(rho[2]+rho[3])
    Pr.D1.givT1 <- rho[5]/(rho[4]+rho[5])
    #proportion marker negative
    p.marker.neg <- mean(neg[event==1 & trt==1])*Pr.D1.trt1 + 
                    mean(neg[event==0 & trt==1])*Pr.D0.trt1 +
                    mean(neg[event==1 & trt==0])*Pr.D1.trt0 +
                    mean(neg[event==0 & trt==0])*Pr.D0.trt0
  
    p.marker.pos <- mean(pos[event==1 & trt==1])*Pr.D1.trt1 + 
                    mean(pos[event==0 & trt==1])*Pr.D0.trt1 +
                    mean(pos[event==1 & trt==0])*Pr.D1.trt0 +
                    mean(pos[event==0 & trt==0])*Pr.D0.trt0
  
    #Average Benefit (B) of no treatment among marker negatives...
 
    #empirical estimate
    #cat(paste("\n", mean(event[trt==0 & neg ==1])))
    #if(!is.finite(mean(event[trt==0 & neg ==1]))) browser()
    B.neg.emp.trt0 <- ifelse(sum(trt==0 & neg==1)>0, expit( logit(mean(event[trt==0 & neg ==1])) + logit(Pr.D1.givT0) - logit(mean(event[trt==0])) ), 0)
    B.neg.emp.trt1 <- ifelse(sum(trt==1 & neg==1)>0, expit( logit(mean(event[trt==1 & neg ==1])) + logit(Pr.D1.givT1) - logit(mean(event[trt==1])) ), 0)

    B.neg.emp <- B.neg.emp.trt1 - B.neg.emp.trt0

    #model based estimate
    n.11 <- sum((event==1)*(trt==1))
    n.01 <- sum((event==0)*(trt==1))
    n.10 <- sum((event==1)*(trt==0))
    n.00 <- sum((event==0)*(trt==0))

    #numerator
    B.neg.mod.num <- (1/n.11)*Pr.D1.trt1*ifelse(sum(neg==1 & event==1 & trt==1)>0, sum(-trt.effect[neg==1 & event==1 & trt==1]), 0) + 
                     (1/n.01)*Pr.D0.trt1*ifelse(sum(neg==1 & event==0 & trt==1)>0, sum(-trt.effect[neg==1 & event==0 & trt==1]), 0) + 
                     (1/n.10)*Pr.D1.trt0*ifelse(sum(neg==1 & event==1 & trt==0)>0, sum(-trt.effect[neg==1 & event==1 & trt==0]), 0) + 
                     (1/n.00)*Pr.D0.trt0*ifelse(sum(neg==1 & event==0 & trt==0)>0, sum(-trt.effect[neg==1 & event==0 & trt==0]), 0) 
 
    #denominator
    B.neg.mod.den <- (1/n.11)*Pr.D1.trt1*sum(neg[event==1 & trt==1]) + 
                     (1/n.01)*Pr.D0.trt1*sum(neg[event==0 & trt==1]) + 
                     (1/n.10)*Pr.D1.trt0*sum(neg[event==1 & trt==0]) + 
                     (1/n.00)*Pr.D0.trt0*sum(neg[event==0 & trt==0])    

    B.neg.mod <- ifelse(B.neg.mod.den>0, B.neg.mod.num/B.neg.mod.den, 0)
 
    #Average Benefit (B) of treatment among marker positives...
    #empirical estimate
     
    B.pos.emp.trt0 <- ifelse(sum(trt==0 & pos==1)>0, expit( logit(mean(event[trt==0 & pos ==1])) + logit(Pr.D1.givT0) - logit(mean(event[trt==0])) ), 0)
    B.pos.emp.trt1 <- ifelse(sum(trt==1 & pos==1)>0, expit( logit(mean(event[trt==1 & pos ==1])) + logit(Pr.D1.givT1) - logit(mean(event[trt==1])) ), 0)

    B.pos.emp <- B.pos.emp.trt0 - B.pos.emp.trt1

    #model based estimate
    n.11 <- sum((event==1)*(trt==1))
    n.01 <- sum((event==0)*(trt==1))
    n.10 <- sum((event==1)*(trt==0))
    n.00 <- sum((event==0)*(trt==0))

    #numerator
    B.pos.mod.num <- (1/n.11)*Pr.D1.trt1*ifelse(sum(pos==1 & event==1 & trt==1)>0, sum(trt.effect[pos==1 & event==1 & trt==1]), 0) + 
                     (1/n.01)*Pr.D0.trt1*ifelse(sum(pos==1 & event==0 & trt==1)>0, sum(trt.effect[pos==1 & event==0 & trt==1]), 0) + 
                     (1/n.10)*Pr.D1.trt0*ifelse(sum(pos==1 & event==1 & trt==0)>0, sum(trt.effect[pos==1 & event==1 & trt==0]), 0) + 
                     (1/n.00)*Pr.D0.trt0*ifelse(sum(pos==1 & event==0 & trt==0)>0, sum(trt.effect[pos==1 & event==0 & trt==0]), 0) 
 
    #denominator
    B.pos.mod.den <- (1/n.11)*Pr.D1.trt1*sum(pos[event==1 & trt==1]) + 
                     (1/n.01)*Pr.D0.trt1*sum(pos[event==0 & trt==1]) + 
                     (1/n.10)*Pr.D1.trt0*sum(pos[event==1 & trt==0]) + 
                     (1/n.00)*Pr.D0.trt0*sum(pos[event==0 & trt==0])    

    B.pos.mod <- ifelse(B.pos.mod.den>0, B.pos.mod.num/B.pos.mod.den, 0)


    #Theta

  
  if(is.null(data[["marker.pos"]])){
    Theta.emp <- B.neg.emp*p.marker.neg
    Theta.mod <- B.neg.mod*p.marker.neg
    
  }else{
    Theta.emp <- B.pos.emp*p.marker.pos
    Theta.mod <- B.pos.mod*p.marker.pos
  }
  

    # Variance of Trt Effects
    
      # calculate sample weights 
       n.cc <- length(event)
       Pr.S <- n.cc/rho[1] # Pr(S=1)
 
       #wi <- rep(0, length(event))

       #wi[event==1 & trt==1] <- Pr.D1.trt1/(mean(event==1 & trt==1)*Pr.S)
       #wi[event==1 & trt==0] <- Pr.D1.trt0/(mean(event==1 & trt==0)*Pr.S)              
       #wi[event==0 & trt==1] <- Pr.D0.trt1/(mean(event==0 & trt==1)*Pr.S)
       #wi[event==0 & trt==0] <- Pr.D0.trt0/(mean(event==0 & trt==0)*Pr.S)
      
       #if(sum(wi > 25) > 0 ) { warning("Warning: there exist sampling weights > 25, Variance Estimate may be unstable")}
        

 
#sum((wi*(trt.effect-d)^2))/sum(wi)

                
                   # mean(trt.effect[event==1 & trt==1]^2)*Pr.D1.trt1 + 
                  # mean(trt.effect[event==0 & trt==1]^2)*Pr.D0.trt1 +
                  # mean(trt.effect[event==1 & trt==0]^2)*Pr.D1.trt0 +
                  # mean(trt.effect[event==0 & trt==0]^2)*Pr.D0.trt0 -
                  # (mean(trt.effect[event==1 & trt==1])*Pr.D1.trt1 + 
                  # mean(trt.effect[event==0 & trt==1])*Pr.D0.trt1 +
                  # mean(trt.effect[event==1 & trt==0])*Pr.D1.trt0 +
                  # mean(trt.effect[event==0 & trt==0])*Pr.D0.trt0)^2 #sum((wi*(trt.effect-d)^2))/sum(wi)

    #Total Gain (TG)

    delta.F <- get.F.stratified.case.control( trt.effect, event, trt, rho = rho)
    
    ooo <- order(trt.effect)

    s.delta.F <- delta.F[ooo]

    TG <- sum( diff(c(0, s.delta.F))*abs( sort(trt.effect) - (Pr.D1.givT0 - Pr.D1.givT1) ) )
    
  #event rates
  ER.trt0.emp = ifelse(sum(trt==0)>0, expit( logit(mean(event[trt==0])) + logit(Pr.D1.givT0) - logit(mean(event[trt==0])) ), 0)
  ER.trt0.mod = mean(data$fittedrisk.t0[event==1 & trt==1])*Pr.D1.trt1 + 
                mean(data$fittedrisk.t0[event==0 & trt==1])*Pr.D0.trt1 +
                mean(data$fittedrisk.t0[event==1 & trt==0])*Pr.D1.trt0 +
                mean(data$fittedrisk.t0[event==0 & trt==0])*Pr.D0.trt0
  
  ER.trt1.emp = ifelse(sum(trt==1)>0, expit( logit(mean(event[trt==1])) + logit(Pr.D1.givT1) - logit(mean(event[trt==1])) ), 0)
  ER.trt1.mod = mean(data$fittedrisk.t1[event==1 & trt==1])*Pr.D1.trt1 + 
                mean(data$fittedrisk.t1[event==0 & trt==1])*Pr.D0.trt1 +
                mean(data$fittedrisk.t1[event==1 & trt==0])*Pr.D1.trt0 +
                mean(data$fittedrisk.t1[event==0 & trt==0])*Pr.D0.trt0
  
  if(is.null(data[["marker.pos"]])){
    #default is trt all
    ER.mkrbased.emp = ER.trt1.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt1.mod - Theta.mod
  }else{
    ER.mkrbased.emp = ER.trt0.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt0.mod - Theta.mod    
  }
  

Var.Delta =   mean(trt.effect[event==1 & trt==1]^2)*Pr.D1.trt1 + 
  mean(trt.effect[event==0 & trt==1]^2)*Pr.D0.trt1 +
  mean(trt.effect[event==1 & trt==0]^2)*Pr.D1.trt0 +
  mean(trt.effect[event==0 & trt==0]^2)*Pr.D0.trt0 -(ER.trt0.emp - ER.trt1.emp)^2

  
  list(     p.neg = p.marker.neg,
            p.pos = p.marker.pos,
            B.neg.emp = B.neg.emp,
            B.neg.mod = B.neg.mod, 
            B.pos.emp = B.pos.emp, 
            B.pos.mod = B.pos.mod, 
            Theta.emp = Theta.emp, 
            Theta.mod = Theta.mod, 
            Var.Delta = Var.Delta, 
            TG        = TG, 
            ER.trt0.emp = ER.trt0.emp, 
            ER.trt0.mod = ER.trt0.mod, 
            ER.trt1.emp = ER.trt1.emp, 
            ER.trt1.mod = ER.trt1.mod, 
            ER.mkrbased.emp = ER.mkrbased.emp,
            ER.mkrbased.mod = ER.mkrbased.mod
  ) 
  
  
}
