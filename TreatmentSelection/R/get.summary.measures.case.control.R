get.summary.measures.case.control <-
function(data, rho, d=0){  
  #rho <- rho[1]
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
  #proportion marker negative
    p.marker.neg <- mean(neg[event==1])*rho[3] + mean(neg[event==0])*(1-rho[3])
    p.marker.pos<- mean(pos[event==1])*rho[3] + mean(pos[event==0])*(1-rho[3])

 #Average Benefit (B) of no treatment among marker negatives...
 
    #empirical estimate

    B.neg.emp.trt0 <-  ifelse( sum(event[trt==0 & neg==1])>0, expit( logit(mean(event[trt==0 & neg==1])) + logit(rho[3]) - logit(mean(event)) ), 0)
    B.neg.emp.trt1 <-  ifelse( sum(event[trt==1 & neg==1])>0, expit( logit(mean(event[trt==1 & neg==1])) + logit(rho[3]) - logit(mean(event)) ), 0) 

    B.neg.emp <- B.neg.emp.trt1 - B.neg.emp.trt0

    #model based estimate
    n.r   <- sum(event)
    n.nor <- sum(1-event)
    #numerator
    B.neg.mod.num <-     (1/n.r)*(rho[3])*ifelse(sum(neg==1 & event==1)>0, sum(-trt.effect[neg==1 & event==1]), 0) + 
                     (1/n.nor)*(1-rho[3])*ifelse(sum(neg==1 & event==0)>0, sum(-trt.effect[neg==1 & event==0]), 0)
    #denominator
    B.neg.mod.den <-     (1/n.r)*(rho[3])*sum(neg[event==1]) + 
                     (1/n.nor)*(1-rho[3])*sum(neg[event==0])

    B.neg.mod <- ifelse(B.neg.mod.den>0, B.neg.mod.num/B.neg.mod.den, 0)
 
    #Average Benefit (B) of treatment among marker positives...
 
    #empirical estimate
    B.pos.emp.trt0 <-  ifelse( sum(event[trt==0 & pos==1])>0, expit( logit(mean(event[trt==0 & pos==1])) + logit(rho[3]) - logit(mean(event)) ), 0)
    B.pos.emp.trt1 <-  ifelse( sum(event[trt==1 & pos==1])>0, expit( logit(mean(event[trt==1 & pos==1])) + logit(rho[3]) - logit(mean(event)) ), 0) 

    B.pos.emp <- B.pos.emp.trt0 - B.pos.emp.trt1

    #model based estimate
    #numerator
    B.pos.mod.num <-     (1/n.r)*(rho[3])*ifelse(sum(pos==1 & event==1)>0, sum(trt.effect[pos==1 & event==1]), 0) + 
                     (1/n.nor)*(1-rho[3])*ifelse(sum(pos==1 & event==0)>0, sum(trt.effect[pos==1 & event==0]), 0)
    #denominator
    B.pos.mod.den <-     (1/n.r)*(rho[3])*sum(pos[event==1]) +
                     (1/n.nor)*(1-rho[3])*sum(pos[event==0])

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
       #wi[event==1] <- rho[1]/(mean(event ==1)*Pr.S)
       #wi[event==0] <- (1-rho[1])/(mean(event==0)*Pr.S)  

       #if(sum(wi > 25) > 0 ) { warning("Warning: there exist sampling weights > 25, Variance Estimate may be unstable")}

     
   #  Var.Delta = mean((trt.effect[event==1])^2)*rho[3] + mean((trt.effect[event==0])^2)*(1-rho[3]) 
                   
                  # mean(trt.effect[event==1]^2)*rho[3] + mean(trt.effect[event==0]^2)*(1-rho[3]) - 
                  #(mean(trt.effect[event==1])*rho[3]   + mean(trt.effect[event==0])*(1-rho[3]))^2  #sum((wi*(trt.effect- d )^2))/sum(wi)

    #Total Gain (TG)


    delta.F <- get.F.case.control( trt.effect, event, trt, rho = rho)
    
    ooo <- order(trt.effect)

    s.delta.F <- delta.F[ooo]
    p0 <- ((mean(1-trt[event==1]))*rho[3])/(1-rho[2])
    p1 <- ((mean(trt[event==1]))*rho[3])/(rho[2])
    TG <- sum( diff(c(0, s.delta.F))*abs( sort(trt.effect) - (p0 - p1) ))
    
  #event rates
  ER.trt0.emp = ifelse( sum(trt==0)>0, expit( logit(mean(event[trt==0])) + logit(rho[3]) - logit(mean(event)) ), 0)
  ER.trt0.mod = mean(data$fittedrisk.t0[event==1])*rho[3] + mean(data$fittedrisk.t0[event==0])*(1-rho[3])

  ER.trt1.emp = ifelse( sum(trt==1)>0, expit( logit(mean(event[trt==1])) + logit(rho[3]) - logit(mean(event)) ), 0)
  ER.trt1.mod = mean(data$fittedrisk.t1[event==1])*rho[3] + mean(data$fittedrisk.t1[event==0])*(1-rho[3])
  
  if(is.null(data[["marker.pos"]])){
    #default is trt all
    ER.mkrbased.emp = ER.trt1.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt1.mod - Theta.mod
  }else{
    ER.mkrbased.emp = ER.trt0.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt0.mod - Theta.mod    
  }
  
 
 
 Var.Delta = mean((trt.effect[event==1])^2)*rho[3] + mean((trt.effect[event==0])^2)*(1-rho[3]) - (ER.trt0.emp - ER.trt1.emp)^2
 
  
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
