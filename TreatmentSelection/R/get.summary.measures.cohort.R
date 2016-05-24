get.summary.measures.cohort <-
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
    #proportion marker negative
    p.marker.neg <- mean(neg)
    p.marker.pos <- mean(pos)
    #Average Benefit (B) of no treatment among marker negatives...
 
    #empirical estimate
    B.neg.emp <- ifelse(length(event[trt==0 & neg==1]) > 0 & length(event[trt==1 & neg==1]) > 0, 
                        mean(event[trt==1 & neg==1]) - mean(event[trt==0 & neg==1]) ,
                       0)

    #model based estimate
    B.neg.mod <- ifelse(sum(neg) >0, -mean(trt.effect[neg==1]), 0)
 
    #Average Benefit (B) of treatment among marker positives...
 
    #empirical estimate
    B.pos.emp <- ifelse(length(event[trt==0 & pos==1]) > 0 & length(event[trt==1 & pos==1]) > 0, 
                        mean(event[trt==0 & pos==1]) - mean(event[trt==1 & pos==1]), 
                       0)

    #model based estimate
    B.pos.mod <- ifelse(sum(pos)>0, mean(trt.effect[pos==1]), 0)

    #Theta
  
  if(is.null(data[["marker.pos"]])){
    Theta.emp <- B.neg.emp*p.marker.neg
    Theta.mod <- B.neg.mod*p.marker.neg
    
  }else{
    Theta.emp <- B.pos.emp*p.marker.pos
    Theta.mod <- B.pos.mod*p.marker.pos
  }
  
    p0.hat <- mean(event[trt==0])
    p1.hat <- mean(event[trt==1])
    Var.Delta <- mean((trt.effect - (p0.hat - p1.hat))^2)

    delta.F <- get.F.cohort( trt.effect, event, trt, rho = NULL)
    
    ooo <- order(trt.effect)

    s.delta.F <- delta.F[ooo]
    
    TG <- sum( diff(c(0, s.delta.F))*abs( sort(trt.effect) - (p0.hat - p1.hat))) 

  #event rates
  ER.trt0.emp = mean(event[trt==0])
  ER.trt0.mod = mean(data$fittedrisk.t0)
  ER.trt1.emp = mean(event[trt==1])
  ER.trt1.mod = mean(data$fittedrisk.t1)
  
  if(is.null(data[["marker.pos"]])){
    #default is trt all
    ER.mkrbased.emp = ER.trt1.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt1.mod - Theta.mod
  }else{
    ER.mkrbased.emp = ER.trt0.emp - Theta.emp 
    ER.mkrbased.mod = ER.trt0.mod - Theta.mod    
  }
    
    
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
