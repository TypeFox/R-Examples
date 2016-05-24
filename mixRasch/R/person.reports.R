`person.reports` <-
function(LC, n.c, ud, treat.extreme, max.iter, conv.crit, steps, as.LCA){

  matcher <- match(ud$people,ud$x.x) 
  matcher.d <- match(ud$people,ud$degen.x.x)

  for(c in 1:n.c){
    namer <- c("theta", "SE.theta", "r", "infit", "in.Z", "outfit", "out.Z") 
    suppressWarnings(LC[[c]]$person.par <- data.frame(cbind(LC[[c]]$person.par$theta[matcher],
                                             LC[[c]]$person.par$SE.theta[matcher], ud$r[matcher],
                                             LC[[c]]$person.par$in.out[matcher,])) )    
    colnames(LC[[c]]$person.par) <- namer 

  if(length(ud$degen.r) > 0){
    skiprest <- FALSE
    ud$degen.r <- ifelse(ud$degen.r == 0, ud$degen.r + treat.extreme, ud$degen.r - treat.extreme)
    degen.theta <- rep(0,length(ud$degen.r))
    if(! as.LCA) for(its in 1:max.iter) {
     odt <- degen.theta
     degen.theta <- deg.theta(LC[[c]], ud$degen.r, ud$degen.x.i, degen.theta, steps)
     if(any(is.na(degen.theta))) { warning("Extreme person estimates did not converge \n")
                                   skiprest <- TRUE
                                   break}
     if(abs(max(odt - degen.theta)) < conv.crit) break
    }  
  if(! skiprest){
   LC[[c]]$person.par$theta[! is.na(matcher.d)] <- degen.theta[na.omit(matcher.d)] 

   Pxji <- array(apply(LC[[c]]$item.par$delta,2,P.xj, th=degen.theta),
                dim=c(steps,length(degen.theta),LC[[c]]$i.stat$n.i))
   LC[[c]]$person.par$SE.theta[! is.na(matcher.d)] <- (-1*d.v(Pxji, ud$degen.r, 
                                matrix(! is.na(ud$degen.x.i), nrow=length(ud$degen.r)))$d2[na.omit(matcher.d)])^(-.5)
  }
  }     
 }

LC
}

