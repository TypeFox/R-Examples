my.loop <- function(penden.env) {
  DD <- get("DD",penden.env)
  n.liste <- matrix(0,1,dim(get("liste",penden.env))[2])
  assign("i",i <- 2,penden.env)
  max.iter <- get("max.iter",penden.env)
  p <- get("p",penden.env)
  #fix.lambda <- get("fix.lambda",penden.env)
   while(get("calc",penden.env)) {
    convert <- logical()
    
    if(i==2) old.ck <- get("liste",penden.env)[(i-1),((4+p):(4+p+DD-1))] + 1e-08   else old.ck <- get("liste",penden.env)[(i-1),((4+p):(4+p+DD-1))]
     if(all(abs((get("ck.val.temp",penden.env)/old.ck)-1)<0.01) | i-1>max.iter) {
      convert <- TRUE
      #if(!fix.lambda) assign("lambda",get("lambda.save",penden.env),penden.env)
    }
    else convert <- FALSE

    if(convert) {
      #assign("liste",liste,penden.env)
      break
    }
    else {
      if(get("no",penden.env)){
        #assign("pen.log.like",0,penden.env)
        #assign("log.like",0,penden.env)
        #assign("AIC",0,penden.env)
        #assign("BIC",0,penden.env)
        break
      }
      else assign("f.hat.val",get("f.hat.val.temp",penden.env),penden.env)
      assign("ck.val",get("ck.val.temp",penden.env),penden.env)
      assign("pen.log.like",get("pen.log.like.temp",penden.env),penden.env)
      #if(!fix.lambda) assign("marg.log.like",get("marg.log.like.temp",penden.env),penden.env)
      assign("log.like",get("log.like.temp",penden.env),penden.env)
      #if(!fix.lambda) assign("lambda",get("lambda.temp",penden.env),penden.env)
      penalty.matrix(penden.env)

      n.liste <- c(get("pen.log.like",penden.env),get("log.like",penden.env),0,get("lambda",penden.env),get("ck.val",penden.env))
      assign("liste",rbind(get("liste",penden.env),n.liste),penden.env)
      #liste[i,1] <- get("pen.log.like",penden.env)
      #liste[i,2] <- get("log.like",penden.env)
      #if(!fix.lambda) liste[i,3] <- get("marg.log.like",penden.env)
      #liste[i,(4:(4+p-1))] <- get("lambda",penden.env)
      #liste[i,((4+p):(4+p+DD-1))] <- get("ck.val",penden.env)
 
      Derv1(penden.env)

      Derv2(penden.env)
 
      if(new.weights(penden.env)=="quadprog fehler") {
        #assign("pen.log.like",0,penden.env)
        #assign("log.like",0,penden.env)
        #assign("AIC",0,penden.env)
        #assign("BIC",0,penden.env)
        break       
      }
      
      #if(!fix.lambda) {
      #  new.lambda(penden.env)
      #  penalty.matrix(penden.env,temp=TRUE)
      #  pen.log.like(penden.env,temp.lambda=TRUE)
      #  marg.likelihood(penden.env,temp=TRUE)
      #}
      assign("i",i <- i+1,penden.env)
      #assign("liste",liste,penden.env)
    }  
  }
}
