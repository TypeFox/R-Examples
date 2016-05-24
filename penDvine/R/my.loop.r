my.loop <- function(penden.env) {
  DD <- get("DD",penden.env)
  eps <- 0.01
  p <- get("p",penden.env)
  n.liste <- matrix(0,1,3+DD+1)
  assign("i",i <- 2,penden.env)
  max.iter <- get("max.iter",penden.env)
  assign("calc",TRUE,penden.env)
  fix.lambda <- get("fix.lambda",penden.env)
  assign("lambda.out",FALSE,penden.env)

  while(get("calc",penden.env)) {
    #browser()
    if(get("wrong.lambda",penden.env)){
      assign("calc",FALSE,penden.env)
      f.hat.val(penden.env,cal=TRUE)
      if(get("no",penden.env)) break
      pen.log.like(penden.env,cal=TRUE)
      Derv1(penden.env)
      Derv2(penden.env)
      marg.likelihood(penden.env,get("pen.log.like",penden.env))
      my.IC(penden.env)
      n.liste <- c(get("pen.log.like",penden.env),get("log.like",penden.env),get("marg.log.like",penden.env),get("lambda",penden.env),get("cAIC",penden.env),get("ck.val",penden.env))
      list <- rbind(get("liste",penden.env),n.liste)
      rownames(list) <- seq(0,(i-1))
      assign("liste",list,penden.env)
      break
    }
    if(get("no",penden.env)) break
    assign("f.hat.val",get("f.hat.val.temp",penden.env),penden.env)
    old.ck <- get("ck.val",penden.env)
    assign("ck.val",get("ck.val.temp",penden.env),penden.env)
    assign("log.like",get("log.like.temp",penden.env),penden.env)
    assign("pen.log.like",get("pen.log.like.temp",penden.env),penden.env)
    if(!fix.lambda) {
      help.lambda <- new.lambda(penden.env)
      if((abs(help.lambda-get("lambda",penden.env))<eps*get("lambda",penden.env))|(i-1)>max.iter|get("lambda.out",penden.env)) {
      #if(all(abs((get("ck.val.temp",penden.env)/old.ck)-1)<0.005) | i-1>max.iter) {
        assign("calc",FALSE,penden.env)
        assign("lambda",help.lambda,penden.env)
        f.hat.val(penden.env,cal=TRUE)
        if(get("no",penden.env)) break
        pen.log.like(penden.env,cal=TRUE)
        Derv1(penden.env)
        Derv2(penden.env)
        marg.likelihood(penden.env,get("pen.log.like",penden.env))
        my.IC(penden.env)
        n.liste <- c(get("pen.log.like",penden.env),get("log.like",penden.env),get("marg.log.like",penden.env),get("lambda",penden.env),get("cAIC",penden.env),get("ck.val",penden.env))
        list <- rbind(get("liste",penden.env),n.liste)
        rownames(list) <- seq(0,(i-1))
        assign("liste",list,penden.env)
        break
      }
      else {
        f.hat.val(penden.env,temp=TRUE)
        if(get("no",penden.env)) break
        pen.log.like(penden.env,temp=TRUE)
        Derv1(penden.env,temp=TRUE,lambda=help.lambda)
        Derv2(penden.env,temp=TRUE,lambda=help.lambda)
        marg.likelihood(penden.env,get("pen.log.like",penden.env),temp=TRUE)
        my.IC(penden.env,temp=TRUE)
        assign("lambda",help.lambda,penden.env)
        new.weights(penden.env,lambda.temp=help.lambda)
        n.liste <- c(get("pen.log.like.temp",penden.env),get("log.like.temp",penden.env),get("marg.log.like.temp",penden.env),help.lambda,get("cAIC.temp",penden.env),get("ck.val",penden.env))
        list <- rbind(get("liste",penden.env),n.liste)
        rownames(list) <- seq(0,(i-1))
        assign("liste",list,penden.env)
      }
    }
    if(fix.lambda) {
      if(all(abs((get("ck.val.temp",penden.env)/old.ck)-1)<0.005) | i-1>max.iter) {
        f.hat.val(penden.env,cal=TRUE)
        if(get("no",penden.env)) break
        pen.log.like(penden.env,cal=TRUE)
        Derv1(penden.env)
        Derv2(penden.env)
        marg.likelihood(penden.env,get("pen.log.like",penden.env))
        my.IC(penden.env)
        n.liste <- c(get("pen.log.like",penden.env),get("log.like",penden.env),get("marg.log.like",penden.env),get("lambda",penden.env),get("cAIC",penden.env),get("ck.val",penden.env))
        list <- rbind(get("liste",penden.env),n.liste)
        rownames(list) <- seq(0,(i-1))
        assign("liste",list,penden.env)
        break
      }
      else {
        f.hat.val(penden.env,temp=TRUE)
        if(get("no",penden.env)) break
        pen.log.like(penden.env,temp=TRUE)
        Derv1(penden.env,temp=TRUE,lambda=get("lambda",penden.env))
        Derv2(penden.env,temp=TRUE,lambda=get("lambda",penden.env))
        marg.likelihood(penden.env,get("pen.log.like",penden.env),temp=TRUE)
        my.IC(penden.env,temp=TRUE)
        n.liste <- c(get("pen.log.like.temp",penden.env),get("log.like.temp",penden.env),get("marg.log.like.temp",penden.env),get("lambda",penden.env),get("cAIC.temp",penden.env),get("ck.val",penden.env))
        list <- rbind(get("liste",penden.env),n.liste)
        rownames(list) <- seq(0,(i-1))
        assign("liste",list,penden.env)
        new.weights(penden.env,lambda.temp=get("lambda",penden.env))
      }
    }
    assign("i",i <- i+1,penden.env)
  }
}
