calculate <- function(penden.env) {
  f.hat.val(penden.env,cal=TRUE)
  if(get("no",penden.env)) {
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("AIC",0,penden.env)
    assign("cAIC",0,penden.env)
    assign("BIC",0,penden.env)
    class(penden.env) <- "pencopula"
    return(penden.env)
  }
  liste <- get("liste",penden.env)
  pen.log.like(penden.env,cal=TRUE)
  Derv1(penden.env)
  Derv2(penden.env)
  p <- get("p",penden.env)
  DD <- get("DD",penden.env)
  assign("i",i <- 1,penden.env)
  liste[i,1] <- get("pen.log.like",penden.env)
  liste[i,2] <- get("log.like",penden.env)
  liste[i,(4:(4+p-1))] <- get("lambda",penden.env)
  liste[i,((4+p):(4+p+DD-1))] <- get("ck.val",penden.env)

  assign("liste",liste,penden.env)
  
  assign("calc",TRUE,penden.env)
 
  if(new.weights(penden.env)=="fehler"){
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("AIC",0,penden.env)
    assign("cAIC",0,penden.env)
    assign("BIC",0,penden.env)
    #obj <- list(penden.env=penden.env)
    class(penden.env) <- "pencopula"
    return(penden.env)
  }

  my.loop(penden.env)
  
  Derv1(penden.env)
  Derv2(penden.env)

  my.IC(penden.env)
}
