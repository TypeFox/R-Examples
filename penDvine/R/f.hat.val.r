f.hat.val <- function(penden.env,cal=FALSE,temp=FALSE) {
  if(cal) {
    fit <- get("tilde.PSI.d.D",penden.env)%*%get("ck.val",penden.env)
    if(all(fit>=0)) {
      assign("f.hat.val",fit,penden.env)
    }
    else{
      assign("no",TRUE,penden.env)
      return(paste("d=",get("d",penden.env),"D=",get("D",penden.env),"lambda=",get("lambda",penden.env)[1],sep=""))
    }
  }
  if(temp) {
    fit <- get("tilde.PSI.d.D",penden.env)%*%get("ck.val.temp",penden.env)
    if(all(fit>=0)) assign("f.hat.val.temp",fit,penden.env)
    else{
      assign("no",TRUE,penden.env)
      return(paste("d=",get("d",penden.env),"D=",get("D",penden.env),"lambda=",get("lambda",penden.env)[1],sep=""))
    }
  }
}
