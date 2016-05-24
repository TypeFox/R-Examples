 pen.log.like <- function(penden.env,cal=FALSE,temp=FALSE) {

  if(cal) {
   assign("pen.log.like",sum(sapply(get("f.hat.val",penden.env),log))-0.5*(t(get("ck.val",penden.env))%*%get("DDD.sum",penden.env)%*%get("ck.val",penden.env)),penden.env)
    assign("log.like",sum(sapply(get("f.hat.val",penden.env),log)),penden.env)
  }
  if(temp) {
    assign("pen.log.like.temp",sum(sapply(get("f.hat.val.temp",penden.env),log))-0.5*(t(get("ck.val.temp",penden.env))%*%get("DDD.sum",penden.env)%*%get("ck.val.temp",penden.env)),penden.env)
    assign("log.like.temp",sum(sapply(get("f.hat.val.temp",penden.env),log)),penden.env)
 }
 
}
