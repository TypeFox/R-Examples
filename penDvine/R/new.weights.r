new.weights <- function(penden.env,lambda.temp=NULL) {
  dd <- get("dd",penden.env)
  p <- get("p",penden.env)
  DD <- get("DD",penden.env)
  #l.A <- length(get("A.Restrict",penden.env)[,1,1])
  l.A <- length(get("T.marg",penden.env)[,1,1])
  vec <- seq(1,(l.A-1))
 
  eps <- 1e-08
  if(get("base",penden.env)=="Bernstein") assign("AA.help",t(get("T.marg",penden.env)[vec,,1]),penden.env)
  if(get("base",penden.env)=="B-spline") assign("AA.help",t(get("A.Restrict",penden.env)[vec,,1]),penden.env)
             
  for(j in 2:p) {
    if(get("base",penden.env)=="Bernstein") assign("AA.help",cbind(get("AA.help",penden.env),t(get("T.marg",penden.env)[vec,,j])),penden.env)
    if(get("base",penden.env)=="B-spline") assign("AA.help",cbind(get("AA.help",penden.env),t(get("A.Restrict",penden.env)[vec,,j])),penden.env)
  }

  meq <- 1+p*(length(vec))
  if(get("base",penden.env)=="B-spline") {
    bvec <- c(rep(0,1+p*length(vec)),-get("ck.val",penden.env)+rep(eps,length(get("ck.val",penden.env))))
    assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),diag(1,DD)),penden.env)
  } 
  if(get("base",penden.env)=="Bernstein") {
    #bvec <- c(rep(0,1+p*length(vec)),-get("ck.val",penden.env)+rep(eps,length(get("ck.val",penden.env))))
    bvec <- c(0,rep(1/(get("K",penden.env)+1),2*length(vec))-t(get("AA.help",penden.env))%*%get("ck.val",penden.env),-get("ck.val",penden.env)+rep(eps,length(get("ck.val",penden.env))))
    assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),diag(1,DD)),penden.env)
  }
  Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
  Derv2(penden.env,temp=TRUE,lambda=lambda.temp)
  aa <- try(obj <- solve.QP(Dmat=-get("Derv2.pen.temp",penden.env),dvec=get("Derv1.pen.temp",penden.env),Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE))
  if(class(aa)=="try-error") {
    assign("wrong.lambda",TRUE,penden.env)
    return("fehler")
  }
  if(any(is.na(obj$solution))&get("i",penden.env)==1) {
     h <- 1
     calc <- TRUE
     lambda.seq <- seq(1,51,length=11)
      while(calc) {
     	lambda.temp <- lambda.seq[h]
        Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
        Derv2(penden.env,temp=TRUE,lambda=lambda.temp)
        aa <- try(obj <- solve.QP(Dmat=-get("Derv2.pen.temp",penden.env),dvec=get("Derv1.pen.temp",penden.env),Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE))
        if(class(aa)=="try-error"|any(is.na(obj$solution))) h <- h+1
        else {
          assign("lambda",lambda.seq[h],penden.env)
          calc <- FALSE
          assign("lambda.change",TRUE,penden.env)
          assign("wrong.lambda",FALSE,penden.env)
        }
     if(h>9) {
        assign("wrong.lambda",TRUE,penden.env)
        return("fehler")
     }
     }
  }
  else {
     assign("lambda.change",FALSE,penden.env)
     assign("wrong.lambda",FALSE,penden.env)
  }
  if(get("i",penden.env)>1& any(is.na(obj$solution))) {
     #browser()
     print("is.na in obj$solution")
     assign("wrong.lambda",TRUE,penden.env)
     return("fehler")
  }
  assign("delta",obj$solution,penden.env)
  assign("ck.val.temp",get("ck.val",penden.env)+obj$solution,penden.env)
  Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
  Derv2(penden.env,temp=TRUE,lambda=lambda.temp)

  f.hat.val(penden.env,temp=TRUE)
  pen.log.like(penden.env,temp=TRUE)
  if(get("no",penden.env)) return("fehler")
  else return("kein fehler")
}
