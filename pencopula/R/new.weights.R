new.weights <- function(penden.env) {
  dd <- get("dd",penden.env)
  p <- get("p",penden.env)
  DD <- get("DD",penden.env)
  calc <- TRUE
  l.A <- length(get("A.Restrict",penden.env)[,1,1])
  if(get("base",penden.env)=="B-spline") vec <- seq(1,(l.A-1))
  if(get("base",penden.env)=="Bernstein") vec <- seq(1,(l.A-1))
  assign("AA.help",t(get("A.Restrict",penden.env)[vec,,1]),penden.env)          
  for(j in 2:p) {
     assign("AA.help",cbind(get("AA.help",penden.env),t(get("A.Restrict",penden.env)[vec,,j])),penden.env)
  }
  prob.val2 <- 1
  meq <- 1+p*(length(vec))

  if(get("base",penden.env)=="Bernstein") {
    bvec <- c(rep(0,1+p*length(vec)),-get("ck.val",penden.env)+rep(1e-05,length(get("ck.val",penden.env))))
    assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),diag(1,DD)),penden.env)
    aa <- try(obj <- solve.QP(Dmat=-get("Derv2.pen",penden.env),dvec=get("Derv1.pen",penden.env),Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE))
    if(class(aa)=="try-error") return("fehler")
    if(any(is.na(obj$solution))) return("fehler")
  }
  
  if(get("base",penden.env)=="B-spline") {
    if(!get("adapt.grid",penden.env)) {
      assign("val",get("tilde.PSI.d.D.knots.start.g.all",penden.env)%*%get("ck.val",penden.env),penden.env)
    #assign("ind.help",seq(1,dim(get("tilde.PSI.d.D.knots.start.g.all",penden.env))[1]),penden.env)
      if(get("add",penden.env)) bvec <- c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g.all",penden.env)%*%(get("ck.val",penden.env)+rep(1e-06,length(get("ck.val",penden.env)))))
      else  bvec <- c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g.all",penden.env)%*%(get("ck.val",penden.env)))
      assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),t(get("tilde.PSI.d.D.knots.start.g.all",penden.env))),penden.env)
      aa <- try(obj <- solve.QP(Dmat=-get("Derv2.pen",penden.env),dvec=get("Derv1.pen",penden.env),Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE))
      if(class(aa)=="try-error") return("fehler")
      if(any(is.na(obj$solution))) return("fehler")
    }

    else {
      assign("val",get("tilde.PSI.d.D.knots.start.g",penden.env)%*%get("ck.val",penden.env),penden.env)
      if(get("i",penden.env)>1) {
        x <- ecdf(get("val",penden.env))
        prob.val2 <- max(x(round(0.02*max(get("val",penden.env)),2)),0.5)
      }
      while(calc) {
        if(get("i",penden.env)>1) assign("ind.help",which(get("val",penden.env)<=quantile(get("val",penden.env),probs=c(prob.val2))),penden.env)
        else assign("ind.help",seq(1,dim(get("tilde.PSI.d.D.knots.start.g",penden.env))[1]),penden.env)
        if(get("add",penden.env)) bvec <- c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g",penden.env)[get("ind.help",penden.env),]%*%(get("ck.val",penden.env)+rep(1e-06,length(get("ck.val",penden.env)))))
        else bvec <- c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g",penden.env)[get("ind.help",penden.env),]%*%(get("ck.val",penden.env)))
        assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),t(get("tilde.PSI.d.D.knots.start.g",penden.env)[get("ind.help",penden.env),])),penden.env)
        
        aa <- try(obj <- solve.QP(Dmat=-get("Derv2.pen",penden.env),dvec=get("Derv1.pen",penden.env),Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE))
        
        if(class(aa)=="try-error") return("fehler")
        if(any(is.na(obj$solution))) {
          prob.val2 <- prob.val2 + 0.25
        }
        else calc <- FALSE
        if(prob.val2>1) calc <- FALSE
      }
      
      if(any(is.na(obj$solution))) {
        if(get("add",penden.env)) bvec.temp <-  c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g.all",penden.env)%*%(get("ck.val",penden.env)+rep(1e-06,length(get("ck.val",penden.env)))))
        else bvec <- c(rep(0,1+p*length(vec)),-get("tilde.PSI.d.D.knots.start.g.all",penden.env)%*%(get("ck.val",penden.env)))
        assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),t(get("tilde.PSI.d.D.knots.start.g.all",penden.env))),penden.env)
        obj2 <- solve.QP(Dmat=-get("Derv2.pen",penden.env),dvec=get("Derv1.pen",penden.env),Amat=get("Amat",penden.env),bvec=bvec.temp,meq=meq,factorized=FALSE)
        if(any(is.na(obj2$solution))) return("fehler")
        else obj <- obj2
      }
    }
  }

  assign("delta",obj$solution,penden.env)
  assign("ck.val.temp",get("ck.val",penden.env)+obj$solution,penden.env)
  f.hat.val(penden.env,temp=TRUE)
  if(get("no",penden.env)) return("fehler")
  pen.log.like(penden.env,temp=TRUE)
}
