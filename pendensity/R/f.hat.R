f.hat <- function(penden.env,ck.temp=NULL) {
  l.y <- matrix(1:get("n",penden.env))
  if(is.null(ck.temp)) ck.temp <- get("ck.temp",penden.env)
  Stand.abw <- get("Stand.abw",penden.env)
  base <- get("base",penden.env)
  Z.index <- get("Z.index",penden.env)

  if(base=="gaussian") {
    if(get("x.null",penden.env)) return(colSums(apply(l.y,1,function(i,ck.temp,Stand.abw,base.den,Z.index) (ck.temp/Stand.abw)*base.den[,i],ck.temp,Stand.abw,get("base.den",penden.env),Z.index)))
    else return(colSums(apply(l.y,1,function(i,ck.temp,Stand.abw,base.den,Z.index) (ck.temp/Stand.abw)[Z.index[i],]*base.den[,i],ck.temp,Stand.abw,get("base.den",penden.env),Z.index)))
  }

  if(base=="bspline") {
    if(get("x.null",penden.env)) return(colSums(apply(l.y,1,function(i,ck.temp,base.den) ck.temp*base.den[,i] ,ck.temp,get("base.den",penden.env))))
    else return(colSums(apply(l.y,1,function(i,ck.temp,base.den,Z.index) ck.temp[Z.index[i],]*base.den[,i] ,ck.temp,get("base.den",penden.env),Z.index)))
  }
} 
