eqlipse<-function(x,y , cov,   wcols = c(1,2) , dof=2, pct=0.05, ... )
{
  ######   add an ellipse to the
  ###  X-Y plot of an earthquake
###   cov = wup$cov[ 2:4, 2:4]

  if(missing(dof) ){ dof = 2 }
   if(missing(wcols) ){  wcols = c(1,2) }
  if(missing( pct) ){  pct = 0.05 }
  if(is.null(dof)) { dof=2 }


  
  theta=seq(from=0, by=.01, to=2*pi);

  r=matrix( rep(0, times=2*length(theta)), ncol=2)

  C=cov[wcols,wcols]

  circ1 = cbind(cos(theta), sin(theta))

  U2 = eigen(solve(C))

####   flip: this is the bizarre matlab way of doing things
  u = U2$vectors[, 2:1]
  lam = rev(U2$values)
###########
  
###  delta=sqrt(qchisq(0.95,dof));

########  95 percent confidence bounds



  pfact = 1-pct/2
  
  delta= qt(pfact, dof)
  
  zed =  diag( (delta/  sqrt(lam)) )
  rmore =  circ1 %*% zed  %*%  t(u)

  ##   plot(m[1]+rmore[,1],m[3]+rmore[,2],type='l', xlim=c(-50, 50) , ylim=c(7, 12) , ann=FALSE);

   ##   if(is.null(border)){ bcol="black" }
  
  polygon(x+rmore[,1],y+rmore[,2], ... );

}

