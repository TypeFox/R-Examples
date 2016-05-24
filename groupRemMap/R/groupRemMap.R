group.remmap <-
function(X, Y,G,lam1,lam2,gamma=0.5, phi0=NULL, C.m=NULL)
{


  ## X:  n (sample size) by p (number of predictors) data matrix; 
  ## Y:  n (sample size) by q (number of responses) data matrix; 
  ## G: group indicator corresponding to each predictor;
  ## lam1: numeric scale; specifying l1 penalty parameter; 
  ## lam2: numeric scale; specifying bridge penalty parameter;
  ## gamma: numeric scale: specifying bridge degree; 
  ## phi0: NULL or numeric matrix (p by q);
  ## C.m:  p (number of predictors) by q (number of responses) data matrix; 


  X=as.matrix(X);
  Y=as.matrix(Y);
  n=nrow(X);
  p=ncol(X);
  q=ncol(Y);

  G_label=unique(G);
  J=length(G_label);
  W=sapply(G_label,function(v){(sum(v==G))^(1-gamma)});

  X.v=as.vector(t(X));
  Y.v=as.vector(t(Y));

  if(is.null(C.m))
  {
    C.v=rep(1, p*q)
  } else
  {
    C.v=as.vector(t(C.m))
  }

  iter.count=0
  RSS=0
  Edebug=rep(0, n*q)

  Phi.m=matrix(0, nrow=p, ncol=q)
  Phi.v=as.vector(Phi.m)
  Phi.vdebug=Phi.v

  
  ###### begin estimation
  if(is.null(phi0))
  {
   junk=.C("groupremmap",
          as.integer(n),
          as.integer(p),
          as.integer(q),
          as.double(X.v),
          as.double(Y.v),
          as.integer(C.v),
          as.integer(G),
          as.integer(J),
          as.integer(G_label),
          as.double(W),
          as.double(lam1[1]),
          as.double(lam2[1]),
          as.double(gamma),
          phi.out=as.double(Phi.v),
          phi.debug=as.double(Phi.vdebug),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)          )  
  } 

 else 
 {
       phi.ini.v=as.vector(t(phi0))  

       junk=.C("groupRemMapIni",
          as.integer(n),
          as.integer(p),
          as.integer(q),
          as.double(X.v),
          as.double(Y.v),
          as.integer(C.v),
          as.integer(G),
          as.integer(J),
          as.integer(G_label),
          as.double(W),
          as.double(lam1[1]),
          as.double(lam2[1]),
          as.double(gamma),
          as.double(phi.ini.v),
          phi.out=as.double(Phi.v),
          phi.debug=as.double(Phi.vdebug),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)          )  
  }


  phi.result=matrix(junk$phi.out, nrow=p, byrow=T)
  E=matrix(junk$Edebug, nrow=n, ncol=q, byrow=T)
  rss.v=apply(E^2,2,sum)
  iter.count=junk$n.iter

  return(list(phi=phi.result, rss.v=rss.v))
}

