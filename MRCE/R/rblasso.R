rblasso <-
function(s, m, om, nlam, n, B0=NULL, soft=NULL, objective=0, tol=1e-5, maxit=500, quiet=TRUE)
{
  p=dim(s)[1]
  q=dim(om)[1]
  if(is.null(B0))
    B0=as.double(rep(0, p*q))
  if(is.null(soft))
   soft=as.double(m)  
 
  objective=as.double(objective)
  soft=as.double(soft)
  s = as.double(s)
  m = as.double(m)
  om = as.double(om)
  nlam=as.double(nlam)
  tol=as.double(tol)
  totalit=0
  mode(n) = "integer"
  mode(p) = "integer"
  mode(q) = "integer"
  mode(maxit) = "integer"
  mode(totalit)="integer"
  dotCoutput=.C("blasso", B=B0, S=s, M=m, Om=om, soft=soft, pin=p,
            qin=q, nin=n, lam=nlam, tol=tol, maxit=maxit, totalit=totalit, objective=objective)
	
  if(!quiet)
  {  
    cat("Total iterations for solving for B was", dotCoutput$totalit, "\n")
  }
  B = matrix(dotCoutput$B, nrow=p, ncol=q)
  return(list(B=B, obj=dotCoutput$objective))
}

