
rblasso <-
function(s, m, om, nlam, tol=1e-5, sbols=1, maxit, quiet=0, warm=0, B0=NULL)
{
  p=dim(s)[1]
  q=dim(om)[1]
  s = as.double(s)
	m = as.double(m)
	om = as.double(om)
	nlam=as.double(nlam)
  tmp=matrix(0, nrow=p, ncol=q)
	tol=as.double(tol*sbols)

  mode(p) = "integer"
  mode(q) = "integer"
  mode(maxit) = "integer"
	mode(quiet) = "integer"
	mode(warm) = "integer"

	if(!warm) tmp=as.double(tmp) else tmp = as.double(B0)

	dotCoutput=.C("blasso",Sin=s, Min=m, Omin=om,pin=p,
            qin=q, lamin=nlam, tol=tol,maxit=maxit, Bout=tmp, warm=warm)

  B = matrix(dotCoutput$Bout, nrow=p, ncol=q)
  return(B)
}
