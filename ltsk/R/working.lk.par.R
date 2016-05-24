working.lk.par <- function(l.query,l.obs,th,xcoord='x',ycoord='y',zcoord='z',vlen=15)
{
  th <- th^2 ## sqaured the threshold to accomodate ANN search
  l.vlen = as.integer(vlen)
  
  ## invoke the local kriging library
  nn <- nrow(l.query)
  out <- .C("lk_main",
	as.double(l.obs[,xcoord]),
	as.integer(nrow(l.obs)),
	as.double(l.obs[,ycoord]),
	as.integer(nrow(l.obs)),
	as.double(l.obs[,zcoord]),
	as.integer(nrow(l.obs)),
	as.double(l.query[,xcoord]),
	as.integer(nrow(l.query)),
	as.double(l.query[,ycoord]),
	as.integer(nrow(l.query)),
	as.double(th),
	as.integer(l.vlen),
	krig=double(nn),
	sigma=double(nn),
	hs=double(nn),
	psill=double(nn),
	nugget=double(nn),
	ms=integer(nn))
  ## collect the result
  l.out <- cbind(1:nn,out$krig,out$sigma,out$hs,out$nugget,out$psill,out$ms)
  colnames(l.out) <- c('null','krig','sigma','Hs','nugget','psill','model')
  ## recode missing values 
  ii <- which(l.out[,2]==-99999)
  jj <- which(l.out[,3]==-99999)
  if(length(ii)>0) l.out[ii,2] <- NA
  if(length(jj)>0) l.out[jj,3] <- NA
  data.frame(l.out[,-1])
}
