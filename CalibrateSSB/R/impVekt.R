

ImpVekt = function(mBrutto,r,
            totPop=dim(mBrutto)[1],N_brutto=TRUE,totPopReturn=FALSE,singularReturn=FALSE,ginvtol=1e-06,w=NULL)
{
  wGross = NULL
  if(singularReturn) totPopReturn=TRUE
       inv = function(x) my_ginv(x,tol=ginvtol)

  if(is.null(w)) w = 1
  sqrt_w = sqrt(w)


  totPop_input = totPop
  totPop = rep(NaN,dim(mBrutto)[2])
  for(i in 1:length(totPop_input))
    totPop[i] = totPop_input[i]
  totPop_input = totPop

  tp = !is.na(totPop)

  wregw = function(tot,x,sqrt_w)
  {
     (tot %*% inv(x*sqrt_w)*sqrt_w)[,]
  }


  if(sum(!tp)>0)
     wGross = wregw(totPop[tp],mBrutto[,tp,drop=FALSE],sqrt_w)
     totPop[!tp] = colSums(mBrutto[,!tp,drop=FALSE]*wGross)

  if(length(sqrt_w)>1)	sqrt_w=sqrt_w[r==1]
  vekt = wregw(totPop,mBrutto[r==1,,drop=FALSE],sqrt_w)


  if(N_brutto)
  {
    vekt_netto = vekt
    vekt = r
    vekt[r==1] = vekt_netto
  }


  vekt = as.vector(unlist(vekt,use.names=FALSE))
  if(totPopReturn)
  {
    totP        = totPop
	  totPop      = totPop_input
    if(singularReturn)
      return(list(vekt=vekt,totPop=totPop,wGross=wGross))
      #singular1=svd(R,nv=0,nu=0)$d, # R gir samme svar som mBrutto
      #singular2=svd(qtw %*%q,nv=0,nu=0)$d))
    return(list(vekt=vekt,totPop=totPop,wGross=wGross))
  }
  vekt
}





glmR2ImpVekt = function(glmR,...)
{
   ImpVekt(model.matrix(glmR),glmR$y==1,...)
}


ImpVektFixed = function(mNetto,totPop,wFixed=NULL,useginv=TRUE,ginvtol=1e-06)
{

  if(useginv)
       inv = function(x) my_ginv(x,tol=ginvtol)
  else inv = function(x) solve(x)

  qr_mNetto = qr(mNetto) # mNetto = Q*R
  Q = qr.Q(qr_mNetto)
  R = qr.R(qr_mNetto)
  invR = inv(R)

  totPop_invR = totPop %*% invR

  vekt = as.vector(totPop_invR %*% t(Q))

  if(is.null(wFixed)) return(vekt)

  z = matrix(wFixed - vekt,length(vekt),1)


  za = z[!is.na(wFixed),,drop=FALSE]
  Qa = Q[!is.na(wFixed),,drop=FALSE]
  Qb = Q[is.na(wFixed),,drop=FALSE]


  d =  - t(Qa)  %*% za


  zb = (Qb %*% inv(t(Qb)%*%Qb)) %*% d


  z[is.na(z)] = zb

  nyVekt = as.vector(z) + vekt

  return(nyVekt)
}

