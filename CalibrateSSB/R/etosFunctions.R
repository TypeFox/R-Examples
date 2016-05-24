
# (5.3.14) in ETOS manual (without "- -" typing error)
# (11.8) in Sarndal-Lundstrom-2005
# Without summing over strata
etos_5.3.14 = function(N,n,g,e)
{
  ge=g*e
  (N^2*(1-n/N)/(n*(n-1))) * (sum(ge^2) - (sum(ge))^2/n ) - (N/n)*(N/n-1)*sum(ge*(g-1)*e)
}


# (5.3.15) in ETOS manual
# (11.9) in Sarndal-Lundstrom-2005
# Without summing over strata
etos_5.3.15 = function(N,n,g,e)
{
  (N/n)^2*sum(g*(g-1)*e^2)
}


# Beregner residualer fra modelmatrise og matrise med y-er
# isTotPop er vektor som sier hvilke variabler det finnes tot-pop for
#          Det legges til FALSE om vektoren er for kort
#          Default er at bare N er kjent i populasjonen
# w er designvekter
etos_e1_e2 = function(mNetto,yMatrix,isTotPop=TRUE,ginvtol=1e-06,w=NULL)
{
  inv = function(x) my_ginv(x,tol=ginvtol)
  if(is.null(w)) w = 1
  sqrt_w = sqrt(w)

  isTotPop_input = isTotPop
  isTotPop = rep(FALSE,dim(mNetto)[2])
  for(i in 1:length(isTotPop_input))
    isTotPop[i] = isTotPop_input[i]

  reswreg = function(y,x,sqrt_w)
    y-(x %*% inv(x*sqrt_w)*sqrt_w) %*% y
  a=NULL
  a$e1=reswreg(yMatrix,mNetto[,isTotPop,drop=FALSE],sqrt(w))
  if(sum(!isTotPop)==0) a$e2=a$e1
  else a$e2 = reswreg(yMatrix,mNetto,sqrt(w))
  a
}


etos_e1_e2_by_lm = function(mNetto,yMatrix,isTotPop=TRUE,lmInfluence=TRUE,w=NULL)
{

  isTotPop_input = isTotPop
  isTotPop = rep(FALSE,dim(mNetto)[2])
  for(i in 1:length(isTotPop_input))
    isTotPop[i] = isTotPop_input[i]

  a=NULL
  m   = lm(yMatrix~mNetto,weights=w)
  a$e1 = resid(m)
  a$e2 = a$e1
  if(lmInfluence){
    a$h1 = lm.influence(m)$hat
    a$h2 = a$h1
  }
  if(!(sum(!isTotPop)==0)){
    mNetto = mNetto[,isTotPop,drop=FALSE]
    m   = lm(yMatrix~mNetto,weights=w)
    a$e1 = resid(m)
    if(lmInfluence) a$h1 = lm.influence(m)$hat
  }
  a
}


etosV1V2 = function(e1,e2,w,samplingWeights=NULL,R=is.finite(e1),  N = sum(w,na.rm=T),n = sum(is.finite(w))){

  if(is.null(samplingWeights))   g=w*(n/N)
  else g = w/samplingWeights

  V1 = etos_5.3.14(N,n,g[R],e1[R])
  V2 = etos_5.3.15(N,n,g[R],e2[R])

  list(V1=V1,V2=V2)
}

