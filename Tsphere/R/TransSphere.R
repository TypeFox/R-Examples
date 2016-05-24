TransSphere <-
function(dat,y,fdr,minlam,maxlam=NULL)
{
  n = nrow(dat); p = ncol(dat);
  xc = meanTranspose(dat)$xcen
  qm = .8; qn = 0.6618683;
  indy = y==1
  ttin = function(x,ind){return(t.test(x[ind],x[!ind],alternative="two.sided",paired=FALSE,var.equal=TRUE)$statistic)}
  tto = apply(xc,1,ttin,indy)
  #central-matching scale value
  qtto = sd(tto[tto>=quantile(tto,(1-qm)/2) & tto<=quantile(tto,1-(1-qm)/2)])
  psi = qn/qtto
  #taking top 500 genes
  or = order(abs(tto),decreasing=TRUE)
  nn = min(500,n)
  addin = or[1:nn]
  x = xc[addin,] 
  #sphering algorithm
  mu1hat = apply(x[,y==1],1,mean)
  mu2hat = apply(x[,y==2],1,mean)
  Smat = matrix(mu1hat,nn,p)
  Smat[,y==2] = matrix(mu2hat,nn,sum(y==2))
  xc1 = x - Smat
  temp = meanTranspose(xc1)
  xc2 = temp$xcen
  if(!is.null(maxlam)){
    cvs = CVcov(xc2,maxlam,minlam,steps=5,do=1)
    optlam = cvs$optlam
  }else{optlam = minlam[1]}
  ct = covTranspose11(xc2,optlam*nn,optlam*p,trace=TRUE)
  Sigi = ct$Sigmaihat
  Delti = ct$Deltaihat
  xn = chol(Sigi)%*%xc2%*%chol(Delti)
  xsph = xn + Smat
  #finding significnat genes
  ttsph = apply(xsph,1,ttin,indy)*psi
  pvsph = rep(0,nn)
  pvsph[ttsph>=0] = 1 - pt(ttsph[ttsph>=0],n-2)
  pvsph[ttsph<0] = pt(ttsph[ttsph<0],n-2)
  op = order(pvsph)
  temp = pvsph[op]<=(c(1:length(pvsph))*fdr/n)
  nsiggene = which(!temp)[1] - 1
  siggenes = addin[op[1:nsiggene]]
  return(list(sig.rows=siggenes,t.stats=ttsph,p.vals=pvsph,x.sphered=xsph))
}

