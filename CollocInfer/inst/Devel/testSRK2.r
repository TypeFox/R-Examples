## Assume that FhNEx.R demo has already been run (at least lines 1 through 31
## and 73 through 76

srkproc = make.SRK2proc()

srktimes = seq(0,20,len=81)
#srktimes = seq(0,20,len=801)
#coefs = eval.fd(srktimes, DEfd)
#coefs = DEfd$coefs
more = proc$more
bvals = list(bvals=diag(rep(1,81)),I=SRK2fns()$SRK2indeces(81))
#bvals = list(bvals=diag(rep(1,801)),I=SRK2fns()$SRK2indeces(801))
#bvals = list(bvals=cbind( rbind(diag(40)%x%matrix(c(1,0),2,1),rep(0,40)), c(rep(0,80),1),c(rep(0,80),1)),I=SRK2fns()$SRK2indeces(81))
coefs=matrix(0,81,2)
#coefs=matrix(0,801,2)
coefs[c(bvals$I[,1], length(srktimes)),]=FhNdata
#coefs[c(bvals$I[,1], length(srktimes)),]=data
coefs[bvals$I[,3],]=(FhNdata[1:40,]+FhNdata[2:41,])/2
#coefs[bvals$I[,3],]=(data[1:400,]+data[2:401,])/2


more$weights = matrix(8/3,81,2)
#more$weights = matrix(8/3,801,2)
more$qpts = srktimes

srkproc$bvals = bvals
srkproc$more = more


srklik = lik
srklik$bvals = t(cbind(rbind(diag(40)%x%matrix(c(1,0),2,1),rep(0,40)),c(rep(0,80),1)))
#srklik$bvals = t(cbind(rbind(diag(400)%x%matrix(c(1,0),2,1),rep(0,400)),c(rep(0,800),1)))
#srklik$bvals =  cbind(diag(41),rep(0,41))

# Now try some optimization

pres = ParsMatchOpt(FhNpars,coefs,srkproc)
npars = pres$pars

# Smoothing can be done more specifically with 'inneropt'

Ires2 = inneropt(FhNdata,times=FhNtimes,FhNpars,coefs,srklik,srkproc)
#Ires2 = inneropt(data,times=seq(0,20,len=401),FhNpars,coefs,srklik,srkproc)
# And we can also do profiling

Ores2 = outeropt(FhNdata,FhNtimes,npars,coefs,srklik,srkproc)









### The rest of this is just finite-differencing checks on my derivatives

f = make.SRK2proc()$fn(coefs,bvals,pars,more)
gp = make.SRK2proc()$dfdp(coefs,bvals,pars,more)
gc = make.SRK2proc()$dfdc(coefs,bvals,pars,more)

Hcc = make.SRK2proc()$d2fdc2(coefs,bvals,pars,more)
Hcp = make.SRK2proc()$d2fdcdp(coefs,bvals,pars,more)

dfdplist = SRK2fns()$dSRK2dp(coefs,bvals,pars,more,return.diffs=TRUE)

d2fdxdplist = SRK2fns()$d2SRK2dxdp(coefs,bvals,pars,more)


ed2fdxdp11 = 0*d2fdxdplist$d2fdx11
ed2fdxdp12 = 0*d2fdxdplist$d2fdx12 
ed2fdxdp21 = 0*d2fdxdplist$d2fdx21
ed2fdxdp22 = 0*d2fdxdplist$d2fdx22 
for(i in 1:ncol(coefs)){
  for(j in 1:ncol(coefs)){
    for(k in 1:length(pars)){
     ed2fdxdp11[,,i,j] = d2fdxdplist$dfdx11[,i,j]*d2fdxdplist$dfdp1[,i,k] 
     ed2fdxdp12[,,i,j] = d2fdxdplist$dfdx12[,i,j]*d2fdxdplist$dfdp1[,i,k]
     ed2fdxdp21[,,i,j] = d2fdxdplist$dfdx21[,i,j]*d2fdxdplist$dfdp2[,i,k]
     ed2fdxdp22[,,i,j] = d2fdxdplist$dfdx22[,i,j]*d2fdxdplist$dfdp2[,i,k]
    }
  }
}


eps = 1e-6

agp = 0*gp
aHcp = 0*Hcp

adfdp1 = 0*dfdplist$dfdp1
adfdp2 = 0*dfdplist$dfdp2 

ad2fdxdp11 = 0*d2fdxdplist$d2fdx11
ad2fdxdp12 = 0*d2fdxdplist$d2fdx12 
ad2fdxdp21 = 0*d2fdxdplist$d2fdx21
ad2fdxdp22 = 0*d2fdxdplist$d2fdx22 

for(i in 1:length(pars)){
  tpars = pars
  tpars[i] = pars[i] + eps

  tf =  make.SRK2proc()$fn(coefs,bvals,tpars,more)
  tgc = make.SRK2proc()$dfdc(coefs,bvals,tpars,more)
  
  trl =  SRK2fns()$dSRK2dp(coefs,bvals,tpars,more,return.diffs=TRUE)
  
  adfdp1[,,i] = (trl$diffs1-dfdplist$diffs1)/eps
  adfdp2[,,i] = (trl$diffs2-dfdplist$diffs2)/eps

  trl2 = SRK2fns()$dSRK2dx(coefs,bvals,tpars,more,return.diffs=FALSE)

  ad2fdxdp11[,,,i] = (trl2$dfdx11-d2fdxdplist$dfdx11)/eps
  ad2fdxdp12[,,,i] = (trl2$dfdx12-d2fdxdplist$dfdx12)/eps
  ad2fdxdp21[,,,i] = (trl2$dfdx21-d2fdxdplist$dfdx21)/eps
  ad2fdxdp22[,,,i] = (trl2$dfdx22-d2fdxdplist$dfdx22)/eps

  agp[i] = (tf-f)/eps
  aHcp[,i] = (tgc-gc)/eps
}



agc = 0*gc
aHcc = 0*Hcc
aHcp2 = 0*Hcp

for(i in 1:length(coefs)){
  tcoefs = coefs
  tcoefs[i] = coefs[i] + eps

  tf =  make.SRK2proc()$fn(tcoefs,bvals,pars,more)
  tgc = make.SRK2proc()$dfdc(tcoefs,bvals,pars,more)
  tgp = make.SRK2proc()$dfdp(tcoefs,bvals,pars,more)

  agc[i] = (tf-f)/eps
  aHcc[,i] = (tgc-gc)/eps
  aHcp2[i,] = (tgp-gp)/eps
}



dfdxlist = SRK2fns()$dSRK2dx(coefs,bvals,pars,more,return.diffs=TRUE)
d2fdx2list = SRK2fns()$d2SRK2dx2(coefs,bvals,pars,more)


adfdx11 = 0*dfdxlist$dfdx11
adfdx12 = 0*dfdxlist$dfdx12
adfdx21 = 0*dfdxlist$dfdx21
adfdx22 = 0*dfdxlist$dfdx22

ad2fdx111 = 0*d2fdx2list$d2fdx111
ad2fdx122 = 0*d2fdx2list$d2fdx122
ad2fdx112 = 0*d2fdx2list$d2fdx112
ad2fdx211 = 0*d2fdx2list$d2fdx211

for(i in 1:ncol(coefs)){
    tcoefs1 = coefs
    tcoefs1[bvals$I[,1],i] = tcoefs1[bvals$I[,1],i] + eps
    tcoefs2 = coefs
    tcoefs2[bvals$I[,3],i] = tcoefs2[bvals$I[,3],i] + eps

    trl1 =  SRK2fns()$dSRK2dx(tcoefs1,bvals,pars,more,return.diffs=TRUE)
    trl2 =  SRK2fns()$dSRK2dx(tcoefs2,bvals,pars,more,return.diffs=TRUE)
    
    adfdx11[,,i] = (trl1$diffs1-dfdxlist$diffs1)/eps
    adfdx21[,,i] = (trl1$diffs2-dfdxlist$diffs2)/eps

    adfdx12[,,i] = (trl2$diffs1-dfdxlist$diffs1)/eps
    adfdx22[,,i] = (trl2$diffs2-dfdxlist$diffs2)/eps
    
    ad2fdx111[,,,i] = (trl1$dfdx11-d2fdx2list$dfdx11)/eps
    ad2fdx122[,,,i] = (trl2$dfdx12-d2fdx2list$dfdx12)/eps
    ad2fdx112[,,,i] = (trl2$dfdx11-d2fdx2list$dfdx11)/eps
    ad2fdx211[,,,i] = (trl1$dfdx21-d2fdx2list$dfdx21)/eps
}



