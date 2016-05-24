
smpl<- function(pr,x){
   n<- nrow(pr)
   y<- rep(NA,n)
   for(i in 1:n) y[i]<- sample(x,size=1,prob=pr[i,])
   y
}

# estimate variances of (genetic) components
gvar<- function (par, ov){
# par: parameter estimates
# ov: ov$v is a list of genetic matrices
   str<- c("AA", "DD", "HH", "AD", "MH", "EE")
   if(!is.null(par)) par<- par[str]
   ov<- ov$v[str]
   vv <- rep(NA, 6)
      names(vv)<- str
   for (i in 1:6){
      v0<- ov[[i]]*par[i]
         v0<- diag(v0)
      vv[i]<- mean(v0)
   }
   vv
}

# estimate QTL variances
qtlVar<- function(lrt,prdat,simulation=FALSE,nsim=25){
# lrt: data frame (a,d,...)
# probs: prDat$pr
# estimated genetic variance-covariance matrix
   vv<- rep(NA,nrow(lrt))
   for(ii in 1:nrow(lrt)){
      tmp<- lrt[ii,]
      prd<- prdat[,,ii]

      if(!simulation){
         tmp1<- sweep(prd,2,c(tmp$a,tmp$d,-tmp$a),"*")
            tmp1<- rowSums(tmp1) # mean
         tmp2<- cbind(tmp$a-tmp1,tmp$d-tmp1,-tmp$a-tmp1)
            tmp2<- tmp2^2
            tmp2<- tmp2*prd
            tmp2<- rowSums(tmp2) # variance

         vv[ii]<- var(tmp1) + mean(tmp2)
      }else{# simulation method -- takes time
         vr<- rep(NA,nsim)
         for(i in 1:nsim){
            vr.<- rep(NA,nrow(prd))
            pr.<- runif(nrow(prd), min=0, max=1)
            idx1<- pr. <= prd[,1]
            idx2<- (!idx1) & (pr. <= prd[,1]+prd[,2])
            idx3<- !idx1 & !idx2
            if(any(idx1)) vr.[idx1]<- tmp$a
            if(any(idx2)) vr.[idx2]<- tmp$d
            if(any(idx3)) vr.[idx3]<- -tmp$a
            vr[i]<- var(vr.)
         }
         vv[ii]<- mean(vr)
      }
   }
   vv
}

