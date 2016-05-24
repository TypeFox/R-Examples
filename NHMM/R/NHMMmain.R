################################################################
## Copyright 2014 Tracy Holsclaw.

## This file is part of NHMM.

## NHMM is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or any later version.

## NHMM is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## NHMM.  If not, see <http://www.gnu.org/licenses/>.
#############################################################



NHMMmain=function(Rgettheta, z, theta, y, yboo, subseqy, subboo, X, betapriorm, betapriorp, K, iters, emdist, burnin, nmix,  W, psipriorm, psipriorp, priors, outdir, outboo,  delta, yrep,  Xp, Wp, ypred, yhold,  ymiss)
{

T=dim(y)[1]
J=dim(y)[2]
L=dim(X)[1]+K
if(!is.null(W)){A=dim(W)[1]}else{A=0}
B=dim(X)[1]
LL=A+K
if(ypred>0 || !is.null(yhold)){pT=dim(Xp)[2]}


Q=iters+burnin
pb = txtProgressBar(min = 0, max = Q, style = 3)
if(delta==FALSE){  mixes=nmix;  delt=0}
if(delta==TRUE) {  mixes=nmix+1; delt=1}
if(emdist=="gamma"){fam=1}
if(emdist=="normal"){fam=2}
if(emdist=="poisson"){fam=3}


#### BICp
BICp=(K-1)*(K+B)      #(K-1)*(K+B) betas (required)
if(A>0){BICp=BICp+(A+K)*J+(K-2)*J}    #psi, gamma (have W)
if(A==0){BICp=BICp+(mixes-1)*K*J}   #ppp (no W)

if(fam==1 || fam==2){ BICp=BICp+K*nmix*J+K*nmix*J }   ## gamma or normal
if(fam==3){ BICp=BICp+K*nmix*J} ## poisson





### alocate and intialize
denzity=matrix(0,K,T)
QQ=array(1/K,dim=c(K,K,T))
beta=matrix(0,K,L)
psi=matrix(0,K+A,J)  #W - A,T,J
mus=matrix(0,T,J)
M=matrix(0,T,J)
gams=c(-Inf,seq(0,1,length=mixes-1),Inf)
gamy=matrix(rep(gams,each=J),mixes+1,J, byrow=TRUE)
gamboo=c(0,rep(1,mixes-1),0)



#########  Set up input variables W for the emissions 
zbin=matrix(0,T,K)
for(k in 1:K) {  zbin[z==k,k]=1 }
Wbin=array(0,dim=c(K+A,T,J))
if(A>0){ Wbin[(K+1):(K+A),,]=W }
Wbin[1:K,,]=getWbin(z,K,J)

ppp= array(0,dim=c(T,J,mixes))
for(j in 1:J){  mus[,j]=matrix(psi[,j],1,LL) %*% Wbin[,,j] }  ### need for rcpp_getppp
ppp=array(rcpp_getppp(gamy, mus),dim=c(T,J,mixes))

#################### set up input variables X for the transitions 
XX=matrix(1,T,L)
XX[1,1:K]=c(1,rep(0,K-1))      ## z[0] is set to state 1
XX[2:T,1:K]=zbin[1:(T-1),1:K]  #beta0  
XX[,(K+1):L]=t(X)         #L inputs


##############  Set up the mixture component latent variable vvv
vvv=matrix(0,T,J)

if(delta==TRUE)  ### assign vvv to non-zero y's
{  v1=rep(1:(mixes-1),each=sum(y>0)/(mixes-1))
   v1=c(rep(1,sum(y>0)-length(v1)),v1)  #makes sum(y>0) and v1 the same length
   v2=v1
   v2[order(y[y>0])]=v1
   vvv[y>0]=v2
}
if(delta==FALSE)
{  v1=rep(0:(mixes-1),each=length(y)/mixes)
   v1=c(rep(0,length(y)-length(v1)),v1)  #makes sum(y>0) and v1 the same length
   v2=v1
   v2[order(y)]=v1
   vvv=matrix(v2,T,J)
}  



############# Saves   #############################
meanQQ=array(0,dim=c(K,K,T))
loglik=numeric(iters)

if(outboo==FALSE)   #beta, theta, z get save to output directory
{  zsave=matrix(0,T,iters)
   betasave=array(0,dim=c(K,L,iters))
   thetasave=array(0,dim=c(2,nmix,K,J,iters)) 
   psisave=array(0,dim=c(K+A,J,iters))
}else{  
   thetasave=NULL
}


yfull=y
ccc=1
ddd=1


if(ypred>0 || !is.null(yhold))
{  zp=rep(1,pT)
   ypred1=matrix(0,pT,J)
   musp=matrix(0,pT,J)
   
   zbinp=matrix(0,pT,K)
   for(k in 1:K) {  zbinp[zp==k,k]=1 }
   Wbinp=array(0,dim=c(K+A,pT,J))
   if(A>0){ Wbinp[(K+1):(K+A),,]=Wp }
   Wbinp[1:K,,]=getWbin(zp,K,J)
   
   pppp= array(0,dim=c(pT,J,mixes))
   for(j in 1:J){  musp[,j]=matrix(psi[,j],1,LL) %*% Wbinp[,,j] }  ### need for rcpp_getppp
   pppp=array(rcpp_getppp(gamy, musp),dim=c(pT,J,mixes))
   
   #################### set up input variables X for the transitions 
   XXp=matrix(1,pT,L)
   XXp[1,1:K]=c(1,rep(0,K-1))      ## z[0] is set to state 1
   XXp[2:pT,1:K]=zbinp[1:(pT-1),1:K]  #beta0  
   XXp[,(K+1):L]=t(Xp)         #L inputs
   
   QQp=array(1/K,dim=c(K,K,pT))
   pls1=matrix(0,pT,J)
   pls2=numeric(pT)
   pls3=numeric(pT)
}


countJ=apply(yboo,2,sum)  #how many missing for each sequence


theK=K

#### MCMC iterations
for(q in 1:Q) 
{  
   if(theK>1){ beta=Rgetbeta(zbin,beta,XX, betapriorm, betapriorp) } #NHMM,  slow but optimized
    
    vvv=rcpp_getvvv(fam, K, mixes, delt, y,c(ppp), c(theta[1,,,]),  c(theta[2,,,]),  z)   
   
   theta=Rgettheta(y,z, priors, theta, nmix, vvv, delt)
   
   
   M=RgetM(A,K, psi, gamy, Wbin, gamboo, vvv )
    if(nmix>1){ gamy=Rgetgams(gamy, vvv,M,nmix,mixes) }
    psi=Rgetpsi(M, psi, Wbin, psipriorm, psipriorp, A, K)
               
           #QQ=array(rcpp_getQQ(K, z, dirprior,  subseqy) ,dim=c(K,K,T))
    if(theK>1){QQ=array(rcpp_getNQQ(beta, XX),dim=c(K,K,T))}
                
    denzity=array(rcpp_getdenzity(A, c(Wbin), psi, gamy, fam, K, mixes, delt, y,c(ppp), c(theta[1,,,]),  c(theta[2,,,])) ,dim=c(K,T,J))
          
          z=Cgetz( z, QQ, denzity, subseqy) #slow 
    #rcpp_getz( z, c(QQ), denzity, subseqy)  #BROKEN 
    for(j in 1:J){  mus[,j]=matrix(psi[,j],1,LL) %*% Wbin[,,j]  }  ### need for rcpp_getppp
    ppp=array(rcpp_getppp(gamy, mus),dim=c(T,J,mixes))
    
    yfull=rcpp_getymiss(fam, K, z, c(ppp), c(theta[1,,,]),  c(theta[2,,,]),mixes, delt, J)
    y[yboo]=yfull[yboo]
    
   
   
    zbin=Cgetzbin(K,z) 
    Wbin[1:K,,]=array(rcpp_getWbin(z,K,J),dim=c(K,T,J))
    XX=CresetX(XX,zbin)   
   
   
   ##print replicates to file
   if(q>(Q-yrep))
   { write(t(yfull), ncolumns=J, paste(outdir,ccc,"-yrep.txt",sep="")) 
     ccc=ccc+1
   }   

   
   if(q>(Q-ypred) || (!is.null(yhold) && q>burnin))
   {
     zbinp=Cgetzbin(K,zp)   ##needs to be done one at a time
     XXp=CresetX(XXp,zbinp) 
     if(theK>1){QQp[,,1:2]=array(rcpp_getNQQ(beta, XXp[1:2,]),dim=c(K,K,2))}  ### really only need QQ[,,1] but dimensionality needs 2     
     zp[1]= rcpp_rmultinom(QQp[z[T],,1])
     
     for(tt in 2:pT)
     {  zbinp=Cgetzbin(K,zp)   ##needs to be done one at a time
        XXp=CresetX(XXp,zbinp) 
        if(theK>1){QQp[,,(tt-1):tt]=array(rcpp_getNQQ(beta, XXp[(tt-1):tt,]),dim=c(K,K,2))} ### really only want tt
        zp[tt]=rcpp_rmultinom(QQp[zp[tt-1],,tt])
     }
     
     Wbinp[1:K,,]=array(rcpp_getWbin(zp,K,J),dim=c(K,pT,J))
     for(j in 1:J){  musp[,j]=matrix(psi[,j],1,LL) %*% Wbinp[,,j]  }  ### need for rcpp_getppp
     pppp=array(rcpp_getppp(gamy, musp),dim=c(pT,J,mixes))
     
     if(q>(Q-ypred) )  #only need to output ypred of these
     {  ypred1=rcpp_getymiss(fam, K, zp, c(pppp), c(theta[1,,,]),  c(theta[2,,,]),mixes, delt, J)
        write(zp, ncolumns=pT, paste(outdir,"zpred.txt",sep=""), append=TRUE) 
        write(t(ypred1), ncolumns=J, paste(outdir,ddd,"-ypred.txt",sep="")) 
        ddd=ddd+1
     }
     
     
     if(!is.null(yhold))  #PLS
     {  denzityp = array(rcpp_getdenzity(A, c(Wbinp), psi, gamy, fam, K, mixes, delt, yhold, c(pppp), c(theta[1,,,]),  c(theta[2,,,])) ,dim=c(K,pT,J))
        for(tt in 1:pT)
        {  pls1[tt,]=denzityp[zp[tt],tt,]  #T by J
        }
        pls1[is.na(yhold)]=mean(pls1[!is.na(yhold)])  #fill in missing values with mean PLS value
        pls2=apply(pls1,1,prod)  #prod for all J
        pls3=pls3+pls2/(Q-burnin) #mean of Q  
      }
   }
   
   

   
  if(q > burnin)  #save iters worth: Q=iters+burnin
  {  
    
     if(ymiss==TRUE)
     {  for(j in 1:J)
        {  if(countJ[j]>0) {write(y[yboo[,j],j], ncolumns=countJ[j] ,paste(outdir,"ymiss-J",j,".txt", sep=""), append=TRUE)}
        }
     }   
     loglik[q-burnin]=Rgetloglik(denzity,z) 
     meanQQ=meanQQ+QQ/iters
     
     if(outboo==FALSE)  
     {   zsave[,q-burnin]=z
         thetasave[,,,,q-burnin]=theta
         betasave[,,q-burnin]=beta
         psisave[,,q-burnin]=psi
     }else{
        for(k in 1:K)  #2,nmix,K,J 
        {  for(j in 1:J)
           {  write(theta[1,,k,j],ncolumns=nmix,paste(outdir,"K",k,"-J",j,"-theta1.txt", sep=""), append=TRUE)
              write(theta[2,,k,j],ncolumns=nmix,paste(outdir,"K",k,"-J",j,"-theta2.txt", sep=""), append=TRUE)
           }
        }   
        
        for(j in 1:J)
        { write(psi[,j],ncolumns=K+A, paste(outdir,"J-",j,"-psi.txt", sep=""), append=TRUE)  #K+A,J
        }
        for(k in 1:K)
        {  write(beta[k,],ncolumns=L,paste(outdir,"K-",k,"-beta.txt", sep=""), append=TRUE)  #K,L
        }
        write(z, ncolumns=T ,paste(outdir,"zsave.txt", sep=""), append=TRUE)
     } 
  }

  theK=length(unique(z))

  
  setTxtProgressBar(pb, q)
} 
close(pb)
    
    PLS=NULL
    if(!is.null(yhold))  #PLS    
    {  PLS=mean(na.omit(log(pls3)))
    }

   if(outboo==FALSE)
   {  fin=list(fam,T,J,K,B,A,iters, burnin, yrep, outboo,outdir, loglik, BICp, meanQQ, betasave, thetasave, psisave, zsave, PLS)
      names(fin)=c("fam","T","J","K","B","A","iters", "burnin","yrep", "outboo","outdir", "loglik", "BICp", "meanQQ", "betasave", "thetasave", "psisave", "zsave","PLS")
   }else{  fin=list(fam,T,J,K,B,A,iters, burnin,yrep,outboo,outdir,loglik, BICp, meanQQ, PLS)
      names(fin)=c("fam","T","J","K","B","A","iters", "burnin","yrep", "outboo","outdir", "loglik", "BICp", "meanQQ", "PLS")
   }
   fin
}


