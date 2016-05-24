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



NHMM_MVNmain=function( z, theta, y, yboo, subseqy, subboo, X, betapriorm, betapriorp, K, iters,  burnin,   W, psipriorm, psipriorp,  priors1, priors2, outdir,outboo,  yrep, Xp, Wp, ypred, yhold, ymiss)
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


#### BICp
BICp=(K-1)*(K+B)      #(K-1)*(K+B) betas (required)
if(A>0){BICp=BICp+(A+K)*J}    #psi, gamma (have W)
if(A==0){BICp=BICp+K*J}   #ppp (no W)

 BICp=BICp+(K*J*J)/2+(.5*J)*K   ## mvn Sigma (symetric)






### alocate and intialize
denzity=matrix(0,K,T)
QQ=array(1/K,dim=c(K,K,T))
beta=matrix(0,K,L)
psi=matrix(0,K+A,J)  #W - A,T,J

for(j in 1:J)
{  for(k in 1:K)
   {  psi[k,j]=mean(y[z==k,j])
   }
}

mus=matrix(0,T,J)



#########  Set up input variables W for the emissions 
zbin=matrix(0,T,K)
for(k in 1:K) {  zbin[z==k,k]=1 }
Wbin=array(0,dim=c(K+A,T,J))
if(A>0){ Wbin[(K+1):(K+A),,]=W }
Wbin[1:K,,]=getWbin(z,K,J)


for(j in 1:J){  mus[,j]=matrix(psi[,j],1,LL) %*% Wbin[,,j] }  ### need for rcpp_getppp


#################### set up input variables X for the transitions 
XX=matrix(1,T,L)
XX[1,1:K]=c(1,rep(0,K-1))      ## z[0] is set to state 1
XX[2:T,1:K]=zbin[1:(T-1),1:K]  #beta0  
XX[,(K+1):L]=t(X)         #L inputs

vvv=matrix(1,T,J)


############# Saves   #############################
meanQQ=array(0,dim=c(K,K,T))
loglik=numeric(iters)

if(outboo==FALSE)   #beta, theta, z get save to output directory
{  zsave=matrix(0,T,iters)
   betasave=array(0,dim=c(K,L,iters))
   thetasave=array(0,dim=c(J,J,K)) 
   psisave=array(0,dim=c(K+A,J,iters))
}else{  
  thetasave=array(0,dim=c(J,J,K)) 
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
   
   #################### set up input variables X for the transitions 
   XXp=matrix(1,pT,L)
   XXp[1,1:K]=c(1,rep(0,K-1))      ## z[0] is set to state 1
   XXp[2:pT,1:K]=zbinp[1:(pT-1),1:K]  #beta0  
   XXp[,(K+1):L]=t(Xp)         #L inputs
   
   QQp=array(1/K,dim=c(K,K,pT))
   pls1=numeric(pT)
   pls2=numeric(pT)
   miss=apply(is.na(yhold),1,sum)>0 #find missing in yhold
   
}


countJ=apply(yboo,2,sum)  #how many missing for each sequence




theK=K



#theta=array(0,dim=c(J,J,K))
thetainv=array(0,dim=c(J,J,K))
cholS=array(0,dim=c(J,J,K))
detS=numeric(K)
denzitypass=array(1,dim=c(K,T,J))


### TRIAL SET UP
#vars=numeric(J)
#for(j in 1:J){vars[j]=var(y[,j])}
#for(i in 1:K){theta[,,i]=diag(J)*4}


#### MCMC iterations ################################################
for(q in 1:Q) 
{  
   if(theK>1){ beta=Rgetbeta(zbin,beta,XX, betapriorm, betapriorp) } #NHMM,  slow but optimized
    
   for(j in 1:J){  mus[,j]=matrix(psi[,j],1,LL) %*% Wbin[,,j]  } 
    theta=RgetthetaMVN(y,z, priors1, priors2, theta, mus) #Sigma
   
    for(k in 1:K)
    {  thetainv[,,k]=solve(theta[,,k])
       cholS[,,k]=t(chol(theta[,,k]))
       detS[k]=det(thetainv[,,k])
    }
   
    psi=RgetpsiMVN(y, z, psi, Wbin, psipriorm, psipriorp, A, K, thetainv)
               
    if(theK>1){QQ=array(rcpp_getNQQ(zbin, beta),dim=c(K,K,T))}
    
   
    denzity=rcpp_getdenzityMVN(A, c(Wbin), psi, K, y,c(thetainv),detS)
    denzitypass[,,1]=denzity
          z=Cgetz( z, QQ, denzitypass, subseqy) #slow 
 
   
   for(j in 1:J){  mus[,j]=matrix(psi[,j],1,LL) %*% Wbin[,,j]  } 
   for(t in 1:T)
   {  yfull[t,]=mvrnorm(1,mus[t,],theta[,,z[t]])
   }
    #yfull=rcpp_getymissMVN(fam, K, z, c(ppp), c(theta[1,,,]),  c(theta[2,,,]),mixes, delt, J)
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
     zbinp=Cgetzbin(K,zp)
     XXp=CresetX(XXp,zbinp) 
     if(theK>1){QQp=array(rcpp_getNQQ(zbinp, beta),dim=c(K,K,pT))}
     
     zp[1]= rcpp_rmultinom(QQp[z[T],,1])
     for(tt in 2:pT)
     {  zp[tt]=rcpp_rmultinom(QQp[zp[tt-1],,tt-1])
     }
     Wbinp[1:K,,]=array(rcpp_getWbin(zp,K,J),dim=c(K,pT,J))
     for(j in 1:J){  musp[,j]=matrix(psi[,j],1,LL) %*% Wbinp[,,j]  }  ### need for rcpp_getppp

     if(q>(Q-ypred) )  #only need to output ypred of these
     {  for(t in 1:pT)
        {  ypred1[t,]=mvrnorm(1,musp[t,],theta[,,zp[t]])
        }
        #ypred1=rcpp_getymissMVN(fam, K, zp, c(pppp), c(theta[1,,,]),  c(theta[2,,,]),mixes, delt, J)
        write(zp, ncolumns=pT, paste(outdir,"zpred.txt",sep=""), append=TRUE) 
        write(t(ypred1), ncolumns=J, paste(outdir,ddd,"-ypred.txt",sep="")) 
        ddd=ddd+1
     }
     
     
     if(!is.null(yhold))  #PLS     K by T
     {  denzityp = rcpp_getdenzityMVN(A, c(Wbinp), psi, K, yhold, c(thetainv),detS)
        for(tt in 1:pT)
        {  pls1[tt]=denzityp[zp[tt],tt]
        }
        pls1[miss]=mean(pls1[!miss])   #fill in missing values with mean PLS value
        pls2=pls2+pls1/(Q-burnin) 
     }
     

   }
   
   
   
  if(q > burnin)  #save iters worth: Q=iters+burnin
  {  
    if(ymiss==TRUE)
    {  for(j in 1:J)
       {  if(countJ[j]>0) {write(y[yboo[,j],j], ncolumns=countJ[j] ,paste(outdir,"ymiss-J",j,".txt", sep=""), append=TRUE)}
       }
    } 
    
     loglik[q-burnin]=Rgetloglik(denzitypass,z) 
     meanQQ=meanQQ+QQ/iters
     
     if(outboo==FALSE)  
     {   zsave[,q-burnin]=z
         thetasave[,,]=theta/(q-burnin)
         betasave[,,q-burnin]=beta
         psisave[,,q-burnin]=psi
     }else{
       thetasave[,,]=theta/(q-burnin)
       # for(k in 1:K)  #2,nmix,K,J 
       # {    write(c(theta[,,k]),ncolumns=J*J,paste(outdir,"K",k,"-theta.txt", sep=""), append=TRUE)
       # }   
        
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
   {  PLS=mean(na.omit(log(pls2)))
   }

if(outboo==TRUE)
{   for(k in 1:K)  #2,nmix,K,J 
    {    write(t(thetasave[,,k]),ncolumns=J,paste(outdir,"K",k,"-theta.txt", sep=""))
    }
}

fam=4  #MVN
   if(outboo==FALSE)
   {  fin=list(fam,T,J,K,B,A,iters, burnin, yrep, outboo,outdir, loglik, BICp, meanQQ, betasave, thetasave, psisave, zsave, PLS)
      names(fin)=c("fam","T","J","K","B","A","iters", "burnin","yrep", "outboo","outdir", "loglik", "BICp", "meanQQ", "betasave", "thetasave", "psisave", "zsave", "PLS")
   }else{  fin=list(fam,T,J,K,B,A,iters, burnin,yrep,outboo,outdir,loglik, BICp, meanQQ, PLS)
      names(fin)=c("fam","T","J","K","B","A","iters", "burnin","yrep", "outboo","outdir", "loglik", "BICp", "meanQQ", "PLS")
   }
   fin
}


