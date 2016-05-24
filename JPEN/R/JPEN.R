# trace function
tr=function(A) {return(sum(diag(A)))};
# Compute lambda
lamvec=function (c, gam, p) 
{
  up = 2 * (c + gam);
  arb=c(1/1024,1/512,1/256,1/128,1/64,1/32,1/16,1/8,1/4,1/2,1);
  arb1=up*arb;
  return(arb1)
}
# data splitiing in k-fold
f.K.fold <- function(Nobs,K=5) {   rs <- runif(Nobs);  id <- seq(Nobs)[order(rs)];   k <- as.integer(Nobs*seq(1,K-1)/K);   k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE);  k[,1] <- k[,1]+1;  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id);  return(l); }
# main estimation function.
jpen=function(S,gam,lam=NULL)  
  {
    p=dim(S)[1];c1=max(abs(S));
    if(length(lam)==0) 
    {lam=2*gam/p}; 
    E=(S+gam*diag(p))/(1 + gam); 
      E1 = sign(E) * pmax(abs(E) - lam/(2 * (1 + gam)), 0);   
      diag(E1) = diag(E); if(lam<=2*gam/p) {return(E1)} else {if(min(eigen(E1)$values)>.0001) {return(E1)} else {lam=2*gam/p;  E = (S + gam * diag(p))/(1 + gam);    E1 = sign(E) * pmax(abs(E) - lam/(2 * (1 + gam)), 0);   diag(E1) = diag(E);    return(E1)
      }
    }
  }
# inverse cov matrix estimation
jpen.inv=function(S,gam,lam=NULL) {p=dim(S)[1];S1=S;if(min(eigen(S)$values) <.0001) {c=max(abs(S1[col(S1)!=row(S1)]));gstr=1/max(eigen(S1)$values);S2=jpen(S1,2*c/p,gstr/(gstr-1)*(c/p-1/gstr))} else {S2=S}; if(length(lam)==0) {lam=2*gam/p};
E=(solve(S2)+gam*diag(dim(S2)[1]))/(1+gam);   E1=sign(E)*pmax(abs(E)-lam/(2*(1+gam)),0);diag(E1)=diag(E);  return(E1); }
# CV based tunning parameter selection
jpen.tune=function (Ytr, gama, lambda = NULL) 
  {
    p = dim(Ytr)[2]
    if (length(lambda) == 0) {
      Str1 = var(Ytr)
      diag(Str1) = 0
      c = max(abs(Str1));
      L = length(gama);
      l2 = length(lamvec(c, gama[1], p))
      errl = matrix(0, ncol = 3, nrow = L * l2)
    } else {
      L = length(gama)
      l2 = length(lambda)
      errl = matrix(0, ncol = 3, nrow = L * l2)
    }
    K = 5
    n = dim(Ytr)[1]
    p = dim(Ytr)[2]
    cv = f.K.fold(n, 5)
    for (l in 1:L) {
      st = (l - 1) * l2 + 1
      en = st + l2 - 1;errl[st:en,2]=gama[l];
      if (length(lambda) == 0) {
        errl[st:en, 1] = c(lamvec(c, gama[l], p))
      } else {
        errl[st:en, 1] = c(lambda)
      }
    }
    for (r in 1:(L * l2)) {
      for (k in 1:K) {
        n1 = floor(n * (1 - 1/log(n)))
        id = sample(1:n, n1, replace = FALSE)
        ytr = Ytr[id, ]
        yts = Ytr[-id, ]
        Str = var(ytr)
        Sts = var(yts)
        lam = errl[r, 1]
        gam = errl[r, 2]
        Est1 = jpen(Str, gam, lam)
        Est1[abs(Est1) < 1e-04] = 0
        if (min(eigen(Est1)$value) < 1e-04) {
          errl[r, 3] = 10^10
        }
        else {
          errl[r, 3] = errl[r, 3] + norm(Est1-Sts,type="F");
        }
      }
    }
    k = which.min(errl[, 3])
    opt = errl[k, -3]
    return(c(opt))
  }
  
# tuning parameter selection for inverse covariance matrix
jpen.inv.tune=function(Ytr,gama,lambda=NULL) { p=ncol(Ytr);if(length(lambda)==0)     {Str1=var(Ytr);diag(Str1)=0;c=max(abs(Str1[row(Str1)!=col(Str1)]));L=length(gama);        l2=length(lamvec(c,gama[1],p)); errl=matrix(0,ncol=3,nrow=L*l2);} else {L=length(gama);l2=length(lambda);errl=matrix(0,ncol=3,nrow=L*l2);}
K=5;n=dim(Ytr)[1];p=dim(Ytr)[2];cv=f.K.fold(n,5); for(l in 1:L) {st=(l-1)*l2+1;en=st+l2-1;errl[st:en,2]=gama[l];if(length(lambda)==0) {errl[st:en,1]=c(lamvec(c,gama[l],p))} else {errl[st:en,1]=c(lambda)}};   for(r in 1:(L*l2)) {for(k in 1:K)  
{ id=cv[[k]]$train; ytr=Ytr[id,];yts=Ytr[cv[[k]]$test,];Str=var(ytr);Sts=var(yts);lam=errl[r,1];gam=errl[r,2];  Est1=jpen.inv(Str,lam,gam);Est1[abs(Est1)<.0001]=0;                                                                                           
if(min(eigen(Est1)$value)<.0001) { errl[r,3]=10^10} else {errl[r,3]=errl[r,3]+log(-sum(dmvnorm(yts, mean = rep(0, dim(Est1)[1]), sigma = solve(Est1),log = TRUE)))}}};k=which.min(errl[,3]);opt=errl[k,-3];   return(c(opt));}; 
  
