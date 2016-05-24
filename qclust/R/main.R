###########################################################
dmm_fn <- function(x,means,vars,pis)
{
  n=dim(x)[1]
  m=length(pis)
  p=dim(means)[2]
  r=matrix(NA,n,m)
  for (i in 1:m)
  {
    if (p!=1)
    {
      r[,i] = pis[i]*mvtnorm::dmvnorm(x, means[i,], vars[,,i])
    }
    else
    {
      r[,i] = pis[i]*dnorm(x, means[i,], sqrt(vars[,,i]))
    }
  }
  return(r)
}
###########################################################
dmm_noweights_fn <- function(x,means,vars,pis)
{
  n=dim(x)[1]
  m=length(pis)
  p=dim(means)[2]
  r=matrix(NA,n,m)
  for (i in 1:m)
  {
    if (p!=1)
    {
      r[,i] = mvtnorm::dmvnorm(x, means[i,], vars[,,i])
    }
    else
    {
      r[,i] = dnorm(x, means[i,], sqrt(vars[,,i]))
    }
  }
  return(r)
}
###########################################################
Lq_fn <- function(x,q)
{
  if (q!=1)
  {
    return( (x^(1-q)-1)/(1-q) )
  }
  else
  {
    return( log(x))
  }
}
###########################################################
Qclust <- function(d, K = NULL, modelNames = NULL,q)
{
  
  mfit=mclust::Mclust(d, K, modelNames)
  
  if (q==1)
  {
    return(mfit)    
  } else {
    if (dim(d)[2]>1)
    {
      mu_initial=as.matrix(t(mfit$parameters$mean))
      var_initial=as.array(mfit$parameters$variance$sigma)
      pi_initial=as.matrix(mfit$parameters$pro,length(mfit$parameters$pro),1)
    } else {
      mu_initial=matrix(mfit$parameters$mean,mfit$G,1)
      var_initial=array(mfit$parameters$variance$sigma,dim=c(1,1,mfit$G))
      pi_initial=as.matrix(mfit$parameters$pro,length(mfit$parameters$pro),1)  
    }
    qfit=qclust_w_initialvalues(d,mu_initial,var_initial,pi_initial,q,tol=0.00001)
    if (dim(d)[2]>1)
    {
      mfit$parameters$mean=t(qfit$means)
      mfit$parameters$variance$sigma=qfit$vars
      mfit$parameters$pro=as.vector(qfit$pis)
    } else {
      mfit$parameters$mean=t(qfit$means)
      mfit$parameters$variance$sigmasq=as.vector(qfit$vars)
      mfit$parameters$pro=as.vector(qfit$pis)
    }
    mfit$LqLikelihood=qfit$lq_likelihood_ls[length(qfit$lq_likelihood_ls)]
    mfit$LqLikelihood_ls=qfit$lq_likelihood_ls
    mfit$z=qfit$z
    mfit$qBIC=-2*mfit$LqLikelihood+log(mfit$n)*mfit$df
    mfit$q=q
    mfit$classification=qfit$class
    return(mfit)
  }
}
###########################################################
qclust_w_initialvalues <- function(d,means_init,vars_init,pis_init,q,tol=0.00001)
{
  # d is n by p matrix
  # means_init is m by p matrix
  # vars_init is p by p by m array
  # pis_init is m by 1 matrix
  
  n=dim(d)[1]
  p=dim(d)[2]
  m=length(pis_init)
  means=means_init
  vars=vars_init
  pis=pis_init
  diff=100
  it_max=500
  it=1
  means_ls=array(rep(means_init,it_max+1),dim=c(dim(means_init),it_max))
  vars_ls=array(rep(vars_init,it_max+1),dim=c(dim(vars_init),it_max))
  pis_ls=array(rep(pis_init,it_max+1),dim=c(length(pis_init),it_max))
  means_ls[,,1]=means
  vars_ls[,,,1]=vars
  pis_ls[,1]=pis
  lq_likelihood_old=-99999999999
  lq_likelihood_new=99999999999
  lq_likelihood_ls=rep(NA, it_max)
  while (it<=it_max && diff > tol)
  {
    comp_pdfs=dmm_fn(d,means,vars,pis)
    overall_pdf=apply(comp_pdfs,1,sum)
    
    # E step
    z=comp_pdfs/(overall_pdf%*%matrix(1,1,m))
    z[is.nan(z)]=1/m
    z=z^q
    
    # M step
    w_M=z*(comp_pdfs^(1-q))
    w_M=w_M/matrix(1,n,1)%*%apply(w_M,2,sum)
    for (i in 1:m)
    {
      means[i,]=apply(d*(w_M[,i]%*%matrix(1,1,p)),2,sum)
      d_centered=d-matrix(1,n,1)%*%means[i,]
      d_centered_weighted=(sqrt(w_M[,i])%*%matrix(1,1,p))*d_centered
      vars[,,i]=t(as.matrix(d_centered_weighted))%*%as.matrix(d_centered_weighted)  
      #if (!is.non.singular.matrix(vars[,,i]))
      #{
      #    vars[1,2,i]=0
      #    vars[2,1,i]=0
      #} 
    }
    pis=apply(z*(dmm_noweights_fn(d,means,vars,pis)^(1-q)),2,sum)^(1/q)
    pis=as.matrix(pis/sum(pis))
    comp_pdfs=dmm_fn(d,means,vars,pis)
    overall_pdf=apply(comp_pdfs,1,sum)
    lq_likelihood_new=sum(Lq_fn(overall_pdf,q))
    
    diff=lq_likelihood_new-lq_likelihood_old
    lq_likelihood_old=lq_likelihood_new
    means_ls[,,it]=means
    vars_ls[,,,it]=vars
    pis_ls[,it]=pis
    lq_likelihood_ls[it]=lq_likelihood_new
    it=it+1
    # cat("it =",it," ")
  }
  it=it-1
  cat("EMLq ends at it =",it,"\n")
  means_ls=means_ls[,,1:it]
  vars_ls=vars_ls[,,,1:it]
  pis_ls=pis_ls[,1:it]
  lq_likelihood_ls=lq_likelihood_ls[1:it]
  class=rep(0,n)
  z=matrix(0,n,m)
  for (i in 1:n)
  {
    class[i]=which(comp_pdfs[i,]==max(comp_pdfs[i,]))[1]
    z[i,]=comp_pdfs[i,]/sum(comp_pdfs[i,])
  }
  output=list(means=means,
              vars=vars,
              pis=pis,
              lq_likelihood_ls=lq_likelihood_ls,
              class=class,
              z=z,
              par_ls=list(means_ls=means_ls,vars_ls=vars_ls,pis_ls=pis_ls)
  )
  return(output)
}
###########################################################