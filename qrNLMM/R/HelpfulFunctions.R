##########################################################################################
#PAUSE BETWEEN PLOTS
##########################################################################################

readkey <- function()
{
  cat ("Press [enter] to continue")
  line <- readline()
}

####################################################################
####################################################################

minbeta = function(y,x,fixed,covar,p,q,nlmodel)
{
  dif = y - nlmodel(x,fixed,rep(0,q),covar)
  return((sum(p*dif[dif>0]) - sum((1-p)*dif[dif<0])))
}

#######################################################################################
#######################################################################################

densbiv=function(j,bi,beta,sigmae,D,p,y1,x1,cov1,q=q,nj=nj,nlmodel=nlmodel)
{
  bi = matrix(bi,nrow = q,ncol = 1)
  dprod = 1
  for(k in 1:nj[j])
  {
    dprod = dprod * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,bi,cov1[k,]),sigma=sigmae,p=p))
  }
  dens=as.numeric(dprod*dmvnorm(x=as.numeric(bi),mean=matrix(rep(0,q),q,1),sigma=D))
  return(dens)
}

####################################################################
####################################################################

MHbi2 = function(j,M,y1,x1,cov1,bi,bibi,d,q,p,nj,beta=beta,sigmae=sigmae,D=D,nlmodel=nlmodel)
{  
  E_bi = bi
  V_bi = bibi - E_bi%*%t(E_bi)
  GEN = matrix(NA,nrow=q,ncol=(M+1))
  count = 1
  GEN[,1]=rmvnorm(n = 1,mean = E_bi,sigma=V_bi)
  
  while(count <= M)
  {
    cand = rmvnorm(n=1,mean=as.vector(GEN[,count]),sigma=V_bi)
    
    c1 = densbiv(j=j,bi=cand,beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=GEN[,count],mean=as.vector(cand),sigma=V_bi)
    c2 = densbiv(j=j,bi=GEN[,count],beta=beta,sigmae=sigmae,D=D,p=p,y1,x1,cov1,q,nj,nlmodel=nlmodel)*dmvnorm(x=as.vector(cand),mean=as.vector(GEN[,count]),sigma=V_bi)
    alfa = c1/c2
    
    if(is.nan(alfa)+0==1) {alfa=0.0001}
    if (runif(1) < min(alfa, 1))
    {
      count = count + 1
      GEN[,count] = cand
    }
  }
  return(GEN[,2:(M+1)])
}

####################################################################
####################################################################

MElim = function(q)
{
  Dtest = matrix(data = 0,ncol = q,nrow = q)
  Dtest[upper.tri(Dtest, diag = T)]=1
  vDtest = as.vector(Dtest)
  Elim = matrix(data = 0,nrow = (q*(1+q)/2),ncol = q^2)
  count=1
  for(i in 1:q^2)
  {
    if(vDtest[i]==1)
    {
      Elim[count,i]=1
      count = count +1
    }
  }
  return(Elim)
}

####################################################################
####################################################################

####################################################################
####################################################################

logveroIS = function(beta,sigmae,D,y,x,covar,nj,bi,bibi,MIS,n,d,q,p,n.covar,nlmodel=nlmodel)
{
  logvero = 0
  bi = matrix(data = bi,nrow = n,ncol = q)
  
  for(j in 1:n)
  {
    y1  = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]#APPROVED
    x1  = matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=1)#APPROVED
    if(n.covar>0){cov1 = matrix(covar[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=n.covar)}
    
    if(q==1){bibij=bibi[j]}else{bibij=bibi[j,,]}
    
    E_bi = bi[j,]
    V_bi = bibij - E_bi%*%t(E_bi)
    Bgen = rmvnorm(n = MIS,mean = E_bi,sigma=matrix(V_bi,q,q))
    sum = 0
    
    for(l in 1:MIS)
    {
      multden = 1
      for(k in 1:nj[j])
      {
        multden = multden * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,Bgen[l,],cov1[k,]),sigma=sigmae,p=p))
      }
      multt = multden*(dmvnorm(x = Bgen[l,],mean=rep(0,q),sigma=D)/dmvnorm(x = Bgen[l,],mean = E_bi,sigma=V_bi))
      sum = sum + multt
    }
    prom    = sum/MIS
    logvero = logvero + log(prom)
  }
  return(logvero)
}

####################################################################
####################################################################

logveroMC = function(beta,sigmae,D,y,x,covar,nj,MC,n,d,q,p,n.covar,nlmodel=nlmodel)
{
  logvero = 0
  
  for(j in 1:n)
  {
    y1  = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]#APPROVED
    x1  = matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j]))],ncol=1)#APPROVED
    if(n.covar>0){cov1 = matrix(covar[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=n.covar)}
    
    Bgen = rmvnorm(n = MC,mean=rep(0,q),sigma=matrix(D,q,q))
    sum = 0
    
    for(l in 1:MC)
    {
      multden = 1
      for(k in 1:nj[j])
      {
        multden = multden * (dALD(y=y1[k],mu=nlmodel(x1[k],beta,Bgen[l,],cov1[k,]),sigma=sigmae,p=p))
      }
      sum = sum + multden
    }
    prom    = sum/MC
    logvero = logvero + log(prom)
  }
  return(logvero)
}

###############################################################################
###############################################################################

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

###############################################################################
###############################################################################

validate_str = function(prueba)
{
  prueba = gsub("\n","   ",prueba)
  prueba = gsub("[[:punct:]]","   ",gsub("[|[[:digit:]|]]","    ",prueba))
  prueba = gsub("\\d","   ",prueba)
  prueba = paste("  ",prueba,"  ",sep = "")
  functions = c(" fixed "," random "," covar "," exp "," x "," log "," sinh "," cosh "," sqrt "," pnorm "," dnorm "," asin "," acos "," atan "," lgamma "," digamma "," trigamma "," psigamma "," sin "," cos "," tan "," gamma ")
  for(i in 1:length(functions))
  {
    prueba = gsub(functions[i],"  ",prueba)
    
  }
  prueba = gsub("[[:blank:]]"," ",prueba)
  prueba = gsub("^ *|(?<= ) | *$","", prueba, perl=T)
  return(prueba)
}

###############################################################################
###############################################################################


group.plot = function(x,y,groups,...){
  if(length(y) != length(groups)) stop("groups does not match with  the provided data. (length(y) != length(groups))")
  if(length(x) != length(groups)) stop("groups does not match with  the provided data. (length(x) != length(groups))")
  if(length(y) != length(x)) stop("the provided data does not match. (length(y) != length(x))")
  nj = c(as.data.frame(table(as.factor(c(groups))))[,2])
  
  Ymatrix = Xmatrix = matrix(data = NA,nrow = length(nj),ncol = max(nj))
  
  for(j in 1:length(nj))
  {
    Xmatrix[j,(1:nj[j])] = x[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
    Ymatrix[j,(1:nj[j])] = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
  }
  matplot(x = t(Xmatrix),y = t(Ymatrix),...)
}

#group.plot(x = Time,y = weight,groups = Plot)


group.lines = function(x,y,groups,...){
  if(length(y) != length(groups)) stop("groups does not match with  the provided data. (length(y) != length(groups))")
  if(length(x) != length(groups)) stop("groups does not match with  the provided data. (length(x) != length(groups))")
  if(length(y) != length(x)) stop("the provided data does not match. (length(y) != length(x))")
  nj = c(as.data.frame(table(as.factor(c(groups))))[,2])
  Ymatrix = Xmatrix = matrix(data = NA,nrow = length(nj),ncol = max(nj))
  
  for(j in 1:length(nj))
  {
    Xmatrix[j,(1:nj[j])] = x[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
    Ymatrix[j,(1:nj[j])] = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
  }
  matlines(x = t(Xmatrix),y = t(Ymatrix),...)
}

#group.lines(x = Time,y = weight,groups = Plot)


group.points = function(x,y,groups,...){
  if(length(y) != length(groups)) stop("groups does not match with  the provided data. (length(y) != length(groups))")
  if(length(x) != length(groups)) stop("groups does not match with  the provided data. (length(x) != length(groups))")
  if(length(y) != length(x)) stop("the provided data does not match. (length(y) != length(x))")
  nj = c(as.data.frame(table(as.factor(c(groups))))[,2])
  Ymatrix = Xmatrix = matrix(data = NA,nrow = length(nj),ncol = max(nj))
  
  for(j in 1:length(nj))
  {
    Xmatrix[j,(1:nj[j])] = x[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
    Ymatrix[j,(1:nj[j])] = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
  }
  matpoints(x = t(Xmatrix),y = t(Ymatrix),...)
}

#group.points(x = Time,y = weight,groups = Plot)