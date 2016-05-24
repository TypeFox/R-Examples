QSAEM_NL = function(y,x,nj,initial,exprNL,covar=NA,p=0.5,precision = 0.0001,M=20,pc=0.25,MaxIter=500,beta,sigmae,D,nlmodel,d,q)
{
  start.time <- Sys.time()
  
  n = length(nj)
  N = sum(nj)
  
  #NEW CODE
  if(all(is.na(covar)==TRUE)){n.covar = 0}else{n.covar = dim(covar)[2]}
  
  delta1 = 0.001
  delta2 = precision
  
  #assimetry
  vp = (1-2*p)/(p*(1-p))
  tp = sqrt(2/(p*(1-p)))
  
  MDel   = MElim(q)
  ndiag  = (q*(1+q)/2)
  npar   = d+1+ndiag
  
  ###################################
  
  exprFX = gsub("|]","",gsub("[|[]","",as.character(exprNL)))
  paste1 = paste("\"fixed",1:d,"\"",sep = "",collapse = ",")
  paste2 = paste("fixed",1:d,sep = "",collapse = ",")
  pasteA = paste("covar",1:n.covar,sep = "",collapse = ",")
  paste3 = ifelse(n.covar==0,"",pasteA)
  paste4 = ifelse(n.covar==0,paste("x",paste2,sep =","),paste("x",paste2,paste3,sep =","))
  paste5 = paste("random",1:q,sep = "",collapse = ",")
  paste6 = paste("deriv( ~ ",exprFX,",c(",paste1,"),function(",paste(paste4,",",paste5,sep = ""),"){})",sep = "")
  dd     = eval(parse(text = paste6))
  
  critval  = 1
  count    = 0
  teta     = c(beta,sigmae,D[upper.tri(D, diag = T)]) 
  tetam    = matrix(data=NA,nrow=npar,ncol=MaxIter)
  
  EPV      = matrix(0,nrow = npar,ncol = MaxIter)
  
  if(pc==1)
  {
    seqq=rep(1,pc*MaxIter)
  }else
    
  {
    seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
    seqq = c(rep(1,MaxIter-length(seqq)),seqq)
  }
  
  SAEM_bb    = array(data=0,dim=c(MaxIter+1,n,q,q))
  SAEM_bb2   = array(data=0,dim=c(MaxIter+1,n,q,q))
  SAEM_bi    = array(data=0,dim=c(MaxIter+1,n,q))
  SAEM_VGi   = array(data=0,dim=c(MaxIter+1,n,npar))
  SAEM_sig   = array(data=0,dim=c(MaxIter+1,n))
  SAEM_Hbeta = array(data=0,dim=c(MaxIter+1,n,d,d))
  
  for(j in 1:n){SAEM_bb[count+1,j,,] = diag(q)}
  
  pb = tkProgressBar(title = "QRNLMM via SAEM", min = 0,max = MaxIter, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",MaxIter,"     -     ",0,"% done",sep = ""))
  
  while(critval < 3)
  {
    count  = count + 1
    
    sumb1  = matrix(data=0,nrow=d,ncol=1)
    sumb2  = matrix(data=0,nrow=d,ncol=d)
    sumD   = matrix(data=0,nrow=q,ncol=q)
    sumsig = 0
    IE     = 0
    
    x = as.matrix(x)
    
    for (j in 1:n)
    { 
      y1  = y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
      x1  = matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=1)
      if(n.covar>0){cov1 = matrix(covar[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=n.covar)}
      v1s = matrix(data=1,nrow=nj[j],ncol=1)
      
      ##########################################################################
      #PASSO E
      ##########################################################################
      
      ui     = matrix(0,nrow = nj[j],ncol = 1)
      Dui    = matrix(0,nrow = nj[j],ncol = nj[j])
      VGi    = matrix(data = 0,nrow = npar,ncol = M)
      Msig   = array(NA,M)
      Mbb    = array(NA,dim = c(q,q,M))
      Mbb2   = array(NA,dim = c(q,q,M))
      Hbeta  = array(NA,dim = c(d,d,M))
      
      bmetro = matrix(MHbi2(j=j,M=M,y1,x1,cov1,bi=as.matrix(SAEM_bi[count,j,]),bibi=as.matrix(SAEM_bb[count,j,,]),d=d,q=q,p=p,nj=nj,beta=beta,sigmae=sigmae,D=D,nlmodel=nlmodel),q,M)
            
      paste6 = paste("beta[",1:d,"]",sep = "",collapse = ",")
      paste7 = ifelse(n.covar==0,"",paste(",",paste("cov1[,",1:n.covar,"]",sep = "",collapse = ","),sep=""))
      paste8 = paste("attr(dd(x1,",paste6,paste7,",",paste(rep(0,q),sep="",collapse = ","),"),\"gradient\")",sep = "")
      grad   = eval(parse(text = paste8))
      
      for(l in 1:M)
      { 
        for(k in 1:nj[j])
        {
          chi = ((y1[k] - nlmodel(x1[k],beta,bmetro[,l],cov1[k,]))^2)/(sigmae*tp^2)
          if(chi==0){chi=10^-100}
          psi = (tp^2)/(4*sigmae)
          
          ui[k]    = Egig(lambda = 0.5,chi = chi,psi = psi,func = "x")
          Dui[k,k] = Egig(lambda = 0.5,chi = chi,psi = psi,func = "1/x")
        }
        
        Msig[l]    = t(y1)%*%Dui%*%y1 - 2*vp*t(y1)%*%v1s + ((tp^4)/4)*t(ui)%*%v1s -
          2*t(y1)%*%Dui%*%nlmodel(x1,beta,bmetro[,l],cov1) +
          2*vp*t(v1s)%*%nlmodel(x1,beta,bmetro[,l],cov1) +
          t(nlmodel(x1,beta,bmetro[,l],cov1))%*%Dui%*%nlmodel(x1,beta,bmetro[,l],cov1)
        
        Mbb[,,l]   = bmetro[,l]%*%t(bmetro[,l])
        
        GG1        = -(1/(sigmae*tp^2))*t(grad)%*%(2*vp*v1s + 2*Dui%*%nlmodel(x1,beta,bmetro[,l],cov1) -2*Dui%*%y1)
        GG2        = -(3/2)*(nj[j])*(1/sigmae) + (1/(2*(tp^2)*(sigmae^2)))*Msig[l]
        GG3        = (1/2)*MElim(q)%*%(kronecker(X = solve(D),Y = solve(D)))%*%as.vector(Mbb[,,l]-D)
        VGi[,l]    = rbind(GG1,GG2,GG3)
        
        Hbeta[,,l]  = -(1/(sigmae*tp^2))*t(grad)%*%Dui%*%grad
      }
      
      E_bi    = apply(bmetro,1,mean)
      E_bb    = apply(Mbb,c(1,2),mean)
      E_sig   = mean(Msig)
      E_VGi   = apply(VGi,1,mean)
      E_Hbeta = apply(Hbeta,c(1,2),mean)
      
      SAEM_bb[count+1,j,,]  = SAEM_bb[count,j,,] + seqq[count]*(E_bb - SAEM_bb[count,j,,])
      SAEM_sig[count+1,j]  = SAEM_sig[count,j] + seqq[count]*(E_sig - SAEM_sig[count,j])    
      SAEM_bi[count+1,j,]  = SAEM_bi[count,j,] + seqq[count]*(E_bi - SAEM_bi[count,j,])
      SAEM_VGi[count+1,j,] = SAEM_VGi[count,j,] + seqq[count]*(E_VGi - SAEM_VGi[count,j,]) 
      SAEM_Hbeta[count+1,j,,] = SAEM_Hbeta[count,j,,] + seqq[count]*(E_Hbeta - SAEM_Hbeta[count,j,,]) 
      
      ##########################################################################
      #PASSO M
      ##########################################################################
      
      #PASSO M sigmae
      sumsig  = sumsig + SAEM_sig[count+1,j]
      
      #PASSO M matriz D 
      sumD = sumD + SAEM_bb[count+1,j,,]
      
      #PASSO M betas
      sumb1 = sumb1 + SAEM_VGi[count+1,j,1:d]#Gradiente
      sumb2 = sumb2 + SAEM_Hbeta[count+1,j,,]#Hessiana
      
      #SUM do prod vector gradiente
      IE = IE + SAEM_VGi[count+1,j,]%*%t(SAEM_VGi[count+1,j,])
    }
    
    ssumb2 = tryCatch(solve(-sumb2), error=function(e) stop("Error trying to compute an inverse of a matrix involved in beta. Frequently, for Nonlinear models this happens when the initial values are not properly."))
    
    beta     =  beta + ssumb2%*%sumb1
    D        = sumD/n
    sigmae   = sumsig/(3*N*tp^2)
    EP     = sqrt(diag(tryCatch(solve(IE), error=function(e) diag(npar))))
    EPV[,count] = EP
    
    param    = teta
    teta     = c(beta,sigmae,D[upper.tri(D, diag = T)])
    criterio = abs(teta-param)/(abs(param)+delta1)
    if(max(criterio) < delta2){critval=critval+1}else{critval=0}
    #PRUEBAS
    #############################################################################
    tetam[,count] = teta
    setTkProgressBar(pb, count, label=paste("Iter ",count,"/",MaxIter,"     -     ",round(count/MaxIter*100,0),"% done",sep = ""))
    
    if  (count == MaxIter){critval=10}
  }
  
  #loglik = logveroMC(beta,sigmae,D,y,x,covar,nj,MC=1000,n=n,d=d,q=q,p=p,nlmodel=nlmodel)
  loglik = logveroIS(beta,sigmae,D,y,x,covar,nj,bi=SAEM_bi[count+1,,],bibi=SAEM_bb[count+1,,,],MIS=500,n=n,d=d,q=q,p=p,n.covar,nlmodel=nlmodel)
  
  AIC    = -2*loglik +2*npar
  BIC    = -2*loglik +log(N)*npar
  HQ     = -2*loglik +2*log(log(N))*npar
  table  = data.frame(beta,EP[1:d],beta/EP[1:d],2*pnorm(abs(beta/EP[1:d]),lower.tail = F))
  rownames(table) = paste("beta",1:d)
  colnames(table) = c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  res     = list(quantile = p,nlmodel = nlmodel,iter = count,criterio = max(criterio),beta = beta,sigmae= sigmae,D = D,EP=EP,table = table,loglik=loglik,AIC=AIC,BIC=BIC,HQ=HQ,time = time.taken)
  conv    = list(teta = tetam[,1:count],EPV = EPV[,1:count])
  obj.out = list(conv=conv,res = res)
  
  if  (count == MaxIter)
  {
    setTkProgressBar(pb, MaxIter, label=paste("MaxIter reached ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    #Sys.sleep(2)
    close(pb)
  }
  else
  {
    setTkProgressBar(pb, MaxIter, label=paste("Convergence at Iter ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    #Sys.sleep(2)
    close(pb)
  }
  
  class(obj.out)  =  "QSAEM_NL"
  return(obj.out)
}

####################################################################
####################################################################