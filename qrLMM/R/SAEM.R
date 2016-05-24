##########################################################################################
#QUANTILE REGRESSION FOR LINEAR MIXED MODEL
##########################################################################################

QSAEM_COM_7 =  function(y,x,z,nj,p,precision=0.0001,MaxIter=300,M=20,pc=0.5,beta=beta,sigmae=sigmae,D=D)
{
  start.time <- Sys.time()
  
  n = length(nj)
  N = sum(nj)
  d = dim(x)[2]
  q = dim(z)[2]
  
  delta1 = 0.001
  delta2 = precision
  
  #assymetry parameters
  vp = (1-2*p)/(p*(1-p))
  tp = sqrt(2/(p*(1-p)))
  
  MDel   = MElim(q)
  ndiag  = (q*(1+q)/2)
  npar   = d+1+ndiag
  
  critval   = 1
  critval2  = 1
  count     = 0
  teta      = c(beta,sigmae,D[upper.tri(D, diag = T)]) 
  tetam     = matrix(data=NA,nrow=npar,ncol=MaxIter)
  
  EPV       = matrix(0,nrow = npar,ncol = MaxIter) 
  
  if(pc==1){
    seqq=rep(1,pc*MaxIter)
  }else{
    seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
    seqq = c(rep(1,MaxIter-length(seqq)),seqq)
  }
  
  
  SAEM_bb   = array(data=0,dim=c(MaxIter+1,n,q,q))
  SAEM_bi   = array(data=0,dim=c(MaxIter+1,n,q))
  
  SAEM_VGi  = array(data=0,dim=c(MaxIter+1,n,npar))
  
  SAEM_ui   = vector("list", n)
  SAEM_Dui  = vector("list", n)
  SAEM_bbZD = vector("list", n)
  SAEM_DZb  = vector("list", n)
  
  for(j in 1:n)
  {
    SAEM_bb[count+1,j,,] = diag(q)
    
    SAEM_ui[[j]]   = array(data=0,dim=c(MaxIter+1,nj[j]))
    SAEM_Dui[[j]]  = array(data=0,dim=c(MaxIter+1,nj[j],nj[j]))
    SAEM_bbZD[[j]] = array(data=0,dim=c(MaxIter+1,q,nj[j]))
    SAEM_DZb[[j]]  = array(data=0,dim=c(MaxIter+1,nj[j]))
  }
  
  pb = tkProgressBar(title = "QRLMM via SAEM", min = 0,max = MaxIter, width = 300)
  setTkProgressBar(pb, 0, label=paste("Iter ",0,"/",MaxIter,"     -     ",0,"% done",sep = ""))
  
  while(critval < 3 && critval2 < 3)
  {
    count  = count + 1
    
    sumb1  = matrix(data=0,nrow=d,ncol=1)
    sumb2  = matrix(data=0,nrow=d,ncol=d)
    sumD   = matrix(data=0,nrow=q,ncol=q)
    sumsig = 0
    IE     = 0
    
    for (j in 1:n)
    { 
      y1=y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
      x1=matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=d)
      z1=matrix(z[(sum(nj[1:j-1])+1):(sum(nj[1:j])),],ncol=q)
      v1s = matrix(data=1,nrow=nj[j],ncol=1)
      
      ##########################################################################
      #PASSO E
      ##########################################################################
      
      ui       = matrix(0,nrow = nj[j],ncol = 1)
      sum_ui   = matrix(0,nrow = nj[j],ncol = 1)
      Dui      = matrix(0,nrow = nj[j],ncol = nj[j])
      sum_Dui  = matrix(0,nrow = nj[j],ncol = nj[j])
      
      sum_bb   = matrix(0,nrow = q,ncol = q)
      sum_bbZD = matrix(0,nrow = q,ncol = nj[j])
      sum_DZb  = matrix(0,nrow = nj[j],ncol = 1)
      
      VGi   = matrix(data = 0,nrow = npar,ncol = M)
      
      bmetro = matrix(data = MHbi2(j=j,M=M,y1,x1,z1,bi=as.matrix(SAEM_bi[count,j,]),bibi=as.matrix(SAEM_bb[count,j,,]),d=d,q=q,p=p,nj=nj,beta=beta,sigmae=sigmae,D=D),nrow = q,ncol = M)
      
      for(l in 1:M)
      { 
        for(k in 1:nj[j])
        {
          chi = (as.numeric(y1[k]-x1[k,]%*%beta-z1[k,]%*%bmetro[,l])^2)/(sigmae*tp^2)
          psi = (tp^2)/(4*sigmae)
          
          ui[k]    = Egig(lambda = 0.5,chi = chi,psi = psi,func = "x")
          Dui[k,k] = Egig(lambda = 0.5,chi = chi,psi = psi,func = "1/x")
        }
        sum_ui   = sum_ui   + ui
        sum_Dui  = sum_Dui + Dui
        sum_bbZD = sum_bbZD + bmetro[,l]%*%t(bmetro[,l])%*%t(z1)%*%Dui
        sum_bb   = sum_bb + bmetro[,l]%*%t(bmetro[,l])
        sum_DZb  = sum_DZb + Dui%*%z1%*%bmetro[,l]
        
        t_G_b       = t(x1)%*%(Dui%*%(y1-x1%*%beta) - Dui%*%z1%*%bmetro[,l] - vp*v1s)
        t_G_sig     = -(3/2)*(nj[j])*(1/sigmae) + (1/(2*(tp^2)*(sigmae^2)))*(t(y1-x1%*%beta-z1%*%bmetro[,l])%*%Dui%*%(y1-x1%*%beta-z1%*%bmetro[,l]) - 2*vp*t(y1-x1%*%beta-z1%*%bmetro[,l])%*%v1s + ((tp^4)/4)*t(ui)%*%v1s)
        t_G_D       = bmetro[,l]%*%t(bmetro[,l])
        
        GG1         = (1/(sigmae*tp^2))*t_G_b
        GG2         = t_G_sig
        GG3         = (1/2)*MElim(q)%*%(kronecker(X = solve(D),Y = solve(D)))%*%as.vector(t_G_D-D)
        VGi[,l]     = rbind(GG1,GG2,GG3)  
      }
      
      E_ui   = sum_ui/M
      E_Dui = sum_Dui/M  
      E_bbZD = sum_bbZD/M
      E_bb   = sum_bb/M
      E_bi   = apply(bmetro,1,mean)
      E_DZb  = sum_DZb/M
      E_VGi  = apply(VGi,1,mean)
      
      SAEM_ui[[j]][count+1,]   = SAEM_ui[[j]][count,] + seqq[count]*(E_ui - SAEM_ui[[j]][count,])
      SAEM_Dui[[j]][count+1,,] = SAEM_Dui[[j]][count,,] + seqq[count]*(E_Dui - SAEM_Dui[[j]][count,,])
      SAEM_bbZD[[j]][count+1,,] = SAEM_bbZD[[j]][count,,] + seqq[count]*(E_bbZD - SAEM_bbZD[[j]][count,,])    
      SAEM_bb[count+1,j,,]  = SAEM_bb[count,j,,] + seqq[count]*(E_bb - SAEM_bb[count,j,,])    
      SAEM_bi[count+1,j,]  = SAEM_bi[count,j,] + seqq[count]*(E_bi - SAEM_bi[count,j,])    
      SAEM_DZb[[j]][count+1,] = SAEM_DZb[[j]][count,] + seqq[count]*(E_DZb - SAEM_DZb[[j]][count,])
      
      SAEM_VGi[count+1,j,] = SAEM_VGi[count,j,] + seqq[count]*(E_VGi - SAEM_VGi[count,j,]) 
      
      ##########################################################################
      #PASSO M
      ##########################################################################
      
      #PASSO M betas
      sumb1 = sumb1 + (t(x1)%*%(SAEM_Dui[[j]][count+1,,]%*%y1 - SAEM_DZb[[j]][count+1,] - vp*v1s))
      sumb2 = sumb2 + (t(x1)%*%SAEM_Dui[[j]][count+1,,]%*%x1)
      
      #PASSO M sigmae
      
      tsig = t(y1-x1%*%beta)%*%SAEM_Dui[[j]][count+1,,]%*%(y1-x1%*%beta) - 
        2*t(y1-x1%*%beta)%*%SAEM_DZb[[j]][count+1,] + 
        tr(z1%*%SAEM_bbZD[[j]][count+1,,])- 
        2*vp*t(y1-x1%*%beta)%*%v1s + 
        2*vp*t(z1%*%SAEM_bi[count+1,j,])%*%v1s + 
        ((tp^4)/4)*t(SAEM_ui[[j]][count+1,])%*%v1s
      
      sumsig  = sumsig + as.numeric(tsig)
      
      #PASSO M matriz D 
      
      sumD = sumD + SAEM_bb[count+1,j,,]
      
      #SUM do prod vector gradiente
      
      IE = IE + SAEM_VGi[count+1,j,]%*%t(SAEM_VGi[count+1,j,])
    }
    
    beta     = solve(sumb2)%*%sumb1 
    D        = sumD/n
    sigmae   = sumsig/(3*N*tp^2)
    EP     = sqrt(diag(solve(IE)))
    EPV[,count] = EP
    
    param     = teta
    teta      = c(beta,sigmae,D[upper.tri(D, diag = T)])
    criterio  = abs(teta-param)/(abs(param)+delta1)
    criterio2 = abs(teta-param)/(EP+0.0001)
    if(max(criterio) < delta2){critval=critval+1}else{critval=0}
    if(max(criterio2) < 0.0002){critval2=critval2+1}else{critval2=0}
    
    #PRUEBAS
    #############################################################################
    tetam[,count] = teta
    setTkProgressBar(pb, count, label=paste("Iter ",count,"/",MaxIter,"     -     ",round(count/MaxIter*100,0),"% done",sep = ""))
    
    if  (count == MaxIter){critval=10}
  }
  
  loglik = logveroIS(beta,sigmae,D,y,x,z,nj,bi=SAEM_bi[count+1,,],bibi=SAEM_bb[count+1,,,],MIS=500,n=n,d=d,q=q,p=p)
  AIC    = -2*loglik +2*npar
  BIC    = -2*loglik +log(N)*npar
  HQ     = -2*loglik +2*log(log(N))*npar
  table  = data.frame(beta,EP[1:d],beta-(1.96*EP[1:d]),beta+(1.96*EP[1:d]),beta/EP[1:d],2*pnorm(abs(beta/EP[1:d]),lower.tail = F))
  rownames(table) = paste("beta",1:d)
  colnames(table) = c("Estimate","Std. Error","Inf CI95%","Sup CI95%","z value","Pr(>|z|)")
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  res     = list(iter = count,criterio = max(criterio,criterio2),beta = beta,sigmae= sigmae,D = D,EP=EP,table = table,loglik=loglik,AIC=AIC,BIC=BIC,HQ=HQ,time = time.taken)
  conv    = list(teta = tetam[,1:count],EPV = EPV[,1:count])
  obj.out = list(conv=conv,res = res)
  
  if  (count == MaxIter)
  {
    setTkProgressBar(pb, MaxIter, label=paste("MaxIter reached ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    Sys.sleep(1)
    close(pb)
  }
  else
  {
    setTkProgressBar(pb, MaxIter, label=paste("Convergence at Iter ",count,"/",MaxIter,"    -    100 % done",sep = ""))
    Sys.sleep(1)
    close(pb)
  }
  
  class(obj.out)  =  "QSAEM_COM_7"
  return(obj.out)
}