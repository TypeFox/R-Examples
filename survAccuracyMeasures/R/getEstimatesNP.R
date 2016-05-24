

getEstimatesNP <- function(data, 
                         cutpoint,  
                         measures,
                         predict.time,
                         CalVar = TRUE, subcohort=FALSE)
{
  
  
  
  

  N = nrow(data); ## cohort size
  data0 <- data #unsorted data
  data = data [order(data$yi),] ## sorted by marker
  
  ck = data$yi;     
  wgtk = data$wi;  
  xk = data$xi; 
  sk = data$si; 
  psk=data$psi; 
  
  c0 <- cutpoint; 
  t0 <- predict.time
  

  nc = length(ck); # sampled size 
  ind.ck = (1:nc)[order(ck)]
  CWk = WGT.FUN(data[,c(1,2)],data0, t0=t0)  # use full cohort xi di to calculate censoring weights
  
  St0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk >= t0))/sum(CWk*wgtk)
  Ft0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk <  t0))/sum(CWk*wgtk)
  Fck = sum.I(ck,">=",ck,wgtk*CWk)/sum(CWk*wgtk)
  St0 = max(St0.Fck)            ## St0     = P(T> t0)
  Ft0 = 1-St0                   ## Ft0     = P(T<=t0)
  FPR.ck= (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)
  TPR.ck= (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
  NPV.ck= St0.Fck/Fck           ## P(T> t0|Y<=ck)
  PPV.ck= (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)
  AUC = sum(TPR.ck*(FPR.ck-c(FPR.ck[-1],0)))
  
  ## acc.ck: accuracy estimates at c
  nm.acc = c("FPR","TPR","NPV","PPV"); 
  acc.ck = data.frame("cutoff"=ck,"FPR"=FPR.ck,"TPR"=TPR.ck,"NPV"=NPV.ck, "PPV"=PPV.ck)    
  

  if(!is.null(c0)){
    tmpind.c0 = sum.I(c0, ">=", ck); 
    acc.c0 = acc.ck[tmpind.c0,]; 
    F.c0 = Fck[tmpind.c0]
  }else{ acc.c0 = acc.ck }
  est= data.frame("AUC" = AUC, acc.c0[-1])

  est  = est[,measures]
  #est= list("AUC" = AUC, "ACC.u0"=acc.uk[tmpind.u0,-c(1,ind0)],"ACC.c0"=acc.c0, "ACC.all" = acc.ck) ##output all cutoff
  if(!CalVar){
   list("estimates" = est)
  }else{
    
    #this portion of the code is not functional
    
    ###### Variance calculation below ########
   # Phi=Phi.C.new.FUN(xk=data$xi,dk=data$di, Ti=data0$xi, Di=data0$di, t0 = t0)
    Phi = NA #this portion of the code is not functional
    ## doing u0 and c0 together
    c.u0 = NULL; acc.u0.temp=NULL; F.c0.b =F.c0
    
    CC = c(c.u0,c0); #nu0 = length(u0); 
    acc.c0.b = rbind(acc.u0.temp,acc.c0); 
    
    U.ACC.c0.tmp = as.list(1:4); 
    U.ACC.c0 =Wexp.c0= as.list(1:4); 
    names(U.ACC.c0.tmp) = names(U.ACC.c0)=nm.acc; n.acc.c=length(U.ACC.c0)
    
    
    I.ck.c0 = 1*(ck>=VTM(CC,nc)); 
    U.ACC.c0.tmp$FPR =  (xk >  t0)*(I.ck.c0-  VTM(acc.c0.b$FPR,nc))/St0      ## exp for FPRhat(c)-FPR(c)
    U.ACC.c0.tmp$TPR =  (xk <= t0)*(I.ck.c0-  VTM(acc.c0.b$TPR,nc))/(1-St0)  ## exp for TPRhat(c)-TPR(c)
    U.ACC.c0.tmp$NPV =  (1-I.ck.c0)*(1*(xk> t0)-VTM(acc.c0.b$NPV,nc))/VTM(F.c0.b,nc)
    U.ACC.c0.tmp$PPV =    I.ck.c0*(1*(xk<=t0)-VTM(acc.c0.b$PPV,nc))/(1-VTM(F.c0.b,nc))
    
    U.AUC = (xk<=t0)/(1-St0)*(1-FPR.ck-AUC)+(xk>t0)/St0*(TPR.ck-AUC)
    
    Wexp.np.AUC = CWk*U.AUC+Phi%*%(wgtk*CWk*U.AUC)/sum(wgtk)
    se.auc = sqrt(Est.Var.CCH.trueweights(N,Wexp.np.AUC,data,data$si, subcohort=FALSE))
    se.u0 = NULL

    se.c0 = NULL
    if (!is.null(c0)) {
      npc = length(U.ACC.c0)	
      for(kk in 1:npc){ 
        
          U.ACC.c0[[kk]] = U.ACC.c0.tmp[[kk]]
        Wexp.c0[[kk]] = CWk*U.ACC.c0[[kk]]+Phi%*%(wgtk*CWk*U.ACC.c0[[kk]])/sum(wgtk)
        se.c0 = c(se.c0,sqrt(Est.Var.CCH.trueweights(N,data.frame(Wexp.c0[[kk]]),data,data$si, subcohort=subcohort)))
      }
      se.c0 = data.frame(matrix(se.c0,nrow=length(c0)))
      names(se.c0) = nm.acc    
    }
   
    se <- data.frame("AUC" = se.auc, se.c0)
    se <- se[,measures] 
    list("estimates" = est,"se" =se) 
  }	
  
  
}