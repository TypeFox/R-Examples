# Paulo Cortez @ 2013

isbest=function(Cur,Best,metric)
{if(length(Cur)==0 || is.na(Cur)) return (FALSE) 
 else{ 
      if(is.list(metric) && !is.null(metric$metric) ) metric=metric$metric
      if(is.function(metric)) return (Cur<Best)
      else { return (switch(metric,
                CRAMERV=,F1=,MCC=,ALIFTATPERC=,NALIFT=,ALIFT=,ACC=,ACCLASS=,KAPPA=,COR=,PRECISION=,TPR=,KENDALL=,SPEARMAN=,TOLERANCE=,
                TNR=,R2=,R22=,NAREC=,NAUC=,TPRATFPR=,AUCCLASS=,AUC=Cur>Best,
                Cur<Best))
           }
     } 
}
worst=function(metric)
{ if(is.list(metric) && !is.null(metric$metric) ) metric=metric$metric
  if(is.function(metric)) return (Inf) # assumption that metric should follow (i.e. "<" is better)
  else{
  return (switch(metric, 
                CRAMERV=,F1=,PRECISION=,TPR=,TNR=,ALIFTATPERC=,NALIFT=,ALIFT=,KAPPA=,ACC=,ACCLASS=,NAUC=,AUCCLASS=,AUC=,TPRATFPR=,TOLERANCE=,NAREC=-0.1, # [0-1 or 100]
                MCC=,COR=,KENDALL=,SPEARMAN=-1.1, # [-1,1]
                q2=1.1,
                R2=-0.1,
                R22=-Inf,
                CE=,BER=100.1, # [0,100], %
                Inf)) # other regression metrics
      }
}

# x - vector of predictions, y - vector of desired values  
# metric - metric or vector of metrics
metrics=function(y,x=NULL,D=0.5,TC=-1,val=NULL,aggregate="no")
{ warning("Deprecated function, please use instead: mmetric(y,x,metric=\"ALL\",D,TC,val,aggregate)")
  return(mmetric(y,x,metric="ALL",D,TC,val,aggregate))
}

is.mmetric=function(metric)
{ M=c("MAEO","MSEO","KENDALL",
      "ACC","CE","BER","KAPPA","CRAMERV","ACCLASS","TPR","TNR","PRECISION","F1","MCC",
      "ACC","CE","BER","KAPPA","CRAMERV","ACCLASS","TPR","TNR","PRECISION","F1","MCC","BRIER","BRIERCLASS","AUC","AUCCLASS","NAUC","TPRATFPR","ALIFT","NALIFT","ALIFTATPERC",
      "SAE","MAE","MdAE","GMAE","MaxAE","RAE","SSE","MSE","MdSE","RMSE","GMSE","HRMSE","RSE","RRSE","ME","COR","q2","R2","Q2","NAREC","TOLERANCE","MAPE","MdAPE","RMSPE","RMdSPE","SMAPE","SMdAPE","SMinkowski3","MMinkowski3","MdMinkowski3"
      ) 
  return(sum(M==metric)>0)
}

mmetric=function(y,x=NULL,metric,D=0.5,TC=-1,val=NULL,aggregate="no")
{ 
 if(is.list(y) && is.null(x)) # mining object
     { # special metrics: sum all?
       if(sum(metric=="CONF")>0) CONF=TRUE else CONF=FALSE
       if(sum(metric=="ROC")>0) ROC=TRUE else ROC=FALSE
       if(sum(metric=="LIFT")>0) LIFT=TRUE else LIFT=FALSE
       if(sum(metric=="REC")>0) REC=TRUE else REC=FALSE
       R=y$runs;res=NULL
       if(is.function(metric)) LM=1 else LM=length(metric)
       if(aggregate=="no") 
       {
        if(CONF||ROC||LIFT||REC){res=vector("list",length=R); for(i in 1:R) res[[i]]=mmetric(y$test[[i]],y$pred[[i]],metric=metric,D=D,TC=TC,val=val)}
        #else if(LM==1) 
        #    { res=vector(length=R); for(i in 1:R) res[i]=mmetric(y$test[[i]],y$pred[[i]],metric=metric,D=D,TC=TC,val=val) }
        else{ 
              res=matrix(nrow=R,ncol=LM)
              for(i in 1:R) { if(i==1){ aux=mmetric(y$test[[i]],y$pred[[i]],metric=metric,D=D,TC=TC,val=val);
                                        res=matrix(nrow=R,ncol=length(aux))
                                        res[i,]=aux
                                      }
                              else res[i,]=mmetric(y$test[[i]],y$pred[[i]],metric=metric,D=D,TC=TC,val=val) 
                            }
              if(length(names(aux))==0) names(aux)=metric[1] # new line
              res=data.frame(res)
              if(LM==1 && metric!="ALL" && !is.factor(y$test[[1]][1])) {C=length(aux);naux=vector(length=C);for(i in 1:C) naux[i]=paste(metric,i,sep="");names(res)=naux;}
              else names(res)=names(aux)
            }
       }
       else if(aggregate=="sum"||aggregate=="mean") 
       {
           # does not work for roc, lift and rec! 
           for(i in 1:R) 
           {
              aux=mmetric(y$test[[i]],y$pred[[i]],metric=metric,D=D,TC=TC,val=val)
              if(LM==1 && CONF)
               {
                if(i==1)
                {
                  res=aux$conf
                }
                else{ res=res+aux$conf }
               }
              else{ if(i==1) res=aux else res=res+aux }
           }
           if(aggregate=="mean") res=res/R
       }
       return(res)
     }
else if(is.function(metric)) return(metric(y,x)) # resmode=1
else if(NCOL(y)>1) # ranking 
  {
    res=NULL;nres=NULL;
    if(length(metric)==1 && metric=="ALL") metric=c("KENDALL","SPEARMAN")
    if(sum(metric=="KENDALL")>0) KENDALL=TRUE else KENDALL=FALSE # -1 to 1
    if(sum(metric=="SPEARMAN")>0) SPEARMAN=TRUE else SPEARMAN=FALSE
    LM=length(metric)
    if(LM==1) rsingle=TRUE else rsingle=FALSE

    C=NCOL(y);Total=NROW(y) 
    if(KENDALL) kendall=0
    if(SPEARMAN) spearman=0
    for(k in 1:Total)
       {
         if(KENDALL) kendall=kendall+cor(y[k,],x[k,],method="kendall")
         if(SPEARMAN) spearman=spearman+cor(y[k,],x[k,],method="spearman")
       }
    if(KENDALL){kendall=kendall/Total;if(rsingle) return(kendall) else{res=c(res,kendall);nres=c(nres,"KENDALL")}} else kendall=NULL
    if(SPEARMAN){spearman=spearman/Total;if(rsingle) return(spearman) else{res=c(res,spearman);nres=c(nres,"SPEARMAN")}} else spearman=NULL

    if(!is.null(res)) {names(res)=nres; 
                       # sort res:
                       I=NULL # for regression, this works perfectly?
                       for(i in 1:LM)
                        {
                          ii=which(nres==metric[i])[1]
                          if(!is.na(ii)) I=c(I,ii) 
                        }
                       res=res[I]
                      }
   return(res)
  }
else if(is.factor(y)) # classification
  {
    res=NULL;nres=NULL;
    if(length(metric)==1 && metric=="ALL")
         { 
           if(is.ordered(y)) metric=c("MAEO","MSEO","KENDALL")
           else if(is.factor(x)) metric=c("ACC","CE","BER","KAPPA","CRAMERV","ACCLASS","TPR","TNR","PRECISION","F1","MCC") 
           else metric=c("ACC","CE","BER","KAPPA","CRAMERV","ACCLASS","TPR","TNR","PRECISION","F1","MCC","BRIER","BRIERCLASS","AUC","AUCCLASS","NAUC","TPRATFPR","ALIFT","NALIFT","ALIFTATPERC")
         }
    LM=length(metric)
    if(sum(metric=="CONF")>0) CONF=TRUE else CONF=FALSE
    if(sum(metric=="ROC")>0) ROC=TRUE else ROC=FALSE
    if(sum(metric=="LIFT")>0) LIFT=TRUE else LIFT=FALSE
    if(sum(metric=="ACC")>0) ACC=TRUE else ACC=FALSE
    if(sum(metric=="CE"|metric=="MER")>0) CE=TRUE else CE=FALSE
    if(sum(metric=="MAEO")>0) MAEO=TRUE else MAEO=FALSE
    if(sum(metric=="MSEO")>0) MSEO=TRUE else MSEO=FALSE
    if(sum(metric=="KENDALL")>0) KENDALL=TRUE else KENDALL=FALSE
    if(sum(metric=="BER")>0) BER=TRUE else BER=FALSE 
    if(sum(metric=="KAPPA")>0) KAPPA=TRUE else KAPPA=FALSE
    if(sum(metric=="CRAMERV")>0) CRAMERV=TRUE else CRAMERV=FALSE 

    if(sum(metric=="ACCLASS")>0) ACCLASS=TRUE else ACCLASS=FALSE
    if(sum(metric=="TPR")>0) TPR=TRUE else TPR=FALSE 
    if(sum(metric=="TNR")>0) TNR=TRUE else TNR=FALSE 
    if(sum(metric=="PRECISION")>0) PRECISION=TRUE else PRECISION=FALSE 
    if(sum(metric=="F1")>0) F1=TRUE else F1=FALSE 
    if(sum(metric=="MCC")>0) MCC=TRUE else MCC=FALSE

    if(!is.factor(x)) # prob
        {
         if(sum(metric=="BRIER")>0) BRIER=TRUE else BRIER=FALSE 
         if(sum(metric=="BRIERCLASS")>0) BRIERCLASS=TRUE else BRIERCLASS=FALSE 
         if(sum(metric=="AUC")>0) AUC=TRUE else AUC=FALSE 
         if(sum(metric=="AUCCLASS")>0) AUCCLASS=TRUE else AUCCLASS=FALSE 
         if(sum(metric=="NAUC")>0) NAUC=TRUE else NAUC=FALSE 
         if(sum(metric=="TPRATFPR")>0) TPRATFPR=TRUE else TPRATFPR=FALSE 
         if(sum(metric=="ALIFT")>0) ALIFT=TRUE else ALIFT=FALSE 
         if(sum(metric=="NALIFT")>0) NALIFT=TRUE else NALIFT=FALSE 
         if(sum(metric=="ALIFTATPERC")>0) ALIFTATPERC=TRUE else ALIFTATPERC=FALSE 
        } 

    if(LM==1) rsingle=TRUE else rsingle=FALSE
    if(CONF||ROC||LIFT) reslist=TRUE # list
    else reslist=FALSE# named vector

    if(TC>0) C=2 else C=length(levels(y[1]))
    Total=length(y)
    if(CONF||ACC||CE||MAEO||MSEO||BER||KAPPA||CRAMERV||KENDALL||ACCLASS||TPR||TNR||PRECISION||F1||MCC) # conf
    { conf=Conf(y,x,D=D,TC=TC,predreturn=TRUE)
      if(CRAMERV||KENDALL) pred=conf$pred
      conf=conf$conf
    } else conf=NULL

    if(ACC||CE||KAPPA){diag=0;diagr=rep(0,C);}
    if(MAEO) maeo=0 else maeo=NULL
    if(MSEO) mseo=0 else mseo=NULL
    if(BER) ber=vector(length=C)
    if(ACCLASS) acclass=rep(0,C) else acclass=NULL
    if(TPR||F1) tpr=rep(0,C) else tpr=NULL
    if(TNR) tnr=rep(0,C) else tnr=NULL
    if(PRECISION||F1) precision=rep(0,C) else precision=NULL
    if(MCC) mcc=rep(0,C) else mcc=NULL
    if(F1) f1=rep(0,C) else f1=NULL
    if(ACC||CE||MAEO||MSEO||KAPPA||BER||ACCLASS||TPR||TNR||PRECISION||F1||MCC)
    {
         for(k in 1:C) 
         { 
           if(ACC||CE||KAPPA) diag=diag+conf[k,k]
           if(KAPPA||BER) sum_conf_k=sum(conf[k,])
           if(KAPPA) diagr[k]=sum_conf_k*(sum(conf[,k])/Total)
           if(MAEO||MSEO) { for(i in 1:C) 
                            {
                             err=conf[k,i]*(i-k)
                             if(MAEO) maeo=maeo+abs(err)
                             if(MSEO) mseo=mseo+(err)^2
                            }
                          }
           if(BER) ber[k]=sum(conf[k,-k])/sum_conf_k
           if(ACCLASS||TPR||TNR||PRECISION||F1||MCC)
           {TP=conf[k,k]
	    FN=0
 	    for(i in 1:C) # iterator?
		if(i!=k) FN=FN+conf[k,i]
            FP=0
 	    for(i in 1:C) # iterator?
		if(i!=k) FP=FP+conf[i,k]
            TN=Total-TP-FN-FP
            if(ACCLASS) acclass[k]=100*(TP+TN)/Total 
            if((TPR||F1) && TP!=0) tpr[k]=100*TP/(FN+TP)
            if(TNR && TN!=0) tnr[k]=100*TN/(TN+FP) 
            if((PRECISION||F1) && TP!=0) precision[k]=100*TP/(TP+FP)
            if(F1 && precision[k]!=0 && tpr[k]!=0) f1[k]=2*((precision[k]*tpr[k])/(precision[k]+tpr[k]))
            if(MCC) { mcc[k]=TP*TN+FP*FN; if(mcc[k]!=0) mcc[k]=mcc[k]/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) else mcc[k]=0 }
           }
         }
    }
    if(ACC||CE) {acc=c(diag/Total)*100; 
                 if(rsingle && ACC) return(acc) else if(ACC){res=c(res,acc);nres=c(nres,"ACC")}
               } else acc=NULL # total accuracy, in percentage
    if(CE) {ce=100-acc; if(rsingle) return(ce) else{res=c(res,ce);nres=c(nres,"CE")}} else ce=NULL
    if(KAPPA) {kappa=100*(diag-sum(diagr))/(Total-sum(diagr)); if(rsingle) return(kappa) else{res=c(res,kappa);nres=c(nres,"KAPPA")}} else kappa=NULL # G
    if(BER) {ber=100*mean(ber);if(rsingle) return(ber) else{res=c(res,ber);nres=c(nres,"BER")}}else ber=NULL # G
    if(MAEO){maeo=maeo/Total;if(rsingle) return(maeo) else{res=c(res,maeo);nres=c(nres,"MAEO")}} 
    if(MSEO){mseo=mseo/Total;if(rsingle) return(mseo) else{res=c(res,mseo);nres=c(nres,"MSEO")}} 
    if(CRAMERV) # single
     {
      #T=try(chisq.test(y,pred)$statistic[[1]],silent=TRUE)
      T=suppressWarnings(try(chisq.test(y,pred)$statistic[[1]],silent=TRUE))
      if(!is.numeric(T)) T=0
      cramerv=sqrt(T/(Total*(C-1)))
      if(rsingle) return(cramerv) else{res=c(res,cramerv);nres=c(nres,"CRAMERV")}
    } else cramerv=NULL
    if(KENDALL && is.ordered(y)) 
     { 
       if(!is.ordered(pred)) pred=ordered(pred,levels=levels(y[1]))
       c=0;d=0;et=0;ep=0
       for(k in 1:(Total-1))
         for(i in (k+1):Total)
           {
             if( (y[k]<y[i] && pred[k]<pred[i]) || (y[k]>y[i] && pred[k]>pred[i]) ) c=c+1
             else if( (y[k]<y[i] && pred[k]>pred[i]) || (y[k]>y[i] && pred[k]<pred[i]) ) d=d+1
             else if( y[k]==y[i] && (pred[k]>pred[i]||pred[k]<pred[i]) ) et=et+1
             else if( pred[k]==pred[i] && (y[k]>y[i]||y[k]<y[i]) ) ep=ep+1
             # else # ignore
           }
       kendall=(c-d)/(sqrt(c+d+et)*sqrt(c+d+ep)) 
       if(rsingle) return(kendall) else{res=c(res,kendall);nres=c(nres,"KENDALL")}
     } else kendall=NULL

    if( (ACCLASS||TPR||TNR||PRECISION||F1||MCC)||(!is.factor(x)&& (BRIERCLASS||AUCCLASS)) ) naux=vector(length=C) else NAUX=FALSE

    if(ACCLASS){if(rsingle) return(acclass) else {res=c(res,acclass); for(i in 1:C) naux[i]=paste("ACCLASS",i,sep=""); nres=c(nres,naux)}}
    if(TPR){if(rsingle) return(tpr) else {res=c(res,tpr); for(i in 1:C) naux[i]=paste("TPR",i,sep=""); nres=c(nres,naux)}}
    if(TNR) {if(rsingle) return(tnr) else {res=c(res,tnr); for(i in 1:C) naux[i]=paste("TNR",i,sep=""); nres=c(nres,naux)}}
    if(PRECISION) {if(rsingle) return(precision) else {res=c(res,precision); for(i in 1:C) naux[i]=paste("PRECISION",i,sep=""); nres=c(nres,naux)}}
    if(F1) {if(rsingle) return(f1) else {res=c(res,f1); for(i in 1:C) naux[i]=paste("F1",i,sep=""); nres=c(nres,naux)}}
    if(MCC) {if(rsingle) return(mcc) else {res=c(res,mcc); for(i in 1:C) naux[i]=paste("MCC",i,sep=""); nres=c(nres,naux)}}
 
    if(!is.factor(x)) # prob
        {
         if(sum(metric=="BRIER")>0) BRIER=TRUE else BRIER=FALSE 
         if(sum(metric=="BRIERCLASS")>0) BRIERCLASS=TRUE else BRIERCLASS=FALSE 
         if(sum(metric=="AUC")>0) AUC=TRUE else AUC=FALSE 
         if(sum(metric=="AUCCLASS")>0) AUCCLASS=TRUE else AUCCLASS=FALSE 
         if(sum(metric=="NAUC")>0) NAUC=TRUE else NAUC=FALSE 
         if(sum(metric=="TPRATFPR")>0) TPRATFPR=TRUE else TPRATFPR=FALSE 
         if(sum(metric=="ALIFT")>0) ALIFT=TRUE else ALIFT=FALSE 
         if(sum(metric=="NALIFT")>0) NALIFT=TRUE else NALIFT=FALSE 
         if(sum(metric=="ALIFTATPERC")>0) ALIFTATPERC=TRUE else ALIFTATPERC=FALSE 

         if(BRIER||BRIERCLASS)
                  { 
                    L=levels(y[1])
                    if(TC>0)
                    {
                     T=as.numeric(y==L[TC])
                     tbrier=sum((T-x[,TC])^2)/Total # Brier=MSE
                     if(BRIERCLASS) brier=c(tbrier,tbrier) 
                    }
                    else # TC== -1, all!
                    {
                     prop=table(y)[]/Total; 
                     tbrier=0;brier=vector(length=C) 
                     for(i in 1:C) 
                     {
                      T=as.numeric(y==L[i])
                      brier[i]=sum((T-x[,i])^2)/Total # Brier=MSE
                      if(prop[i]>0) tbrier=tbrier+prop[i]*brier[i]
                     }
                    }
                    if(BRIER) { if(rsingle) return(tbrier) else {res=c(res,tbrier);nres=c(nres,"BRIER")}}
                    if(BRIERCLASS){if(rsingle) return(brier) else{res=c(res,brier); for(i in 1:C) naux[i]=paste("BRIERCLASS",i,sep=""); nres=c(nres,naux)}}
                  } else {brier=NULL;tbrier=NULL}
         if(ROC||AUC||AUCCLASS||NAUC||TPRATFPR)
                { if(TC<1) TC2=C else TC2=TC
                  if(C>2) {roc=ROCcurve(y,x);tauc=roc$auc;auc=vector(length=C); for(i in 1:C) auc[i]=roc$roc[[i]]$auc;}
                     else {roc=twoclassROC(y,x[,TC2],Positive=levels(y[1])[TC2]);tauc=roc$auc;auc=c(tauc,tauc) }
                  if(AUC){if(rsingle) return(tauc) else{res=c(res,tauc);nres=c(nres,"AUC")}}
                  if(AUCCLASS){if(rsingle) return(auc) else{res=c(res,auc);for(i in 1:C) naux[i]=paste("AUCCLASS",i,sep="");nres=c(nres,naux)}}
                  if(NAUC){ if(is.null(val)) val2=1 
                            else if(is.list(val)) val2=val[[which(metric=="NAUC")[1]]]
                            else if(length(val)>1) val2=val[which(metric=="NAUC")[1]]
                            else val2=val
                            if(val2>1) val2=1
                            if(C>2) roc2=partialcurve(roc$roc[[TC2]]$roc,val2) else roc2=partialcurve(roc$roc,val2)
                            nauc=curvearea(roc2,val2)
                            if(rsingle) return(nauc) else{res=c(res,nauc);nres=c(nres,"NAUC")}
                          } else nauc=NULL
                  if(TPRATFPR) {  if(is.null(val)) val2=1 
                                  else if(is.list(val)) val2=val[[which(metric=="TPRATFPR")[1]]]
                                  else if(length(val)>1) val2=val[which(metric=="TPRATFPR")[1]]
                                  else val2=val
                                  if(val2>1) val2=1
                                  if(C>2) roc2=partialcurve(roc$roc[[TC2]]$roc,val2)
                                  else roc2=partialcurve(roc$roc,val2)
                                  if(is.vector(roc2)) tpratfpr=roc2[2] else tpratfpr=roc2[nrow(roc2),2]
                                  if(rsingle) return(tpratfpr) else{res=c(res,tpratfpr);nres=c(nres,"TPRATFPR")}
                          } else tpratfpr=NULL
                } else {roc=NULL;auc=NULL;tauc=NULL;nauc=NULL;tpratfpr=NULL}
         if(LIFT||ALIFTATPERC||ALIFT||NALIFT)
                {
                  if(TC<1) TC2=C else TC2=TC
                  lift=LIFTcurve(y,x,TC=TC2) # does not work for more than 2 classes
                  if(ALIFT) {alift=lift$area;if(rsingle) return(alift) else{res=c(res,alift);nres=c(nres,"ALIFT")}}
                  else alift=NULL
                  if(NALIFT) { 
                               if(is.null(val)) val2=1 
                               else if(is.list(val)) val2=val[[which(metric=="NALIFT")[1]]]
                               else if(length(val)>1) val2=val[which(metric=="NALIFT")[1]]
                               else val2=val
                               if(val2>1) val2=1
                               lift2=partialcurve(lift$alift,val2); nalift=curvearea(lift2,val2)
                               if(rsingle) return(nalift) else{res=c(res,nalift);nres=c(nres,"NALIFT")}
                               } else nalift=NULL
                  if(ALIFTATPERC)
                             {
                               if(is.null(val)) val2=1 
                               else if(is.list(val)) val2=val[[which(metric=="ALIFTATPERC")[1]]]
                               else if(length(val)>1) val2=val[which(metric=="ALIFTATPERC")[1]]
                               else val2=val
                               if(val2>1) val2=1
                               lift2=partialcurve(lift$alift,val2)
                               if(is.vector(lift2)) aliftatperc=lift2[2] else aliftatperc=lift2[nrow(lift2),2]
                               if(rsingle) return(aliftatperc) else{res=c(res,aliftatperc);nres=c(nres,"ALIFTATPERC")}
                             } else aliftatperc=NULL 
                } else {lift=NULL;alift=NULL;nalift=NULL;aliftatperc=NULL}
        }
    else {roc=NULL;tauc=NULL;auc=NULL;brier=NULL;tbrier=NULL;nauc=NULL;tpratfpr=NULL;lift=NULL;alift=NULL;nalift=NULL;aliftatperc=NULL}

    if(!is.null(res)) { names(res)=nres;
                        I=NULL # for classification: problem with multi-class!!! XXX
                        Lnres=length(nres)
                        nmetric=vector(length=Lnres) 
                        i=1;k=1;stop=FALSE
                        while(!stop)
                        {
                         m=metric[i] 
                         if(m!="ACCLASS"&&m!="TPR"&&m!="TNR"&&m!="PRECISION"&&m!="F1"&&m!="MCC"&&m!="AUCCLASS"&&m!="BRIERCLASS"){nmetric[k]=m;k=k+1;}
                         else {
                               for(j in 1:C) {nmetric[k]=paste(metric[i],j,sep="");k=k+1;}
                              }
                         i=i+1;
                         if(i>LM) stop=TRUE
                        }
                        I=NULL
                        for(i in 1:Lnres)
                        {
                          ii=which(nres==nmetric[i])[1]
                          if(!is.na(ii)) I=c(I,ii) 
                        }
                        res=res[I]
                      }

    if(reslist) {res=list(res=res,conf=conf,roc=roc,lift=lift)}
    return(res)
 } # y factor
 else # regression
 {
    # absolute measures:
    res=NULL;nres=NULL;LM=length(metric)
    if(length(metric)==1 && metric=="ALL") metric=c("SAE","MAE","MdAE","GMAE","MaxAE","RAE","SSE","MSE","MdSE","RMSE","GMSE","HRMSE","RSE","RRSE","ME","COR","q2","R2","Q2","NAREC","TOLERANCE","MAPE","MdAPE","RMSPE","RMdSPE","SMAPE","SMdAPE","SMinkowski3","MMinkowski3","MdMinkowski3")

    LM=length(metric)
    if(sum(metric=="SAE")>0) SAE=TRUE else SAE=FALSE
    if(sum(metric=="MAE")>0) MAE=TRUE else MAE=FALSE
    if(sum(metric=="MdAE")>0) MdAE=TRUE else MdAE=FALSE
    if(sum(metric=="GMAE")>0) GMAE=TRUE else GMAE=FALSE
    if(sum(metric=="MaxAE")>0) MaxAE=TRUE else MaxAE=FALSE
    if(sum(metric=="RAE")>0) RAE=TRUE else RAE=FALSE
    # Square measures:
    if(sum(metric=="SSE")>0) SSE=TRUE else SSE=FALSE
    if(sum(metric=="MSE")>0) MSE=TRUE else MSE=FALSE
    if(sum(metric=="MdSE")>0) MdSE=TRUE else MdSE=FALSE
    if(sum(metric=="RMSE")>0) RMSE=TRUE else RMSE=FALSE
    if(sum(metric=="GMSE")>0) GMSE=TRUE else GMSE=FALSE
    if(sum(metric=="HRMSE")>0) HRMSE=TRUE else HRMSE=FALSE
    if(sum(metric=="RSE")>0) RSE=TRUE else RSE=FALSE
    if(sum(metric=="RRSE")>0) RRSE=TRUE else RRSE=FALSE
    if(sum(metric=="ME")>0) ME=TRUE else ME=FALSE
    # Cor measures:
    if(sum(metric=="COR")>0) COR=TRUE else COR=FALSE
    if(sum(metric=="q2")>0) q2=TRUE else q2=FALSE

    if(sum(metric=="R2")>0) R2=TRUE else R2=FALSE
    if(sum(metric=="R22")>0) R22=TRUE else R22=FALSE
    if(sum(metric=="Q2")>0) LQ2=TRUE else LQ2=FALSE
    # REC measures: 
    if(sum(metric=="REC")>0) REC=TRUE else REC=FALSE
    if(sum(metric=="NAREC")>0) NAREC=TRUE else NAREC=FALSE
    if(sum(metric=="TOLERANCE")>0) TOLERANCE=TRUE else TOLERANCE=FALSE
    # forecasting measures:
    if(sum(metric=="MdAPE")>0) MdAPE=TRUE else MdAPE=FALSE
    if(sum(metric=="RMSPE")>0) RMSPE=TRUE else RMSPE=FALSE
    if(sum(metric=="RMdSPE")>0) RMdSPE=TRUE else RMdSPE=FALSE
    if(sum(metric=="MAPE")>0) MAPE=TRUE else MAPE=FALSE
    if(sum(metric=="SMAPE")>0) SMAPE=TRUE else SMAPE=FALSE
    if(sum(metric=="SMdAPE")>0) SMdAPE=TRUE else SMdAPE=FALSE
    if(sum(metric=="MRAE")>0) MRAE=TRUE else MRAE=FALSE
    if(sum(metric=="MdRAE")>0) MdRAE=TRUE else MdRAE=FALSE
    if(sum(metric=="GMRAE")>0) GMRAE=TRUE else GMRAE=FALSE
    if(sum(metric=="THEILSU2")>0) THEILSU2=TRUE else THEILSU2=FALSE
    if(sum(metric=="MASE")>0) MASE=TRUE else MASE=FALSE
    # Minkowski errors
    if(sum(metric=="SMinkowski3")>0) SMINKOWSKI3=TRUE else SMINKOWSKI3=FALSE
    if(sum(metric=="MMinkowski3")>0) MMINKOWSKI3=TRUE else MMINKOWSKI3=FALSE
    if(sum(metric=="MdMinkowski3")>0) MdMINKOWSKI3=TRUE else MdMINKOWSKI3=FALSE
    if(LM==1) rsingle=TRUE else rsingle=FALSE
    if(REC) reslist=TRUE # list
    else reslist=FALSE# named vector

    if(SAE||MAE||MdAE||GMAE||MaxAE||RAE||MAPE||SMAPE||SMdAPE||MdAPE||ME||SSE||MSE||MdSE||RMSE||RSE||RRSE||GMSE||R22||LQ2||MRAE||MdRAE||GMRAE||RMSPE||RMdSPE||THEILSU2||MASE||SMINKOWSKI3||MMINKOWSKI3||MdMINKOWSKI3) err=y-x
    if(SAE||MAE||MdAE||GMAE||MaxAE||RAE||SMAPE||SMdAPE||SMINKOWSKI3||MMINKOWSKI3||MdMINKOWSKI3) eabs=abs(err)
    if(SAE||RAE) {sae=sum(eabs); if(rsingle && SAE) return(sae) else if(SAE){res=c(res,sae);nres=c(nres,"SAE")}} else sae=NULL
    if(MAE) {mae=mean(eabs);     if(rsingle) return(mae) else{res=c(res,mae);nres=c(nres,"MAE")}} else mae=NULL
    if(MdAE) {mdae=median(eabs); if(rsingle) return(mdae) else{res=c(res,mdae);nres=c(nres,"MdAE")}} else mdae=NULL
    if(GMAE) {gmae=prod(eabs)^(1/(length(eabs)));  if(rsingle) return(gmae) else{res=c(res,gmae);nres=c(nres,"GMAE")}} else gmae=NULL
    if(MaxAE) {maxae=max(eabs);  if(rsingle) return(maxae) else{res=c(res,maxae);nres=c(nres,"MaxAE")}} else maxae=NULL
 
    if(SSE||MSE||MdSE||RMSE||GMSE||RSE||RRSE||R22||LQ2||THEILSU2) esqr=(err)^2
    if(SSE||RSE||RRSE||R22||LQ2) {sse=sum(esqr);if(rsingle && SSE) return(sse) else if(SSE){res=c(res,sse);nres=c(nres,"SSE")}} else sse=NULL
    if(MSE||RMSE||THEILSU2) {mse=mean(esqr);if(rsingle && MSE) return(mse) else if(MSE){res=c(res,mse);nres=c(nres,"MSE")}} else mse=NULL
    if(RMSE||THEILSU2) {rmse=sqrt(mse);if(rsingle && RMSE) return(rmse) else if(RMSE){res=c(res,rmse);nres=c(nres,"RMSE")}} else rmse=NULL 
    if(MdSE) {mdse=median(esqr);if(rsingle) return(mdse) else{res=c(res,mdse);nres=c(nres,"MdSE")}} else mdse=NULL
    if(GMSE) {gmse=prod(esqr)^(1/(length(esqr)));if(rsingle) return(gmse) else{res=c(res,gmse);nres=c(nres,"GMSE")}} else gmse=NULL
    if(HRMSE) {hrmse=sqrt( mean((1-(x/y))^2) ) ;if(rsingle) return(hrmse) else{res=c(res,hrmse);nres=c(nres,"HRMSE")}} else hrmse=NULL

    if(ME) { me=mean(err);if(rsingle && ME) return(me) else{res=c(res,me);nres=c(nres,"ME")}} else me=NULL

    if(RAE||RSE||RRSE||R22||LQ2) {ymean=mean(y)}

    if(RAE) {rae=100*sae/sum(abs(y-ymean));if(rsingle) return(rae) else{res=c(res,rae);nres=c(nres,"RAE")}} else rae=NULL
    if(RSE||RRSE||R22||LQ2) {sum_ym_esqr=sum((y-ymean)^2)}

    if(RSE) {rse=100*sse/sum_ym_esqr;if(rsingle) return(rse) else{res=c(res,rse);nres=c(nres,"RSE")}} else rse=NULL
    if(RRSE){rrse=100*sqrt(sse/sum_ym_esqr);if(rsingle) return(rrse) else {res=c(res,rrse);nres=c(nres,"RRSE")}} else rrse=NULL
    if(R22) {r22=1-sse/sum_ym_esqr;if(rsingle) return(r22) else {res=c(res,r22);nres=c(nres,"R22")}} else r22=NULL
# problem with this formulation, check better:
    if(LQ2) {Q2=sse/sum_ym_esqr;if(rsingle) return(Q2) else {res=c(res,Q2);nres=c(nres,"Q2")}} else Q2=NULL

    if(MAPE||MdAPE||RMSPE||RMdSPE) pe=err/y
    if(MAPE||MdAPE) ape=abs(pe)
    if(MAPE) {mape=100*mean(ape);if(rsingle) return(mape) else {res=c(res,mape);nres=c(nres,"MAPE")}} else mape=NULL
    if(MdAPE) {mdape=100*median(ape);if(rsingle) return(mdape) else {res=c(res,mdape);nres=c(nres,"MdAPE")}} else mdape=NULL

    if(RMSPE||RMdSPE) pe2=pe^2
    if(RMSPE) {rmspe=sqrt(100*mean(pe2));if(rsingle) return(rmspe) else {res=c(res,rmspe);nres=c(nres,"RMSPE")}} else rmspe=NULL
    if(RMdSPE) {rmdspe=sqrt(100*median(pe2));if(rsingle) return(rmdspe) else {res=c(res,rmdspe);nres=c(nres,"RMdSPE")}} else rmdspe=NULL

    if(SMAPE||SMdAPE) map=eabs/(abs(x)+abs(y)) 
    if(SMAPE) {smape=200*mean(map); if(rsingle) return(smape) else {res=c(res,smape);nres=c(nres,"SMAPE")}} else smape=NULL
    if(SMdAPE){smdape=200*median(map); if(rsingle) return(smdape) else {res=c(res,smdape);nres=c(nres,"SMdAPE")}} else smdape=NULL

    # same val for all: randomwalk, see Hyndman paper
    if(MRAE||MdRAE||GMRAE||THEILSU2) { if(!is.null(val)) { if(is.list(val))  val2=val[[which(metric=="MRAE"|metric=="MdRAE"|metric=="GMRAE"|metric=="THEILSU2")[1]]]
                                                           else val2=val
                                                           if(length(val2)==1) val2=c(val2,y[1:(length(y)-1)])
                                                           errb=y-val2
                                                           RESULT=TRUE
                                                         } else RESULT=FALSE
                                     } else RESULT=FALSE
    if(RESULT && (MRAE||MdRAE||GMRAE)) abs_rt=abs(err/errb)
    if(MRAE) { if(RESULT) mrae=mean(abs_rt) else mrae=NA
               if(rsingle) return(mrae) else {res=c(res,mrae);nres=c(nres,"MRAE")}
             } else mrae=NULL
    if(MdRAE){ if(RESULT) mdrae=median(abs_rt) else mdrae=NA
               if(rsingle) return(mdrae) else {res=c(res,mdrae);nres=c(nres,"MdRAE")}
             } else mdrae=NULL
    if(GMRAE){ if(RESULT) gmrae=prod(abs_rt)^(1/(length(abs_rt))) else gmrae=NA
               if(rsingle) return(gmrae) else {res=c(res,gmrae);nres=c(nres,"GMRAE")}
             } else gmrae=NULL
    if(THEILSU2){ if(RESULT) theilsu2=rmse/sqrt(mean(errb^2)) else theilsu2=NA
                  if(rsingle) return(theilsu2) else {res=c(res,theilsu2);nres=c(nres,"THEILSU2")}
                } else theilsu2=NULL

    if(MASE){ if(is.list(val)){val2=val[[which(metric=="MASE")[1]]]}
              else if(length(val)>1) val2=val else val2=NULL
              if(!is.null(val2)){mmase=mean(abs(diff(val2)));mase=mean(abs(err/mmase))} else mase=NA
              if(rsingle) return(mase) else {res=c(res,mase);nres=c(nres,"MASE")}
            } else mase=NULL

    if(SMINKOWSKI3){sminkowski3=sum(eabs^3);if(rsingle) return(sminkowski3) else {res=c(res,sminkowski3);nres=c(nres,"SMinkowski3")}} else sminkowski3=NULL
    if(MMINKOWSKI3){mminkowski3=mean(eabs^3);if(rsingle) return(mminkowski3) else {res=c(res,sminkowski3);nres=c(nres,"MMinkowski3")}} else mminkowski3=NULL 
    if(MdMINKOWSKI3){mdminkowski3=median(eabs^3);if(rsingle) return(mdminkowski3) else {res=c(res,sminkowski3);nres=c(nres,"MdMinkowski3")}} else mdminkowski3=NULL

    if(COR||q2||R2){cor=suppressWarnings(cor(y,x));if(is.na(cor)) cor=0;
                if(rsingle && COR) return(cor) else {res=c(res,cor);nres=c(nres,"COR")}
               } else cor=NULL
    if(R2) {r2=cor^2; if(rsingle) return (r2) else {res=c(res,r2);nres=c(nres,"R2")}} else r2=NULL
    if(q2){q2=1-cor^2;if(rsingle) return(q2) else {res=c(res,q2);nres=c(nres,"q2")}} else q2=NULL

    if(REC||NAREC||TOLERANCE) { if((NAREC||TOLERANCE)&& is.null(val)) val=1 
                                rec=RECcurve(y,x)
                              } else rec=NULL
    if(NAREC){ 
               if(is.list(val)) val2=val[[which(metric=="NAREC")[1]]]
               else if(length(val)>1) val2=val[which(metric=="NAREC")[1]]
               else val2=val
               if(rec[nrow(rec),1]>val2) {rec2=partialcurve(rec,val2);narec=curvearea(rec2,val2)} 
               else {val2=rec[nrow(rec),1];narec=curvearea(rec,val2)}
               if(rsingle) return(narec) else {res=c(res,narec);nres=c(nres,"NAREC")}
             } else narec=NULL
    if(TOLERANCE){ # warning last REC metric, changes rec 
                   if(is.list(val)) val2=val[[which(metric=="TOLERANCE")[1]]]
                   else if(length(val)>1) val2=val[which(metric=="TOLERANCE")[1]]
                   else val2=val
                   if(rec[nrow(rec),1]>val2) rec=partialcurve(rec,val2)
                   if(is.vector(rec)) tolerance=rec[2] else tolerance=rec[nrow(rec),2]
                   if(rsingle) return(tolerance) else {res=c(res,tolerance);nres=c(nres,"TOLERANCE")}
                 } else tolerance=NULL
   # regression return: 
   if(!is.null(res)) {names(res)=nres; 
                      # sort res:
                      I=NULL # for regression, this works perfectly?
                      for(i in 1:LM)
                        {
                          ii=which(nres==metric[i])[1]
                          if(!is.na(ii)) I=c(I,ii) 
                        }
                      res=res[I]
                     }
   if(reslist) {res=list(res=res,rec=rec)}
   return(res)
 }
}

#----------------------------------------
# RECurve by Paulo Cortez, 2006@
#
# following the article of Bi & Bennett 2003:
# J. Bi and K. Bennett, Regression Error Characteristic curve
# In Proceedings of 20th Int. Conf. on Machine Learning (ICML),  
# Washington DC, USA, 2003.
#
# vector.error - vector with the residuals or errors
#                vector.error = y (desired) - x (predicted)
#              or vector.error= y and x= predictions (2nd mode)
RECcurve=function(vector.error,x=NULL)
{
#print(vector.error)
 if(!is.null(x)) vector.error=(vector.error-x)
 Correct=0; Eprev=0; 
 ES<-sort(abs(vector.error))
#print(ES)
 M<-length(ES)+1; M1=M-1;
 X<-matrix(nrow=M,ncol=2)
 M<-length(ES); k=1;i=1;notstop=TRUE;
 while(notstop)
  { a=0; while( (i+a)<M && ES[(i+a+1)]==ES[(i+a)] ) a=a+1;
    if(a>0) {i=i+a-1; Correct=Correct+a-1;}
#cat(" >> i:",i,"a:",a,"k:",k,"prev:",Eprev,"ESi:",ES[i],"\n")
    if(Eprev<ES[i])
      { X[k,1]<-Eprev; X[k,2]<-Correct/M1; Eprev<-ES[i]; k=k+1;}
    Correct=Correct+1
    i=i+1;
    if(i>M1) notstop=FALSE;
  }
  X[k,1]<-ES[M]
  X[k,2]<-Correct/M1
#print(X)
#cat("M:",M,"k:",k,"Cor:",Correct,"\n")
  #X=na.omit(X) #X[,2]<-100*X[,2] # put in percentage
  return (X[(1:k),])
}

# ----------------------------------------------------------

# convert matrix or data.frame into factor with major class 
majorClass=function(x,L)
{
 if(is.vector(x)) return (factor(L[which.max(x)],levels=L))
 else 
 { NX=nrow(x)
   y=vector(length=NX)
   for(i in 1:NX) y[i]=L[which.max(x[i,])]
   return (factor(y,levels=L))
 }
}

# target - vector of factor with the desired values 
# predict - vector of factor with the predicted values
# D - decision thresold
# TC - target concept class, -1 not used
# note: when TC=-1 then majorclass is used, D is not considered!
Conf=function(target,pred,D=0.5,TC=-1,predreturn=FALSE)
{
 L=levels(target[1])
 if(is.vector(pred)) # numeric predictions equal to classes 
 { if(length(L)>2) pred=factor(pred,levels=L)
   else  { 
           if(TC==1) LB=c("TRUE","FALSE") else LB=c("FALSE","TRUE")
           pred=factor(pred>D,levels=LB); target=factor((target==L[TC]),levels=LB)
         }
 }
 else if(is.factor(pred)) # pred factor
 {
   LP=levels(pred[1])
   if(length(LP)<length(L)) levels(pred)=L # only if pred has less classes than target
   if(TC>0)
     {
      pred=factor(pred==L[TC],levels=c("FALSE","TRUE"))
      target=factor((target==L[TC]),levels=c("FALSE","TRUE"))
     }
 }
 else if(!is.factor(pred)) #if(is.matrix(pred) || is.data.frame(pred))  # probabilities
  { if(TC>0) { pred=factor(pred[,TC]>D,levels=c("FALSE","TRUE")); target=factor((target==L[TC]),levels=c("FALSE","TRUE"));}
    else pred=majorClass(pred,L)
  }
 if(predreturn) return(list(conf=table(target,pred),pred=pred)) else return(table(target,pred))
}

# ------------------------------------------------------------------
LIFTcurve<-function(y,x,TC)
{
 if(TC>0 && !is.vector(x) ) { x=x[,TC];} else TC=2
 NC=NCOL(x)
 if(is.factor(y[1])) POS=levels(y[1])[TC] else POS=1
 if(NC==2) x=x[,2]
 NR=NROW(x); if(NR>100) NR=100
 alift=twoclassLift(y,x,Positive=POS,STEPS=NR,type=3)
 return (list(alift=alift,area=curvearea(alift,1)))
}

# ------------------------------------------------------------------
# call of the ROC function: 
# - calls multiROCcurve: if x is matrix or data.frame!
# - calls ROCcurve: else.
# y - vector of factor or numeric (0,1) with the desired values 
# x - vector or matriz of numeric with the predicted values
# TC - target class
ROCcurve<-function(y,x,TC=-1) #,method="int")
{
 if(TC>0 && !is.vector(x) ) { x=x[,TC];} else TC=2
 C=NCOL(x)
 if(C>2) #return (multiROC(y,x)) #,method=method))
 {
  # Provost and Domingos AUC formulation for Multiclass problems
  # y - vector of factor with the desired values 
  # x - matriz of numeric with the predicted values
  #     Note: the sum of x[i,] should be 1 for all i!!!
  ROC=vector("list",C)
  # prevalence of each class:
  SUM=length(y)
  Lev=levels(y[1])
  p=table(y)[]/SUM
  aux=0.0
  for(i in 1:C)
   { #print(paste("i:",i))
     R=twoclassROC(y,x[,i],Positive=Lev[i]) #,method=method)
     ROC[[i]]=R
     if(p[i]>0) aux=aux+R$auc*p[i]
     #cat("i:",i,"p:",p[i],"auc:",R$auc,"aux:",aux,"\n")
   }
  #ROC=list(roc=ROC,auc=aux)
  return(list(roc=ROC,auc=aux))
  #return (ROC) # use: ob$roc[[i]]$roc or ob$roc[[i]]$auc to access individual rocs for each class i
 }
 else{
       if(is.factor(y[1])) POS=levels(y[1])[TC] else POS=1
       if(C==2) x=x[,2]
       #cat(" >> POS:",POS,"lev:",levels(y),"\n")
       return (twoclassROC(y,x,Positive=POS)) #,method=method))
     }
}

# ------------------------------------------------------------------
# practical efficient method for ROC and AUC value
# algorithm 2 of Fawcett 2003, algorithm 1 of Fawcett 2006
# notes: use only with 2 classes
#        this is 2nd implementation, where the <-c(,) was replaced
#        by a much faster [,1]<- and [,2]<- instructions 
#
# y - vector of factor/numeric with the desired values 
# x - vector of numeric with the predicted values
# Positive - a label or number that corresponds to a TRUE/positive class value
## method = "int" - interpolate between 2 points, "pes" - pessimistic Fawcett point, "opt" - optimistic Fawcett point
# ------------------------------------------------------------------
twoclassROC<-function(y, x, Positive=1) #, method="int")
{
 #print(method)
# if(is.factor(y)) {y=as.numeric(y)-1;Positive=1;}
#cat("Y:\n")
#print(summary(y))
#cat("X:\n")
#YY<<-y; XX<<-x
#print(summary(x))
#cat("<<< Posititive:",Positive,"\n")

  Xsize<-length(y)
  Pos<-sum(y[]==Positive) # total actual positives
#PP<<-Positive
#cat("Pos:",Pos,"\n")
  Neg<-Xsize-Pos          # total actual negatives

  Ord<-order(x,decreasing=TRUE) # very fast sort of vector
 
  FP<- 0
  TP<- 0
  FPprev<- 0
  TPprev<- 0
  A<-0  
  fprev<- -Inf

#cat(" --- AUC:",A,"\n")
  R<-matrix(ncol=2,nrow=(Xsize+1))
  k<-1
  for(i in 1:Xsize) 
     {
       if (x[Ord[i]]!=fprev)  
            { 
              if(FP>0) R[k,1]<-FP/Neg else R[k,1]=0
              if(TP>0) R[k,2]<-TP/Pos else R[k,2]=0 # the ROC point
#cat("k:",k,"FP:",FP,"TP:",TP,"Neg:",Neg,"Pos:",Pos,"\n")
              if(!is.na(R[k,1])) 
              {k<-k+1
               A<-A+trap_area(FP,FPprev,TP,TPprev) #,method) # compute the AUC
               fprev <- x[Ord[i]]
               FPprev<- FP
               TPprev<- TP
              }
            }
       if (y[Ord[i]]==Positive) TP<-TP+1 # test[i]
       else FP<-FP+1 
     }
  if(FP>0) R[k,1]<-FP/Neg else R[k,1]=0
  if(TP>0) R[k,2]<-TP/Pos else R[k,2]=1 # the ROC point
#cat("k:",k,"FP:",FP,"TP:",TP,"Neg:",Neg,"Pos:",Pos,"\n")
  if(FP==0 && k<(Xsize+1)) {k=k+1; R[k,]=c(1,1)}
#cat(" --- AUC:",A,"pos:",Pos,"neg:",Neg,"TP",TP,"FP",FP,"\n")
#cat(" ---: FPprev:",FPprev,"TPprev:",TPprev,"trap:",trap_area(1,FPprev,Pos,Pos),"\n")
  if(Neg>0) A<-A+trap_area(Neg,FPprev,Pos,TPprev) #,method)
  else A<-A+trap_area(1,FPprev,Pos,Pos) #,method)
  if(Pos>0 && Neg>0) A<-A/((1.0*Pos)*Neg) 
  else if(Neg>0) A<-A/((1.0*Neg))
  else A<-A/((1.0*Pos))

#RR<<-R;RR=na.omit(RR);print(RR[])
 
  ROC<-list(roc=R[(1:k),],auc=A)
#cat(" --- AUC:",A,"\n")
  return (ROC)
}
# ------------------------------------------------------------------
# internal R function used by ROCcurve: do not use this
# ------------------------------------------------------------------
trap_area<-function(X1,X2,Y1,Y2) #,method="int")
{ 
  return ( (abs(X1-X2)) * ((Y1+Y2)/2) )
  #if(method=="int") return ( (abs(X1-X2)) * ((Y1+Y2)/2) )
  #else if(method=="opt") return ( (abs(X1-X2)) * ((Y1+Y1)/2) ) 
  #else if(method=="pes") return ( (abs(X1-X2)) * ((Y2+Y2)/2) ) 
}
# ------------------------------------------------------------------
xmiddle_point<-function(X1,X2,Y1,Y2,X3)
{ 
 m=(Y2-Y1)/(X2-X1); b=Y1-m*X1;
 return (m*X3+b)
}
#-------------------------------------------------------------------
# vertical averaging of ROC curves, algorithm 3 from Fawcett 2006
# samples - number of FP samples
# ROCS list with length(ROCS) ROC curves, each ROC is [,1] frp and [,2] tpr
#
# returns tpravg with samples+1 rows and 3+nrocs columns: fpr, tpr, mean, confint95, tpr_roc[[1]],...,tpr_roc[[nrocs]]
vaveraging<-function(samples,ROCS,min=0,max=1)
{
  #SAMPLES<<-samples;ROCS<<-ROCS;MIN<<-min;MAX<<-max;
  s=1
  nrocs=length(ROCS)
  if(is.character(max)) 
  {
    fprsamples=1:length(ROCS[[1]][,1])
    for(i in 1:nrocs)
      {
        ROCS[[i]][,1]=fprsamples
        ROCS[[i]]=apply(ROCS[[1]][],2,as.numeric)
      }
  }
  else fprsamples=seq(min,max,length.out=samples)


#cat("min:",min,"max:",max,"\n")
  tpravg=matrix(ncol=(3+nrocs),nrow=length(fprsamples))
  for(k in fprsamples)
  {
    #tprsum=0
    tprsum=rep(0,nrocs)
    for(i in 1:nrocs)
    {
#cat("k:",k,"i:",i,"\n")
     #tprsum=tprsum+TPR_FOR_FPR(k,ROCS[[i]],nrow(ROCS[[i]]))
     tprsum[i]=tprsum[i]+TPR_FOR_FPR(k,ROCS[[i]],nrow(ROCS[[i]]))
#cat("tprsum[",i,"]=",tprsum[i],"\n")
    }
    #tpravg[s,]=c(k,tprsum/nrocs)
#cat("conf:\n")
#TPR<<-tprsum
    tpravg[s,]=c(k,mean(tprsum),conflevel(tprsum),tprsum)
#cat("conf done\n")
    s=s+1
  } 
  return(tpravg)
}
# internal R functions used by vaveraging: do not use this
TPR_FOR_FPR<-function(fprsample,ROC,npts)
{
 # error here, think later...
 #RRR<<-ROC
 #NPTS<<-npts
 #FPR<<-fprsample
 i=1
 NR=nrow(ROC)
 while(i<npts && ROC[(i+1),1]<=fprsample) i=i+1;
 if(i<NR && ROC[i,1]==fprsample) return (ROC[i,2])
 else if(i<NR) return (INTERPOLATE(ROC[i,],ROC[(i+1),],fprsample))
 else return (ROC[i,2])
# else return (INTERPOLATE(ROC[(NR-1),],ROC[NR,],fprsample))
}
INTERPOLATE<-function(ROCP1,ROCP2,X)
{
 #cat("rocp1:",ROCP1,"rocp2:",ROCP2,"x:",X,"\n")
 slope=(ROCP2[2]-ROCP1[2])/(ROCP2[1]-ROCP1[1])
 return (ROCP1[2]+slope*(X-ROCP1[1]))
}
# 95% confidence interval according to a t-student distribution
conflevel=function(x,level=0.95)
{
 RES=try( (t.test(x,conf.level=level)$conf[2]-t.test(x,conf.level=level)$conf[1])/2 , silent=TRUE)
 if(class(RES)=="numeric") return(RES) else return (0)
}

# partial curve (roc, rec, ...)
partialcurve=function(Curve,threshold=1) #,method="int")
{
 I=which(Curve[,1]<=threshold)
 IND=I[length(I)]
 if(Curve[IND,1]==threshold) M=Curve[(1:IND),]
 else 
 {
   M=Curve[(1:(IND+1)),]
   M[(IND+1),]=c(threshold,xmiddle_point(Curve[IND,1],Curve[(IND+1),1],Curve[IND,2],Curve[(IND+1),2],threshold))
   ##if(method=="int") M[(IND+1),]=c(threshold,xmiddle_point(Curve[IND,1],Curve[(IND+1),1],Curve[IND,2],Curve[(IND+1),2],threshold))
   ##else if(method=="pes") M[(IND+1),]=c(threshold,Curve[IND,2])
   ##else if(method=="opt") M[(IND+1),]=c(threshold,Curve[(IND+1),2])
 }
 return(M)
 #return(list(curve=M,auc=rocarea(M,threshold),tprfpr=M[(IND+1),2]))
}

# area of a curve using trapesoidal method
# examples: auc of a roc curve, area of rec, etc...
curvearea<-function(Curve,threshold=1.0) #,method="int")
{ 
  if(is.vector(Curve)) return (0)
  else 
  { if(Curve[nrow(Curve),1]>threshold) Curve=partialcurve(Curve,threshold) #,method=method)
    A=0
    for(i in 2:nrow(Curve)) 
       {
        A=A+trap_area(Curve[i,1],Curve[(i-1),1],Curve[i,2],Curve[(i-1),2]) #,method=method)
#cat("A",A,"T",threshold,"\n")
       }
#cat("A",A,"A/T",A/threshold,"\n")
    if(threshold>0) return (A/threshold)
    else return (0)
  }
}

#-------------------------------------------------------------------

# tolerance of a rec curve
tolerance<-function(REC,tol=0.5)
{
 stop=FALSE; i=1;N=nrow(REC)
 while(i<N && REC[i,1]< tol ) {i=i+1;}

 if(i==N || REC[i,1]==tol) return (REC[i,2])
 else if(i>1) return ( xmiddle_point(REC[(i-1),1],REC[i,1],REC[(i-1),2],REC[i,2],tol) )
}

# mean and confidence interval using t.test
meanint<-function(x,level=0.95)
{
 if(is.matrix(x)||is.data.frame(x))
 {
  C=ncol(x); M=rep(0,C); Conf=rep(0,C);
  for(i in 1:C)
  { M[i]=mean(x[,i])
    Conf[i]=conflevel(x[,i],level=level)
  }
 }
 else
 {
  M=mean(x)
  Conf=conflevel(x,level=level)
 }
 return(list(mean=M,int=Conf))
}
# --------

# m - is matrix or data.frame
#mpairwise=function(m,p.adj="bonf",paired=TRUE)
#{
# NC=NCOL(m);NR=NROW(m)
# x=vector(length=NC*NR); g=x;
# for(i in 1:NC) { ini=(i-1)*NR+1;end=ini+NR-1;
#                  x[ini:end]=m[,i];g[ini:end]=rep(i,NR);}
# g=factor(g)
# P=pairwise.t.test(x,g,p.adjust.method =p.adj,paired=paired)
# return(P)
#}


# experimental stuff:
# type 1 - normal lift
# type 2 - cumulative
# type 3 - cumulative percentage
twoclassLift<-function(y, x, Positive=1,STEPS=10,type=3)
{
#YY<<-y;XX<<-x;PP<<-Positive;
#y=YY;x=XX;Positive=PP;STEPS=10;type=3
  #STEPS=11;Positive="setosa";type=3
  DIV=STEPS^2
  Xsize<-length(y)
  APos<-sum(y[]==Positive) # total actual positives
  if(is.data.frame(x)) x=x[,1] 
  Ord<-order(x,decreasing=TRUE) # very fast sort of vector
  ALL=APos/Xsize
  if(type==3){R=matrix(ncol=2,nrow=STEPS+1)
              R[1,]=c(0,0)
             }
  else R=matrix(ncol=2,nrow=STEPS)
  Portion=Xsize/STEPS
  for(i in 1:STEPS) 
     {
       N=round(i*Portion)
       if(type==1 || type==3) IND=1:N # total actual positives
       else if(type==2) IND=(N-STEPS+1):N
       Pos=sum(y[Ord[IND]]==Positive)
       if(type==3) { R[(i+1),1]=i*STEPS/DIV; R[(i+1),2]=Pos/APos }
       else {R[i,1]=i*STEPS/DIV; R[i,2]=Pos/(ALL*N)}
     }
  return (R)
}

# yaggregate
yaggregate=function(y,N=1)
{if(N==1) return(mean(y)) 
 else if(N==3){r=range(y);return(c(r[1],mean(y),r[2]));} 
 else{y=quantile(y,seq(0,1,length.out=N));attr(y,"names")=NULL;return(y)}
}

# ---- I do not use these functions anymore, kept here just for backup purposes:
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if(FALSE){
#---------------------------------------------------------
# Brier Score: computes the brier score (SSE) for probabilitic output models ([0,1]) 
# y - target factor
# x - numeric predictions (in 0 to 1 probabilities) 
#---------------------------------------------------------
Brier=function(y,x,TC=-1)
{
  L=levels(y[1]);N=length(y); 
  if(TC>0){ NL=1;x=x[,TC];L=L[TC];MSE=1;} else{NL=length(L); PROP=table(y)[]/N; TBRIER=0;MSE=rep(FALSE,(NL+1));} 
  for(i in 1:NL) 
  {
      T=as.numeric(y==L[i])
      if(is.vector(x)) MSE[i]=Sse(T,x)/N 
      else { MSE[i]=Sse(T,x[,i])/N # Brier=MSE
             if(PROP[i]>0) TBRIER=TBRIER+PROP[i]*MSE[i]
           }
  }
  if(!is.vector(x)) MSE[(NL+1)]=TBRIER
  return (MSE)
}

# new: quantile functions:
q1=function(x){return(quantile(x,0.25)[])}
q3=function(x){return(quantile(x,0.75)[])}

# several metric functions (for fast access)
# sum of squared errors
# x - vector of predictions, y - vector of desired values  
#--- 2013
# classification: 
# y - vector of factor with the desired values 
# x - vector of factor with the predicted values
# D - decision thresold
# TC - target concept class, -1 not used

# ------------------------------------------------------------------
# Provost and Domingos AUC formulation for Multiclass problems
# y - vector of factor with the desired values 
# x - matriz of numeric with the predicted values
#     Note: the sum of x[i,] should be 1 for all i!!!
multiROC<-function(y,x) #,method="int")
{
 C=NCOL(x) # number of classes
 ROC<-vector("list",C)
 # prevalence of each class:
 SUM<-length(y)
 Lev<-levels(y[1])
 p=table(y)[]/SUM
 aux=0.0
 for(i in 1:C)
   { #print(paste("i:",i))
     R<-twoclassROC(y,x[,i],Positive=Lev[i]) #,method=method)
     ROC[[i]]=R
     if(p[i]>0) aux=aux+R$auc*p[i]
     #cat("i:",i,"p:",p[i],"auc:",R$auc,"aux:",aux,"\n")
   }
  ROC<-list(roc=ROC,auc=aux)
  return (ROC) # use: ob$roc[[i]]$roc or ob$roc[[i]]$auc to access individual rocs for each class i
}


Accuracy<-function(y,x,D=0.5,TC=-1)
{
 conf=Conf(y,x,TC=TC,D=D)
 D=0; for(i in 1:NCOL(conf))D=D+conf[i,i]
 return (100*D/length(y))
}

Kappa=function(y,x,D=0.5,TC=-1)
{
 conf=Conf(y,x,D,TC);Total=sum(conf);Diag=0;DiagR=0
 for(i in 1:NCOL(conf)) 
    {
      Diag<-Diag+conf[i,i]
      DiagR<-DiagR+(sum(conf[i,])*(sum(conf[,i])/Total))
    }
 return (100*(Diag-DiagR)/(Total-DiagR))
}

# future: compare CRAMERV to flat => variability measure?
# check chisq.test errors... => must have at least 2 levels, etc...
CRAMERV=function(target,pred,D=0.5,TC=-1)
{
 L=levels(target[1])
 if(is.vector(pred)) 
 { if(length(L)>2) pred=factor(pred,levels=L)
   else  { 
           if(TC==1) LB=c("TRUE","FALSE") else LB=c("FALSE","TRUE")
           pred=factor(pred>D,levels=LB); target=factor((target==L[TC]),levels=LB)
         }
 }
 else if(is.factor(pred) && TC>0) 
 {
   pred=factor(pred==L[TC],levels=c("FALSE","TRUE"))
   target=factor((target==L[TC]),levels=c("FALSE","TRUE"))
 }
 else if(!is.factor(pred)) #if(is.matrix(pred) || is.data.frame(pred))  # probabilities
  { if(TC>0) { pred=factor(pred[,TC]>D,levels=c("FALSE","TRUE")); target=factor((target==L[TC]),levels=c("FALSE","TRUE"));}
    else pred=majorClass(pred,L)
  }

 T=chisq.test(target,pred)
 T=T$statistic[[1]]
 k=length(levels(target[1]))
 N=length(target)
 return (sqrt(T/(N*(k-1))))
}



# new rminer 1.3, embrechts hint, http://www.nipsfsc.ecs.soton.ac.uk/evaluation/
# Balanced Error Rate (BER)
metricBER=function(y,x,D=0.5,TC=-1)
{
 conf=Conf(y,x,TC=TC,D=D)
 NC=NCOL(conf)
 RES=vector(length=NC)
 for(i in 1:NC) { 
                  RES[i]=sum(conf[i,-i])/sum(conf[i,])
                }
 return (100*sum(RES)/NC)
}

f1_score=function(y,x,D,TC)
{
 M=metrics(y,x,D)
 Precision=M$precision[TC]
 Recall=M$tpr[TC]
 f1=rep(0,length(TC))
 if(length(TC)>1) I=TC else I=1
 for(i in I) if(Precision[i]!=0 && Recall[i]!=0) f1[i]=2*((Precision[i]*Recall[i])/(Precision[i]+Recall[i]))
 return(f1)
}

Tbrier=function(y,x,TC=-1){ B=Brier(y,x,TC=TC);return(B[length(B)])}

SMinkowski=function(y,x,q=1){ return (sum( abs(y-x)^q )) } # minkowski loss function, Bishop 2006, pattern recognition...
MMinkowski=function(y,x,q=1){ return (mean( abs(y-x)^q )) } 
MdMinkowski=function(y,x,q=1){ return (median( abs(y-x)^q )) }

Me=function(y,x) { return (sum(y-x)) }
Sse=function(y,x) { return (sum((y-x)^2)) }
Mse=function(y,x) { return (mean((y-x)^2)) }
Mdse=function(y,x) { return (median((y-x)^2)) }
Gmse=function(y,x){ return (Gmean((y-x)^2)) }
Rse=function(y,x,ymean=mean(y)) { return (100* sum((y-x)^2)/sum((y-ymean)^2)) } # relative squared error
Rmse=function(y,x) { return (sqrt(mean((y-x)^2))) }
Rrse=function(y,x) { return ( 100*sqrt(sum((y-x)^2)/sum((y-mean(y))^2)) ) }
Sae=function(y,x) { return (sum(abs(y-x))) }
Mae=function(y,x) { return (mean(abs(y-x)))}
Gmae=function(y,x){ return (Gmean(abs(y-x)))}
Rae=function(y,x,ymean=mean(y)) { return (100*sum(abs(y-x))/sum(abs(y-ymean))) } # also known as CumRAE, RelMAE
#Smape=function(y,x) { return ( 100*mean(abs(y-x)/(x+y)))} # wikipedia def.

# forecasting specific, R. Hyndman IJF 2006 "Another Look at Measures of Forecast Accuracy" definition:
Mdae=function(y,x) { return (median(abs(y-x)))}
Mape=function(y,x) { return (100*mean(abs((y-x)/y)))} # R. Hyndman IJF 2006 def:
Mdape=function(y,x) { return (100*median(abs((y-x)/y)))}
Rmspe=function(y,x) { return (sqrt(100*mean(((y-x)/y)^2)))}
Rmdspe=function(y,x) { return (sqrt(100*median(((y-x)/y)^2)))}
#Smape=function(y,x) { return (100*mean(abs(y-x)/((abs(x)+abs(y))/2)))} # def. of http://www.neural-forecasting-competition.com/motivation.htm
Smape=function(y,x) { return (200*mean(abs(y-x)/(abs(x)+abs(y))))}
Smdape=function(y,x) { return (200*median(abs(y-x)/(abs(x)+abs(y))))}

# relative errors:
# b= val2mark naive method forecasts.

# random walk, with or without drift
# ts - time series in samples
# H - number of forecasts, horizon
randomwalk=function(ts,H,drift=TRUE)
{ if(drift) drift=mean(diff(ts)) else drift=0
  return (rep(ts[length(ts)],H)+1:H*drift)
}
Gmean=function(x){return(prod(x)^(1/(length(x))))} # auxiliar geometric mean

Mrae=function(y,x,ts,b=randomwalk(ts,length(y))){if(is.null(ts)&&is.null(b)) return(NA) else return(mean(abs((y-x)/(y-b))))}
Mdrae=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (median(abs((y-x)/(y-b)))) }
Gmrae=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (Gmean(abs((y-x)/(y-b)))) }
TheilsU2=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (Rmse(y,x)/Rmse(y,b)) } # theils'U or U2
#Mase=function(y,x,ts) { N=length(ts); K=1/(N-1); SUM=sum(abs(ts[2:N]-ts[1:(N-1)])); return ( mean( abs( (y-x)/(K*SUM) )) ) }
# ts - time series in samples 
Mase=function(y,x,ts) { N=length(ts); MEAN=mean(abs(ts[2:N]-ts[1:(N-1)])); return ( mean(abs((y-x)/MEAN)) ) } # faster variant?

# wikipedia:
#coefficient of determination, most general definition:
pressR2=function(y,x,ymean=mean(y)) { return (1-sum((y-x)^2)/sum((y-ymean)^2)) }
pressR22=function(y,x,ymean=mean(y)) { return (sum((x-mean(x))^2)/sum((y-ymean)^2)) }
Ss=function(y){return (sum((y-mean(y))^2))} # sum of squares

# NEW rminer 1.3?
# mark embrechts metrics:
Test_q2=function(y,x){1-Correlation(y,x)^2}
Test_Q2=function(y,x,ymean=mean(y)) {return (sum((y-x)^2)/sum((y-ymean)^2))}

Correlation=function(y,x) { COR=suppressWarnings(cor(x,y)); if(is.na(COR)) COR=0; return(COR)}
Narec=function(y,x,val=1){ if(is.null(val)) val=1; R=RECcurve(y,x);if(R[nrow(R),1]>val)R=partialcurve(R,val) else val=R[nrow(R),1];return(curvearea(R,val))}
Tolerance=function(y,x,val=1){ if(is.null(val)) val=1; R=RECcurve(y,x);if(R[nrow(R),1]>val)R=partialcurve(R,val);if(is.vector(R)) return (R[2]) else return(R[nrow(R),2])}
Nauc=function(y,x,val=1,TC=-1){ if(is.null(val)) val=1; RR=ROCcurve(y,x,TC=TC); RR=partialcurve(RR$roc,val); return(curvearea(RR,val))}
Nalift=function(y,x,val=1,TC=-1){ if(is.null(val)) val=1; RR=LIFTcurve(y,x,TC=TC); RR=partialcurve(RR$alift,val); return(curvearea(RR,val))}
Tprfpr=function(y,x,val,TC=-1){ if(is.null(val)) val=0.01; RR=ROCcurve(y,x,TC=TC); RR=partialcurve(RR$roc,val);if(is.vector(RR)) return (RR[2]) else return(RR[nrow(RR),2])}
Aliftperc=function(y,x,val,TC=-1){ if(is.null(val)) val=0.1; RR=LIFTcurve(y,x,TC=TC); RR=partialcurve(RR$alift,val);if(is.vector(RR)) return (RR[2]) else return(RR[nrow(RR),2])}
}
# end of:---- I do not use these functions anymore, kept here just for backup purposes -----
