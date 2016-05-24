setClass("SSM",representation(Info="list",Bootstrap="list",show="matrix"))
setMethod("show", "SSM",
          function(object) {
            cat(paste("\n------------------------------------------"))
            cat(paste("\n-          Superior Set of Models        -"))
            cat(paste("\n------------------------------------------\n"))
            show(object@show)
            cat(paste("\nDetails"))
            cat(paste("\n------------------------------------------\n"))
            cat(paste("\nNumber of eliminated models\t:\t"))
            cat(object@Info$n_elim)
            cat(paste("\nStatistic\t:\t"))
            cat(object@Info$statistic)
            cat(paste("\nElapsed Time\t:\t"))
            print(object@Info$elapsed.time)
          }
)
boot.block<-function(x,v,n,k){
  startIndexes=sample(1:(n-k),v+1)
  blocks=do.call(c,lapply(startIndexes,function(p,k) p+seq_len(k),k=k ))
  return(blocks[1:n])
}
d_b_i_mean.fun<-function(x,M,model.names){      
  export=sapply(1:M,function(l) sum(x[paste(model.names[l],model.names[-l],sep=".")])/(M-1))
  return(export)
}  
LossVaR<-function(realized, evaluated, which='asymmetricLoss', type='normal', delta=25,tau){
  if(delta<=0) stop("delta must be greater than zero")
  if(which!='asymmetricLoss') stop('Not valid choice of "which". Valid choices are "asymmetricLoss"')
  if(which=='asymmetricLoss' & all(type!=c('normal','differentiable'))) stop('Not valid choice of "type". Valid choices are "normal" and "differentiable" ')
  
  if(tau<0 | tau>1) stop('tau must be in (0,1)')
  
  if(is.matrix(evaluated)){m=ncol(evaluated)}else{m=1}
  if(!is.matrix(evaluated)){evaluated=as.matrix(evaluated)}
  if(!is.matrix(realized)){realized=as.matrix(realized)}
  if(ncol(realized)>1) stop('realized must of one column')
  if(nrow(realized)!=nrow(evaluated)) stop('the number of realized and evaluated observation must be the same')
  
  if(which=='asymmetricLoss'){
    if(type=='normal') {
      dummies=matrix(0,nrow=length(realized),ncol=m)
      for(j in 1:m) dummies[realized<evaluated[,j],j]=1
      loss=(tau-dummies)*(realized-evaluated)
    }
    if(type=="differentiable"){
      loss=(tau-(1+exp(delta*(realized-evaluated)))^-1 )*(realized-evaluated)
    }
  }
  return(loss)
}
LossVol<-function(realized, evaluated, which='SE1'){
  
  if(all(which!=c('SE1','SE2','QLIKE','R2LOG','AE1','AE2'))) stop('Not valid choice of "which". Valid choices are "SE1","SE2","QLIKE","R2LOG","AE1","AE2" ')
  
  if(is.matrix(evaluated)){m=ncol(evaluated)}else{m=1}
  if(!is.matrix(evaluated)){evaluated=as.matrix(evaluated)}
  if(!is.matrix(realized)){realized=as.matrix(realized)}
  
  if(ncol(realized)>1) stop('realized must of one column')
  if(nrow(realized)!=nrow(evaluated)) stop('the number of realized and evaluated observation must be the same')
    
  if(which=="SE1") loss=(evaluated-realized)^2
  if(which=="SE2") loss=(evaluated^2-realized^2)^2
  if(which=="QLIKE") loss=log(evaluated^2)+realized^2*evaluated^-2
  if(which=="R2LOG") loss=(log(realized^2*evaluated^-2))^2
  if(which=="AE1") loss=abs(evaluated-realized)
  if(which=="AE2") loss=abs(evaluated^2-realized^2)
  return(loss)
}
LossLevel<-function(realized, evaluated, which='SE'){
  
  if(all(which!=c('SE','AE'))) stop('Not valid choice of "which". Valid choices are "SE","AE" ')
    
  if(is.matrix(evaluated)){m=ncol(evaluated)}else{m=1}
  if(!is.matrix(evaluated)){evaluated=as.matrix(evaluated)}
  if(!is.matrix(realized)){realized=as.matrix(realized)}
  if(ncol(realized)>1) stop('realized must of one column')
  if(nrow(realized)!=nrow(evaluated)) stop('the number of realized and evaluated observation must be the same')
  
  if(which=="SE") loss=(evaluated-realized)^2
  if(which=="AE") loss=abs(evaluated-realized)
  return(loss)
}
MCSprocedure<-function(Loss,alpha=0.15,B=5000,cl=NULL,ram.allocation=T,statistic="Tmax",k=NULL,min.k=3,verbose=T){
  
  time.start=Sys.time()
  
  Loss=as.matrix(Loss)
  
  M_start=ncol(Loss)
  
  #no points in colnames(Loss)
  colnames(Loss)=gsub(".","_",colnames(Loss),fixed=T)
  if(is.null(colnames(Loss))) colnames(Loss)=paste("model",1:ncol(Loss),sep="_")
  if(is.null(cl)) cl_number=0
  if(!is.null(cl)) max.cores=length(cl)
  
  B=round(B)
  if(B<1000) cat(paste("Warning: B is small"))
  
  if(any(is.na(Loss))) stop("NAs in Loss are not allowed")
  if(any(abs(Loss)==Inf)) stop("Inf in Loss are not allowed")
  
  
  repeat{
    
    
    M=ncol(Loss)
    model.names=colnames(Loss)
    
    col.names.d=do.call(c,lapply(1:M,function(x) paste(model.names[x],model.names[-x],sep=".")))
    
    d=do.call(cbind,lapply(1:M,function(x) Loss[,x]-Loss[,-x]))
    
    colnames(d)=col.names.d
    
    d_ij_mean=colMeans(d)
    
    d_i_mean=sapply(1:M,function(x) sum(d_ij_mean[paste(model.names[x],model.names[-x],sep=".")])/(M-1))
    names(d_i_mean)=model.names
    
    block.length=NULL
    
    foo=expand.grid(1:M,1:M)[,2:1]
    foo=foo[foo[,1]!=foo[,2],]
    
    index=col.names.d[foo[,1]<foo[,2]]
    
    if(ncol(d)>2){
      
      d_cut=as.list(as.data.frame(d[,index]))
      
      if(!is.null(cl)) k=max(na.omit(as.numeric(parSapply(cl,d_cut,function(x) try(ar(x)$order,silent=T)))))
      
      if(is.null(cl)) k=max(na.omit(as.numeric(sapply(d_cut,function(x) try(ar(x)$order,silent=T)))))
      
    }else{
      k=ar(d[,1])$order
    }
    
    
    if(k<min.k){k=min.k}
    
    v=ceiling(nrow(d)/k)
    
    n=nrow(d)
    
    if(!is.null(cl)) indexes_b=parLapply(cl,1:B,boot.block,v=v,n=n,k=k)
    if(is.null(cl)) indexes_b=lapply(1:B,boot.block,v=v,n=n,k=k)
    
    if(!is.null(cl)){
      if(ram.allocation){
        weigth.in.cl=(as.numeric(object.size(d))/1e6)*3.2
        
        used.memory=sum(sapply(ls(),function(x) object.size(get(x))))/1e6
        
        cl_number=min(as.numeric(ceiling((memory.limit()*0.8-used.memory)/weigth.in.cl)-1),max.cores)
      }else{cl_number=max.cores}
      
    }else{
      cl_number=1
    }
    
    d=as.data.frame(d)
    
    if(cl_number>1){
      d_ij_avg_resampled=parLapply(cl[1:cl_number],indexes_b,function(x,d,B){ colSums((d[x,]))/B},d=d,B=B)
    }else{        
      d_ij_avg_resampled=lapply(indexes_b,function(x,d,B){ colSums((d[x,]))/B},d=d,B=B)
    }
    
    if(!is.null(cl)) d_b_i_mean=do.call(rbind,parLapply(cl,d_ij_avg_resampled,d_b_i_mean.fun,M=M,model.names=model.names))   
    if(is.null(cl)) d_b_i_mean=do.call(rbind,lapply(d_ij_avg_resampled,d_b_i_mean.fun,M=M,model.names=model.names))   
    
    d_b_i_var=colSums(((d_b_i_mean-d_i_mean)^2)/B)
    names(d_b_i_var)=model.names
    
    d_ij_avg_resampled=do.call(rbind,d_ij_avg_resampled)
    
    d_var_ij=matrix(colSums(((d_ij_avg_resampled-d_ij_mean)^2))/B,nrow=1)
    colnames(d_var_ij)=names(d_ij_mean)
    
    TR=max((abs(d_ij_mean))/(d_var_ij^0.5))
    
    TM=max(d_i_mean/(d_b_i_var^0.5))
    
    Tb_R=sapply(1:B,function(i) max(abs(d_ij_avg_resampled[i,]-d_ij_mean)/(d_var_ij^0.5)))
    
    Tb_M=sapply(1:B,function(i) max((d_b_i_mean[i,]-d_i_mean)/(d_b_i_var^0.5)))
    
    Pr=length(which(TR<Tb_R))/B
    
    Pm=length(which(TM<Tb_M))/B
    
    #definisco le tre elimination rule
    
    v_i_M=d_i_mean/(d_b_i_var^0.5)
    v_i_M_order=order(v_i_M,decreasing=T)
    
    
    TR_all=((d_ij_mean))/(d_var_ij^0.5)
    
    v_i_R=sapply(1:M,function(x) max(TR_all[1,paste(model.names[x],model.names[-x],sep=".")]))
    names(v_i_R)=model.names
    v_i_R_order=order(v_i_R,decreasing=T)
    
    Pm_h0_mk=NULL
    
    for (i in 1:M) {
      if(i==1) {model_temp=model.names}else{model_temp=model.names[-v_i_M_order[1:(i-1)]]}
      TM_temp=max(d_i_mean[model_temp]/(d_b_i_var[model_temp]^0.5))
      Pm_h0_mk=c(Pm_h0_mk,length(which(TM_temp<Tb_M))/B)
    }
    
    
    mcs_Pm=sapply(1:M,function(x) max(Pm_h0_mk[1:x]))
    
    #For TR statistic
    
    combine.names=names(d_ij_mean)
    
    Pr_h0_mk=NULL
    
    for (i in 1:M) {
      remove=NULL
      if(i==1) {model_temp=combine.names}else{
        remove=do.call(c,lapply(1:(i-1),function(x)  which(gsub(model.names[v_i_R_order[x]],"",combine.names)!=combine.names)))
        model_temp=combine.names[-remove]
      }
      if(i<M){
        TR_temp=max((abs(d_ij_mean[model_temp]))/(d_var_ij[1,model_temp]^0.5))
        Pr_h0_mk=c(Pr_h0_mk,length(which(TR_temp<Tb_R))/B)}else{
          Pr_h0_mk=c(Pr_h0_mk,1)
        }
    }
    
    
    mcs_Pr=sapply(1:M,function(x) max(Pr_h0_mk[1:x]))
    
    matrix_show=matrix(,M,7,dimnames=list(model.names,c("Rank_M","v_M","MCS_M","Rank_R","v_R","MCS_R","Loss")))
    
    
    matrix_show[,"v_M"]=v_i_M
    matrix_show[names(sort(v_i_M)),"Rank_M"]=1:M
    matrix_show[model.names[v_i_M_order],"MCS_M"]=mcs_Pm
    matrix_show[,"v_R"]=v_i_R
    matrix_show[names(sort(v_i_R)),"Rank_R"]=1:M
    matrix_show[model.names[v_i_R_order],"MCS_R"]=mcs_Pr
    matrix_show[,"Loss"]=colMeans(Loss)
    
    rm(list="indexes_b")
    
    if(!is.null(cl)) clusterEvalQ(cl,{rm(list=ls());gc()})
    
    gc()
    
    if(statistic=="Tmax") p2test=Pm
    if(statistic=="TR") p2test=Pr
    
    if(p2test>alpha|all(d_var_ij==0)){
      if(verbose){
        cat(paste("\n###########################################################################################################################\n"))
        cat(paste("Superior Set Model created\t:\n"))
        show(matrix_show)
        cat(paste("p-value\t:\n"))
        print(p2test)
        cat(paste("\n###########################################################################################################################"))
      }
      break
    }else{
      if(statistic=="Tmax") eliminate=which(v_i_M==max(v_i_M))
      if(statistic=="TR") eliminate=which(v_i_R==max(v_i_R))
      
      if(verbose) cat(paste("\nModel",model.names[eliminate],"eliminated",Sys.time()))
      Loss=as.matrix(Loss[,-eliminate])
      colnames(Loss)=model.names[-eliminate]
      # show(paste("Set rimanente:"))
      #  show(matrix_show[-which(v_i_M==max(v_i_M)),])
      #  show("p-values del test:")
      #  show(cbind(Pr,Pm))
      if(ncol(Loss)==1) {
        if(verbose){
          cat(paste("\n###########################################################################################################################\n"))
          cat(paste("Superior Set Model created\t:\n"))
          matrix_show=matrix(matrix_show[-eliminate,],nrow=1,dimnames=list(colnames(Loss),colnames(matrix_show)))
          show(matrix_show)
          cat(paste("p-value\t:\n"))
          show(p2test)
          cat(paste("\n###########################################################################################################################"))
        }
        break
      }
    }
    
    
  }    
  
  
  elapsed.time=Sys.time()-time.start
  n_elim=M_start-nrow(matrix_show)
  out=new("SSM",show=matrix_show,Info=list(model.names=model.names,elapsed.time=elapsed.time,
                                           statistic=statistic,n_elim=n_elim,mcs_pvalue=p2test,
                                           alpha=alpha,B=B,k=k),Bootstrap=list(TR=list( Stat=TR,BootDist=Tb_R),Tmax=list(Stat=TM,BootDist=Tb_M)))
  return(out)
  
  
}
