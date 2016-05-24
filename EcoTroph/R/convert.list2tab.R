#--------------------------------------------------------------------------------------------
# convert.list2tab
# Convert the list object returned by the create.ETdiagnosis function into data.frame objects
#--------------------------------------------------------------------------------------------

convert.list2tab=function(diagn.list){

n.fleet=length(diagn.list[[1]]$Catches) 
name.fleet=names(diagn.list[[1]]$Catches)
TL=names(diagn.list[[1]]$BIOM_MF)

fleet.of.interest=diagn.list[['fleet.of.interest']]
if(!is.null(fleet.of.interest)){diagn.list=diagn.list[-length(diagn.list)]}

nam.liste=names(diagn.list)
N=length(nam.liste)

if(!n.fleet==1){
  if(!is.null(fleet.of.interest)){
    mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)
    for(i in 1:n.fleet){if(i==1){m.=mf[,1]}else{m.=c(m.,mf[,1])}}
    mf=as.numeric(unique(sort(m.)))
    for(i in 1:n.fleet){
      if(name.fleet[i]%in%fleet.of.interest){nm=mf}else{nm=rep(1,length(mf))}
      if(i==1){nm.interest=nm}else{nm.interest=paste(nm.interest,'_',nm,sep='')}
    }
  }
}

#-------------------------------------------------------------------------------
#---------------------------- 1. ET_Main_mf ------------------------------------
#-------------------------------------------------------------------------------
name.var=c('B','B_acc','P','P_acc','Kin','Kin_acc',
           'Fish_mort','Fish_mort_acc',#'N_loss','N_loss_acc',
           'F_loss','F_loss_acc','Y')

var=c('BIOM_MF','BIOM_MF_acc','Prod_MF','Prod_MF_acc','Kin_MF','Kin_MF_acc',
      'Fish_mort','Fish_mort_acc',#'N_loss','N_loss_acc',
      'F_loss','F_loss_acc','Catches.tot')

if(n.fleet==1){

  ET_Main_mE=list()
  for(v in 1:length(var)){
    #v=3
    #print(v)
    for(i in 1:N){
      #i=1
      xlist=diagn.list[[nam.liste[i]]]
      Mf=as.numeric(xlist$mf) 
      et=c(Mf,xlist[[var[v]]])
      if(i==1){Et=et}else{Et=cbind(Et,et)}
    }
    Et=cbind(c('mE',TL),Et)
    Et=data.frame(Et)
    #row.names(Et)=NULL;
    #colnames(Et)=NULL
    ET_Main_mE[[name.var[v]]]=Et
  }
  
}else{# n.fleet>1
  
  if(is.null(fleet.of.interest)){
    
    mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)
    for(i in 1:n.fleet){if(i==1){m.=mf[,1]}else{m.=c(m.,mf[,1])}}
    mf=as.numeric(unique(sort(m.)))
    for(n in 1:n.fleet){
      if(n==1){nm=mf}else{nm=paste(nm,'_',mf,sep='')}
    }
    nm=as.character(nm)
    
    if(length(names(diagn.list)[!names(diagn.list)%in%nm])>length(nm)){ #same.mE=F

      ET_Main_mE=list()
      
      for(v in 1:length(var)){
        #v=3
        #print(v)
        for(i in 1:N){
          #i=1
          xlist=diagn.list[[nam.liste[i]]]
          Mf=as.numeric(xlist$mf) 
          et=c(Mf,xlist[[var[v]]])
          if(i==1){Et=et}else{Et=cbind(Et,et)}
        }
        Et=cbind(c(paste('mE.fleet',1:n.fleet,sep=''),TL),Et)
        #Et=data.frame(Et)
        #row.names(Et)=NULL;colnames(Et)=NULL
        ET_Main_mE[[name.var[v]]]=Et
      }
      
   #   for(f in 1:n.fleet){
    #    for(i in 1:N){
    #      #i=1
    #      xlist=diagn.list[[nam.liste[i]]]
    #      Mf=xlist$mf  
    #      et=c(Mf,xlist$Catches[[name.fleet[f]]])
    #      if(i==1){Et=et}else{Et=cbind(Et,et)}
    #    }
    #    Et=cbind(c(paste('mE.fleet',1:n.fleet,sep=''),TL),Et)
        #Et=data.frame(Et)
        #row.names(Et)=NULL;colnames(Et)=NULL
    #    ET_Main_mE[[paste('Y_fleet',f,sep='')]]=Et
    #  }
      
    }else{ #same.mE=T

      ET_Main_mE=list()
      
      for(v in 1:length(var)){
        #v=3
        #print(v)
        for(i in 1:N){
          #i=1
          xlist=diagn.list[[nam.liste[i]]]
          Mf=as.numeric(xlist$mf)[1] 
          et=c(Mf,xlist[[var[v]]])
          if(i==1){Et=et}else{Et=cbind(Et,et)}
        }
        Et=cbind(c('mE',TL),Et)
        #Et=data.frame(Et)
        #row.names(Et)=NULL;colnames(Et)=NULL
        ET_Main_mE[[name.var[v]]]=Et
      }
      
      #for(f in 1:n.fleet){
      #  for(i in 1:N){
          #i=1
      #    xlist=diagn.list[[nam.liste[i]]]
      #    Mf=xlist$mf[1]  
      #    et=c(Mf,xlist$Catches[[name.fleet[f]]])
      #    if(i==1){Et=et}else{Et=cbind(Et,et)}
      #  }
     #   Et=cbind(c('mE',TL),Et)
        #Et=data.frame(Et)
        #row.names(Et)=NULL;colnames(Et)=NULL
    #    ET_Main_mE[[paste('Y_fleet',f,sep='')]]=Et
     # }   
      
    }
  }else{ # !is.null(fleet.of.interest)
    
    mE=matrix(data=unlist(strsplit(nm.interest,'_')),ncol=n.fleet,nrow=length(nm.interest),byrow=T)
    mE=sort(unique(as.numeric(mE)))
    ET_Main_mE=list()
    
    for(v in 1:length(var)){
     # print(v)
      for(i in 1:length(mE)){
        #i=1
        for(k in 1:n.fleet){if(name.fleet[k]%in%fleet.of.interest){nm.=as.character(mE[i])}else{nm.='1'}
                            if(k==1){Nm.=nm.}else{Nm.=paste(Nm.,'_',nm.,sep='')}}
        xlist=diagn.list[[Nm.]]
        et=c(mE[i],xlist[[var[v]]])
        if(i==1){Et=et}else{Et=cbind(Et,et)}
      }
      Et=cbind(c('mE.fleet.of.interest',TL),Et)
      #Et=data.frame(Et)
      #row.names(Et)=NULL;colnames(Et)=NULL
      ET_Main_mE[[name.var[v]]]=Et
    }
    
  #  for(f in 1:n.fleet){
  #    for(i in 1:length(mE)){
        #i=1
  #      xlist=diagn.list[[paste(mE[i],'_1',sep='')]]
  #      et=c(mE[i],xlist$Catches[[name.fleet[f]]])
  #      if(i==1){Et=et}else{Et=cbind(Et,et)}
  #    }
  #    Et=cbind(c('mE.fleet.of.interest',TL),Et)
  #    #Et=data.frame(Et)
      #row.names(Et)=NULL;
  #    colnames(Et)=NULL
  #    ET_Main_mE[[paste('Y_fleet',f,sep='')]]=Et
  #  }
    
  } 
}   

#-------------------------------------------------------------------------------
#---------------------------- 2. ET_Main_diagnose ------------------------------
#-------------------------------------------------------------------------------

name.diagnose=names(diagn.list[[1]]$ET_Main_diagnose)

if(n.fleet==1){
  
  var2=c('TOT_B','TOT_B_acc','TOT_P','TOT_P_acc','Y',
         'TL_TOT_B','TL_TOT_B_acc','TL_Y')
  
  name.col2=c('mE',var2)
  
  for(i in 1:N){
    #i=1
    xlist=diagn.list[[nam.liste[i]]]
    MF=xlist$mf
    
    for(v in var2){
      e.=xlist$ET_Main_diagnose[[v]];if(v==var2[1]){et=e.}else{et=c(et,e.)} 
    }
    
    et=c(MF,et)
    if(i==1){Et=et}else{Et=rbind(Et,et)}
  }
  ET_Main_diagnose=Et
  colnames(ET_Main_diagnose)=name.col2
  row.names(ET_Main_diagnose)=NULL
  ET_Main_diagnose=cbind(colnames(ET_Main_diagnose),t(ET_Main_diagnose))
  #row.names(ET_Main_diagnose)=NULL;colnames(ET_Main_diagnose)=NULL
  
}else{# n.fleet>1
  if(is.null(fleet.of.interest)){
    
    mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)
    for(i in 1:n.fleet){if(i==1){m.=mf[,1]}else{m.=c(m.,mf[,1])}}
    mf=as.numeric(unique(sort(m.)))
    for(n in 1:n.fleet){
      if(n==1){nm=mf}else{nm=paste(nm,'_',mf,sep='')}
    }
    nm=as.character(nm)
    
    if(length(names(diagn.list)[!names(diagn.list)%in%nm])>length(nm)){#same.mE=F
      var2=c('TOT_B','TOT_B_acc','TOT_P','TOT_P_acc','Y',
             name.diagnose[substring(name.diagnose,1,2)=='Y_'],
             'TL_TOT_B','TL_TOT_B_acc','TL_Y',
             name.diagnose[substring(name.diagnose,1,5)=='TL_Y_'])
      
      name.col2=c(paste('mE.fleet',1:n.fleet,sep=''),var2)
      
      for(i in 1:N){
        #i=1
        xlist=diagn.list[[nam.liste[i]]]
        MF=xlist$mf
        
        for(v in var2){
          e.=xlist$ET_Main_diagnose[[v]];if(v==var2[1]){et=e.}else{et=c(et,e.)} 
        }
        
        for(f in 1:n.fleet){m.=MF[f];if(f==1){Mf=m.}else{Mf=c(Mf,m.)}}
        
        et=c(Mf,et)
        if(i==1){Et=et}else{Et=rbind(Et,et)}
      }
      ET_Main_diagnose=Et
      colnames(ET_Main_diagnose)=name.col2
      row.names(ET_Main_diagnose)=NULL
      ET_Main_diagnose=cbind(colnames(ET_Main_diagnose),t(ET_Main_diagnose))
      #row.names(ET_Main_diagnose)=NULL;colnames(ET_Main_diagnose)=NULL      
      
    }else{# same.mE=T
      var2=c('TOT_B','TOT_B_acc','TOT_P','TOT_P_acc','Y',
             name.diagnose[substring(name.diagnose,1,2)=='Y_'],
             'TL_TOT_B','TL_TOT_B_acc','TL_Y',
             name.diagnose[substring(name.diagnose,1,5)=='TL_Y_'])
      
      name.col2=c('mE',var2)
      
      for(i in 1:length(nm)){
        #i=1
        xlist=diagn.list[[nm[i]]]
        MF=xlist$mf
        
        for(v in var2){
          e.=xlist$ET_Main_diagnose[[v]];if(v==var2[1]){et=e.}else{et=c(et,e.)} 
        }
        
        #for(f in 1:n.fleet){m.=MF[f];if(f==1){Mf=m.}else{Mf=c(Mf,m.)}}
        
        et=c(mf[i],et)
        if(i==1){Et=et}else{Et=rbind(Et,et)}
      }
      ET_Main_diagnose=Et
      colnames(ET_Main_diagnose)=name.col2
      row.names(ET_Main_diagnose)=NULL
      ET_Main_diagnose=cbind(colnames(ET_Main_diagnose),t(ET_Main_diagnose))
      #row.names(ET_Main_diagnose)=NULL;colnames(ET_Main_diagnose)=NULL
      
    } 
  }else{#!is.null(fleet.of.interest)
    var2=c('TOT_B','TOT_B_acc','TOT_P','TOT_P_acc','Y',
           name.diagnose[substring(name.diagnose,1,2)=='Y_'],
           'TL_TOT_B','TL_TOT_B_acc','TL_Y',
           name.diagnose[substring(name.diagnose,1,5)=='TL_Y_']
    )
    
    mE=matrix(data=unlist(strsplit(nm.interest,'_')),ncol=n.fleet,nrow=length(nm.interest),byrow=T)
    mE=unique(sort(as.numeric(mE)))
    ET_Main_diagnose=data.frame(matrix(data=NA,ncol=length(mE),nrow=length(var2)))
    
    for(i in 1:length(mE)){
      #i=1
      for(k in 1:n.fleet){if(name.fleet[k]%in%fleet.of.interest){nm.=as.character(mE[i])}else{nm.='1'}
                          if(k==1){Nm.=nm.}else{Nm.=paste(Nm.,'_',nm.,sep='')}}
      xlist=diagn.list[[Nm.]]$ET_Main_diagnose
      for(v in var2){
        vv=xlist[[v]];if(v==var2[1]){vv.=vv}else{vv.=c(vv.,vv)}
      }
      ET_Main_diagnose[,i]=vv. 
    }
    
    ET_Main_diagnose=cbind(c('mE.fleet.of.interest',var2),rbind(mE,ET_Main_diagnose))
    #colnames(ET_Main_diagnose)=NULL;row.names(ET_Main_diagnose)=NULL
  }
}    
tab=ET_Main_mE
tab[['ET_Main_diagnose']]=data.frame(ET_Main_diagnose)
return(tab)
}