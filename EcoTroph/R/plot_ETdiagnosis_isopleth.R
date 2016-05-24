plot_ETdiagnosis_isopleth=function(x,fleet1=NULL,fleet2=NULL,var=NULL,n.level=NULL,relative=NULL,name.fleet1=NULL,name.fleet2=NULL,color=NULL,ask=interactive()){
 diagn.list<-x
par(ask=ask)
 # diagn.list=Liste#; fleet1='catch.1';fleet2='catch.2'
 # var=NULL;levels=NULL;name.fleet1=NULL;name.fleet2=NULL;color=NULL;relative=NULL;n.level=NULL
  fleet.of.interest=diagn.list[['fleet.of.interest']]
  if(!is.null(fleet.of.interest)){diagn.list=diagn.list[-length(diagn.list)]}
  #diagn.list=diagn.list[-length(diagn.list)]
  if(is.null(var)){var=c('TOT_B','TOT_B_acc','TOT_P','TOT_P_acc','Y','Y_fleet1','Y_fleet2','TL_TOT_B',
                         'TL_TOT_B_acc','TL_Y','TL_Y_fleet1','TL_Y_fleet2')}
  if(is.null(relative)){relative=F}
  if(relative){var=paste('R_',var[!var%in%c('TL_TOT_B','TL_TOT_B_acc','TOT_P','TOT_P_acc','TL_Y','TL_Y_fleet1','TL_Y_fleet2')],sep='')}
  if(is.null(fleet.of.interest)){
      if(is.null(name.fleet1)){x.lab='fleet 1'}else{x.lab=name.fleet1}
      if(is.null(name.fleet2)){y.lab='fleet 2'}else{y.lab=name.fleet2}
  }else{x.lab='fleet(s) of interest';y.lab='other fleet(s)' }
  if(is.null(n.level)){n.level=7}
  if(is.null(color)) color=rainbow((n.level-1))#rev(gray(seq(0,1,1/(levels-1))))
  # find each combination of diagnlist to conserve
  
  # fleets
  fleets=names(diagn.list[[1]][['Catches']]); n.fleet=length(fleets); names(fleets)=1:n.fleet
  if(is.null(fleet.of.interest)){
    if(is.null(fleet2)){fleet2=fleets[!fleets %in% fleet1]}
    neutral.fleet=as.character(fleets[!fleets %in% c(fleet1,fleet2)])
  }else{
    fleet1=fleet.of.interest
    fleet2=fleets[!fleets%in%fleet1]
    neutral.fleet=NULL
  }
  FL=list(fleet1=fleet1,fleet2=fleet2)
  
  # Mul_eff
  mf=matrix(data=unlist(strsplit(names(diagn.list),'_')),ncol=n.fleet,nrow=length(diagn.list),byrow=T)
  for(i in 1:n.fleet){if(i==1){m.=mf[,1]}else{m.=c(m.,mf[,1])}}
  mf=sort(as.numeric(unique(m.)));names(mf)=1:length(mf)

  for(i in 1:length(mf)){
    for(j in 1:length(mf)){
      for(k in 1:length(fleets)){
        if(fleets[k]%in%fleet1){Mf=mf[i]}
        if(fleets[k]%in%fleet2){Mf=mf[j]}
        if(fleets[k]%in%neutral.fleet){Mf=1}#;Mf2=NULL
        if(k==1){n.=Mf}else{n.=paste(n.,'_',Mf,sep='')}
      }
      if(j==1){n..=n.;nl=paste(mf[i],'_',mf[j],sep='')}else{n..=c(n..,n.);nl=c(nl,nl=paste(mf[i],'_',mf[j],sep=''))}
    }
    if(i==1){nam=n..;nl.=nl}else{nam=c(nam,n..);nl.=c(nl.,nl)}
  }
  # renommer la liste pour trouver les bonnes combinaisons
  N=cbind(nam,nl.);row.names(N)=NULL;N=data.frame(N)
  ll=diagn.list[names(diagn.list)%in%nam]
  #nll=names(ll)
  #for(n in 1:length(nll)){names(ll[n])=as.character(N[N$nam==nll[n],'nl.'])}
  names(ll)<-N[match(N$nam,names(ll)),]$nl.
  
  # create a matrix for each variable and then plot an ispoleth
  x=mf;y=mf
  nmf1=as.numeric(names(mf[mf==1]))
#  nv=0
  var.mat=list()
  for(v in var){
   # print(v)
  #  nv=nv+1
  #  v='Y'
    #v=var[11];v
    M=matrix(data=NA,ncol=length(mf),nrow=length(mf))
    
    if(!v%in%c('Y_fleet1','Y_fleet2','R_Y_fleet1','R_Y_fleet2','TL_Y_fleet1','TL_Y_fleet2')){
      for(i in 1:length(mf)){# fleet1
        for(j in 1:length(mf)){# fleet2
          M[i,j]=ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[v]]
        }
      }
    }else{
      if(v%in%c('Y_fleet1','Y_fleet2')){
        for(i in 1:length(mf)){
          for(j in 1:length(mf)){
            for(fl in 1:length(FL[[substring(v,3,8)]])){
              m=ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[paste('Y_',strsplit(FL[[substring(v,3,8)]][fl],'catch.')[[1]][2],sep='')]]
              if(fl==1){m.=m}else{m.=m.+m}}
            M[i,j]=m.
          }
        }
      }
      if(v%in%c('R_Y_fleet1','R_Y_fleet2')){
        for(i in 1:length(mf)){
          for(j in 1:length(mf)){
            for(fl in 1:length(FL[[substring(v,5,10)]])){
              m=ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[paste('Y_',strsplit(FL[[substring(v,5,10)]][fl],'catch.')[[1]][2],sep='')]]
              if(fl==1){m.=m}else{m.=m.+m}}
            M[i,j]=m.
          }
        }
        M=M/M[nmf1,nmf1]
      }
      if(v%in%c('TL_Y_fleet1','TL_Y_fleet2')){
        for(i in 1:length(mf)){
          for(j in 1:length(mf)){
            for(fl in 1:length(FL[[substring(v,6,11)]])){
              n=ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[paste('Y_',strsplit(FL[[substring(v,6,11)]][fl],'catch.')[[1]][2],sep='')]]*ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[paste('TL_Y_',strsplit(FL[[substring(v,6,11)]][fl],'catch.')[[1]][2],sep='')]]
              d=ll[[paste(mf[i],'_',mf[j],sep='')]][['ET_Main_diagnose']][[paste('Y_',strsplit(FL[[substring(v,6,11)]][fl],'catch.')[[1]][2],sep='')]]
              if(fl==1){n.=n;d.=d}else{n.=n.+n;d.=d.+d}
              if(fl==length(FL[[substring(v,6,11)]])){m.=n./d.}}
            M[i,j]=m.
          }
        }
      }  
    }
    
    # save matrix
    row.names(M)=mf;colnames(M)=mf
    var.mat[[v]]=M
  
    par(mar=c(5, 4.5, 4, 1))
   #x11()
    if(is.null(fleet.of.interest)){v.=v}else{if(!v%in%c('Y_fleet1','Y_fleet2','R_Y_fleet1','R_Y_fleet2','TL_Y_fleet1','TL_Y_fleet2')){v.=v}else{
      if(v=='TL_Y_fleet1'){v.='TL_Y of fleet(s) of interest'};if(v=='TL_Y_fleet2'){v.='TL_Y of other fleet(s)'}
      if(v=='Y_fleet1'){v.='Y of fleet(s) of interest'};if(v=='Y_fleet2'){v.='Y of other fleet(s)'}
      if(v=='R_Y_fleet1'){v.='R_Y of fleet(s) of interest'};if(v=='R_Y_fleet2'){v.='R_Y of other fleet(s)'}
    }}
   filled.contour(x, y, M, nlevels = n.level, col = color,
   #filled.contour(x, y, t(M), nlevels = n.level, col = color,
   asp = 1, xlab = paste("mE of",x.lab), ylab = paste("mE of",y.lab), axes = TRUE, main = v.)
  #savePlot(paste('D:/Documents and Settings/gatti/Bureau/iso/iso',nv,'.png',sep=''),type='png');dev.off()
  #  points(1,1,type = "p", pch=21, bg="white")
  }
  for(i in 1:length(var)){var.mat[[i]]=t(var.mat[[i]])}
  return(var.mat)
}