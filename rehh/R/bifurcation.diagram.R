bifurcation.diagram<-function(haplohh,mrk_foc,all_foc=1,nmrk_l=10,nmrk_r=10,limhapcount=10,refsize=0.1,linecol="blue",main_leg=NA,xlab_leg="Position"){

if(!(is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")}
if(nmrk_l<0 | nmrk_r<0){stop("nmrk_l and nmrk_r must be positive or null")}
if(limhapcount<1){stop("limhapcount must be >1")}

if(is.na(main_leg)){
 if(all_foc==1){main_leg=paste(haplohh@snp.name[mrk_foc], " (Ancestral Allele)",sep="")}
 if(all_foc==2){main_leg=paste(haplohh@snp.name[mrk_foc], " (Derived Allele)",sep="")}
}
#Checks
if(mrk_foc-nmrk_l<0){
stop("Too much markers on the left")
}
if(mrk_foc+nmrk_r>haplohh@nsnp){
stop("Too much markers on the right")
}

#Recup haplos on the Right
if(nmrk_r>0){
  haplo_r=haplohh@haplo[haplohh@haplo[,mrk_foc]==all_foc,][,(mrk_foc+1):(mrk_foc+nmrk_r)]
  haplo_r=haplo_r[rowSums(haplo_r==0)==0,]
  if(nrow(haplo_r)<limhapcount){stop("Number of available haplotypes on the right lower than limhapcount")}
haplo_der_r=matrix(0,0,4) ; det_haplo_r=matrix(0,1,4) #; haplo_anc=matrix(0,0,4)
colnames(haplo_der_r)=c("CODE_HAP","NhaploDer","HAP_DER1","HAP_DER2") #; colnames(haplo_anc)=c("CODE_HAP","HAP_ANC","n","mrk")
colnames(det_haplo_r)=c("CODE_HAP","HAPLO","mrk","n")
det_haplo_r[1,]=as.matrix(c(1,"",0,nrow(haplo_r)))

for(mrk in 1:nmrk_r){
  nhaplo_tot=nrow(det_haplo_r)
  if(mrk==1){
    tmp=haplo_r[,1]
     }else{
    tmp=apply(haplo_r[,1:mrk],1,paste,collapse="")
     }
  list_hap_0=matrix(det_haplo_r[as.numeric(det_haplo_r[,3])==mrk-1,],ncol=4) #haplo au mrk n-1
  for(i in 1:nrow(list_hap_0)){
   hap_der=list()
   hap_der[[1]]=paste(list_hap_0[i,2],"1",sep="") ; hap_der[[2]]=paste(list_hap_0[i,2],"2",sep="")
   nhap_der=c(sum(tmp==hap_der[[1]]),sum(tmp==hap_der[[2]]))
   tmp_nhap_der=sum(nhap_der>0) ; tmp_hap_der=rep(0,2) ; tmp_cnt=0
   for(all in 1:2){
     if(nhap_der[all]>0){
        nhaplo_tot=nhaplo_tot+1
        det_haplo_r=rbind(det_haplo_r,c(nhaplo_tot,hap_der[[all]],mrk,nhap_der[all])) #; haplo_anc=rbind(haplo_anc,c(nhaplo_tot,list_hap_0[i,1],nhap_der[all],mrk))
        tmp_cnt=tmp_cnt+1 ; tmp_hap_der[tmp_cnt]=nhaplo_tot
   }
  }
        haplo_der_r=rbind(haplo_der_r,c(list_hap_0[i,1],tmp_nhap_der,tmp_hap_der))
}
}

rownames(haplo_der_r)=haplo_der_r[,1] # ; rownames(haplo_anc)=haplo_anc[,1] 
haplo_der_r=haplo_der_r[,-1] #; haplo_anc=haplo_anc[,-1]  
#calcul coordonnees des haplo
coord_r=matrix(0,nrow(det_haplo_r),2) ; rownames(coord_r)=det_haplo_r[,1] ; colnames(coord_r)=c("X","Y")
for(i in nmrk_r:0){
  tmp_haplo=det_haplo_r[as.numeric(det_haplo_r[,3])==i,1]
  coord_r[tmp_haplo,1]=haplohh@position[i+mrk_foc]
  if(i==nmrk_r){
    coord_r[tmp_haplo,2]=1:length(tmp_haplo)/2
     }else{
    for(hap in tmp_haplo){
     if(as.numeric(haplo_der_r[hap,1])==1){coord_r[hap,2]=coord_r[haplo_der_r[hap,2],2]}
     if(as.numeric(haplo_der_r[hap,1])==2){coord_r[hap,2]=mean(coord_r[haplo_der_r[hap,2:3],2])}
     coord_r[hap,1]=haplohh@position[i+mrk_foc]
  }
}
}
}

#Recup haplos on the Left
if(nmrk_l>0){
  haplo_l=haplohh@haplo[haplohh@haplo[,mrk_foc]==all_foc,][,(mrk_foc-1):(mrk_foc-nmrk_l)]
  haplo_l=haplo_l[rowSums(haplo_l==0)==0,]
  if(nrow(haplo_l)<limhapcount){stop("Number of available haplotypes on the left lower than limhapcount")}
haplo_der_l=matrix(0,0,4) ; det_haplo_l=matrix(0,1,4) #; haplo_anc=matrix(0,0,4)
colnames(haplo_der_l)=c("CODE_HAP","NhaploDer","HAP_DER1","HAP_DER2") #; colnames(haplo_anc)=c("CODE_HAP","HAP_ANC","n","mrk")
colnames(det_haplo_l)=c("CODE_HAP","HAPLO","mrk","n")
det_haplo_l[1,]=as.matrix(c(1,"",0,nrow(haplo_l)))

for(mrk in 1:nmrk_l){
  nhaplo_tot=nrow(det_haplo_l)
  if(mrk==1){
    tmp=haplo_l[,1]
     }else{
    tmp=apply(haplo_l[,1:mrk],1,paste,collapse="")
     }
  list_hap_0=matrix(det_haplo_l[as.numeric(det_haplo_l[,3])==mrk-1,],ncol=4) #haplo au mrk n-1
  for(i in 1:nrow(list_hap_0)){
   hap_der=list()
   hap_der[[1]]=paste(list_hap_0[i,2],"1",sep="") ; hap_der[[2]]=paste(list_hap_0[i,2],"2",sep="")
   nhap_der=c(sum(tmp==hap_der[[1]]),sum(tmp==hap_der[[2]]))
   tmp_nhap_der=sum(nhap_der>0) ; tmp_hap_der=rep(0,2) ; tmp_cnt=0
   for(all in 1:2){
     if(nhap_der[all]>0){
        nhaplo_tot=nhaplo_tot+1
        det_haplo_l=rbind(det_haplo_l,c(nhaplo_tot,hap_der[[all]],mrk,nhap_der[all])) #; haplo_anc=rbind(haplo_anc,c(nhaplo_tot,list_hap_0[i,1],nhap_der[all],mrk))
        tmp_cnt=tmp_cnt+1 ; tmp_hap_der[tmp_cnt]=nhaplo_tot
   }
  }
        haplo_der_l=rbind(haplo_der_l,c(list_hap_0[i,1],tmp_nhap_der,tmp_hap_der))
}
}

rownames(haplo_der_l)=haplo_der_l[,1] # ; rownames(haplo_anc)=haplo_anc[,1] 
haplo_der_l=haplo_der_l[,-1] #; haplo_anc=haplo_anc[,-1]  
#calcul coordonnees des haplo
coord_l=matrix(0,nrow(det_haplo_l),2) ; rownames(coord_l)=det_haplo_l[,1] ; colnames(coord_l)=c("X","Y")
for(i in nmrk_l:0){
  tmp_haplo=det_haplo_l[as.numeric(det_haplo_l[,3])==i,1]
  coord_l[tmp_haplo,1]=haplohh@position[mrk_foc-i]
  if(i==nmrk_l){
    coord_l[tmp_haplo,2]=1:length(tmp_haplo)/2
     }else{
    for(hap in tmp_haplo){
     if(as.numeric(haplo_der_l[hap,1])==1){coord_l[hap,2]=coord_l[haplo_der_l[hap,2],2]}
     if(as.numeric(haplo_der_l[hap,1])==2){coord_l[hap,2]=mean(coord_l[haplo_der_l[hap,2:3],2])}
     coord_l[hap,1]=haplohh@position[mrk_foc-i]
  }
}
}
}

#PLOT
if(nmrk_l>0 & nmrk_r>0){
  coord_l[,2]=coord_l[,2] + coord_r[1,2] - coord_l[1,2]
  dum_coord=rbind(coord_l,coord_r)
}
if(nmrk_l==0){dum_coord=coord_r}
if(nmrk_r==0){dum_coord=coord_l}

  plot(dum_coord,pch="",yaxt="n",bty="n",xlab=xlab_leg, main=main_leg,ylab="")
  abline(v=haplohh@position[mrk_foc],lty=2,col=linecol)

if(nmrk_r>0){
  lwd_adjust=as.numeric(det_haplo_r[,4])*refsize/(as.numeric(det_haplo_r[1,1]))
  names(lwd_adjust)=rownames(coord_r)
  for(i in (nmrk_r-1):0){
   tmp_haplo=det_haplo_r[as.numeric(det_haplo_r[,3])==i,1]
    for(hap in tmp_haplo){
     x0=coord_r[hap,1] ; y0=coord_r[hap,2] 
      for(j in 1:as.numeric(haplo_der_r[hap,1])){
        tmp_lwd=lwd_adjust[haplo_der_r[hap,1+j]]
        x1=coord_r[haplo_der_r[hap,1+j],1] ; y1=coord_r[haplo_der_r[hap,1+j],2]
        lines(c(x0,x1),c(y0,y1),lwd=tmp_lwd,lty=1,col=linecol)
     }
  }
}
}

if(nmrk_l>0){
 lwd_adjust=as.numeric(det_haplo_l[,4])*refsize/(as.numeric(det_haplo_l[1,1]))
 names(lwd_adjust)=rownames(coord_l)
 for(i in (nmrk_l-1):0){
  tmp_haplo=det_haplo_l[as.numeric(det_haplo_l[,3])==i,1]
    for(hap in tmp_haplo){
     x0=coord_l[hap,1] ; y0=coord_l[hap,2] 
      for(j in 1:as.numeric(haplo_der_l[hap,1])){
        tmp_lwd=lwd_adjust[haplo_der_l[hap,1+j]]
        x1=coord_l[haplo_der_l[hap,1+j],1] ; y1=coord_l[haplo_der_l[hap,1+j],2]
        lines(c(x0,x1),c(y0,y1),lwd=tmp_lwd,lty=1,col=linecol)
     }
  }
}
}

}
