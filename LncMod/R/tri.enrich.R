tri.enrich <-
function(tri,GOterms,background,inter.thr=2,GOterms.mark=NULL,correction="BH"){
  #filter the GOterms genes to make sure that all the Goterms genes are in the background
    allgoterm<-names(GOterms)
tmpindex<-numeric()
    for(ff in 1:length(allgoterm)){
      tmp<-intersect(GOterms[[allgoterm[ff]]],background)
      if(length(tmp)!=0){
        assign(allgoterm[ff],tmp)
tmpindex<-c(tmpindex,ff)
      }
    }
allgoterm<-allgoterm[tmpindex]
  
  #an object named "goterm_freq" was constructed to contain the gene number in goterms
  goterm_freq<-as.data.frame(matrix(nrow=length(allgoterm),ncol=2))
  colnames(goterm_freq)<-c("goterm","gotarnum")
  goterm_freq[,1]<-allgoterm
  for(kk in 1:length(allgoterm)){goterm_freq[kk,2]<-length(get(allgoterm[kk]))}
  
  #an object named "lnc_freq" was constructed to contain the target gene number of an lncRNA(or more lncRNAs)
  #predict the goterms enriched by the target genes of an lncRNA(or more lncRNAs)
  lncs<-sort(as.character(unique(tri[,1])))
  lnc_freq<-as.data.frame(matrix(ncol=2,nrow=length(lncs)))
  lnc_freq[,1]<-lncs
  colnames(lnc_freq)<-c("lnc","lnctarnum")
  inter<-list()
  for(j in 1:length(lncs)){
    tmplnctar<-unique(tri[tri[,1]==lncs[j],3])
    tmplnctar_num<-length(tmplnctar)
    tmpinter<-numeric()
    for(k in allgoterm){
      tmpinter<-append(tmpinter,length(intersect(tmplnctar,get(k))))
    }#for k(GOterm name)
    inter[[j]]<-tmpinter
    lnc_freq[j,2]<-length(tmplnctar)
  }#for j(the index of lncRNA)
  interd<-unlist(inter)
  lncd<-sort(rep(lncs,length(allgoterm)))
  gotermd<-rep(allgoterm,length(lncs))
  result<-as.data.frame(cbind(lncd,gotermd))
  result<-cbind(result,interd)
  colnames(result)<-c("lnc","goterm","internum")
  result<-result[result[,"internum"]>=inter.thr,]
  result<-merge(result,goterm_freq)
  result<-merge(result,lnc_freq)
  
  #p-value adjusted by the method assigned
  backnum<-length(background)
  pvalue<-apply(result[,c("lnctarnum","gotarnum","internum")],1,function(x){return(1-phyper(x[3],x[1],backnum-x[1],x[2]))})
  result<-result[,c("lnc","goterm","lnctarnum","gotarnum","internum")]
  result<-cbind(result,pvalue)
  fdr<-p.adjust(result[,"pvalue"],method=correction)
  result<-cbind(result,fdr)
  colnames(result)<-c("modulator","GOterm","mtarnum","GOtarnum","internum","P_value","fdr")
  if(length(GOterms.mark)!=0){
    colnames(GOterms.mark)<-c("GOtermterm","mark")
    result<-merge(result,GOterms.mark)
result<-result[,c("modulator","GOterm","mtarnum","GOtarnum","internum","P_value","fdr","mark")]
    return(result)
}else{
    return(result)
  }
}
