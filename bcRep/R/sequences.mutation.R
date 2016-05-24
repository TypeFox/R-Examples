## Julia Bischof
## 10-09-2015

sequences.mutation<-function(mutationtab=NULL,summarytab=NULL,sequence=c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"),
                             functionality=FALSE,junctionFr=FALSE, rsRatio=FALSE,...){
  if(length(mutationtab)==0){
    stop('--> 7_V-REGION-mutation-and-AA-change-table file is missing')
  }
  if(length(sequence)!=1 || !(sequence %in% c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"))){
    stop('--> Sequence information (f.e. V, FR1, CDR2,...) is missing')
  }
  out.list<-list()
  if(sequence=='V'){
    column<-paste(sequence,'_REGION',sep="")  
  }else{
    column<-paste(sequence,'_IMGT',sep="")
  }
  mut.tab<-cbind(apply(data.frame(mutationtab[,column]),1,function(x){length(which(gregexpr("[|]",x)[[1]]>0))}),
                 apply(data.frame(mutationtab[,column]),1,function(x){length(which(gregexpr(",",x)[[1]]>0))}))
  mut.tab<-cbind(mut.tab,mut.tab[,1]-mut.tab[,2])
  if(length(summarytab) == 0 || sequence!="V") {
    mut.tab <- data.frame(mut.tab)
    colnames(mut.tab)<-c("number_mutations","number_replacement","number_silent")
    if(rsRatio==T){
      mut.tab <- cbind(mut.tab,mut.tab[,2]/mut.tab[,3])
      mut.tab[which(mut.tab[,4] %in% c("Inf","NaN")),4]<-0
      colnames(mut.tab)<-c("number_mutations","number_replacement","number_silent","RS_ratio")
      rownames(mut.tab)<-mutationtab$Sequence_ID
    }
  }
  else if(length(summarytab) != 0 && sequence=="V"){
    mut.tab <- data.frame(summarytab$V_REGION_identity_nt, mut.tab)
    colnames(mut.tab)<-c("V_REGION_identity_nt","number_mutations","number_replacement","number_silent")
    if(rsRatio==T){
      mut.tab <- cbind(mut.tab,mut.tab[,3]/mut.tab[,4])
      mut.tab[which(mut.tab[,5] %in% c("Inf","NaN")),5]<-0
      colnames(mut.tab)<-c("V_REGION_identity_nt","number_mutations","number_replacement","number_silent","RS_ratio")
      rownames(mut.tab)<-mutationtab$Sequence_ID
    }
  }
  if(functionality==F && junctionFr == F){
    out.list<-data.frame(mut.tab,row.names=NULL)
  }else{
    out.list<-c(out.list,list(mut.tab))
    names(out.list)<-"Number_of_mutations"
  }

   if(functionality==TRUE){
    nrmut.prod.mut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0),grep("^productive",mutationtab$Functionality)))
    nrmut.unprod.mut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0),grep("^unproductive",mutationtab$Functionality)))
    nrmut.prod.nomut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0),grep("^productive",mutationtab$Functionality)))
    nrmut.unprod.nomut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0),grep("^unproductive",mutationtab$Functionality)))
    nrmut.other.mut<-length(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0))-nrmut.prod.mut-nrmut.unprod.mut
    nrmut.other.nomut<-length(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0))-nrmut.prod.nomut-nrmut.unprod.nomut
    func.tab<-t(data.frame(nrmut.prod.mut,nrmut.unprod.mut,nrmut.other.mut,nrmut.prod.nomut,nrmut.unprod.nomut,nrmut.other.nomut))
    rownames(func.tab)<-c("mutation_in_productive_sequences","mutations_in_unproductive_sequences","mutations_in_sequences_with_unknown_functionality","no_mutation_in_productive_sequences","no_mutations_in_unproductive_sequences","no_mutations_in_sequences_with_unknown_functionality")
    colnames(func.tab)<-"proportion"
    func.tab[,1]<-func.tab[,1]/nrow(mutationtab)
    out.list<-c(out.list,list(data.frame(func.tab,check.names=F)))
    names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"Functionality")
  }
  if(junctionFr==TRUE){
    nrmut.inframe.mut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0),grep("^in-frame",summarytab$JUNCTION_frame)))
    nrmut.outframe.mut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0),grep("^out-of-frame",summarytab$JUNCTION_frame)))
    nrmut.inframe.nomut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0),grep("^in-frame",summarytab$JUNCTION_frame)))
    nrmut.outframe.nomut<-length(intersect(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0),grep("^out-of-frame",summarytab$JUNCTION_frame)))
    nrmut.other.mut<-length(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]>0))-nrmut.inframe.mut-nrmut.outframe.mut
    nrmut.other.nomut<-length(which(mut.tab[,which(colnames(mut.tab)=="number_mutations")]==0))-nrmut.inframe.nomut-nrmut.outframe.nomut
    jf.tab<-t(data.frame(nrmut.inframe.mut,nrmut.outframe.mut,nrmut.other.mut,nrmut.inframe.nomut,nrmut.outframe.nomut,nrmut.other.nomut))
    rownames(jf.tab)<-c("mutation_in_in-frame_sequences","mutations_in_out-of-frame_sequences","mutations_in_sequences_with_unknown_junction_frame","no_mutation_in_in-frame_sequences","no_mutations_in_out-of-frame_sequences","no_mutations_in_sequences_with_unknown_junction_frame")
    jf.tab[,1]<-jf.tab[,1]/nrow(mutationtab)
    colnames(jf.tab)<-"proportion"
    out.list<-c(out.list,list(data.frame(jf.tab,check.names=F)))
    names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"Junction_frame")
  }else if(length(summarytab)==0 && junctionFr==TRUE){
    stop("--> 1_Summary file is missing, no Junction frame information available")
  }
  return(out.list)
}


