# Julia Bischof
# 12-08-2015

#library(stringdist)

sequences.distance<-function(sequences=NULL, groups=NULL, method=c("levenshtein","cosine","q-gram","jaccard","ja-wi","dam-le","hamming","osa","lcs"),divLength=FALSE){
  if(length(sequences)==0){
    stop("--> Sequences are missing")
  }
  if(length(method)!=1){
    stop("--> Please choose one method")
  }
  
  if(length(grep(" |,|;|_|-|/",sequences))>0 && (length(groups)==0 || length(unique(groups))==1)){
    sequences<-unlist(apply(data.frame(sequences),1,function(x){strsplit(x,split=" |,|;|_|-|/")[[1]]}))
    sequences<-sequences[which(sequences!="")]
  }else if(length(grep(" |,|;|_|-|/",sequences))>0 && length(unique(groups))>1){
    temp<-data.frame(groups,sequences)
    groups<-unlist(apply(temp,1,function(x){rep(x[1],(length(gregexpr(",|;|_|-|/",x[2])[[1]])+1))}))
    sequences<-unlist(apply(data.frame(sequences),1,function(x){strsplit(x,split=" |,|;|_|-|/")[[1]]}))
    sequences<-sequences[which(sequences!="")]
    rm(temp)  }
  
  if(method=="levenshtein"){
    dist<-"lv" # Levenshtein
  }else if(method=="cosine"){
    dist<-"cosine" # cosine
  }else if(method=="q-gram"){
    dist<-"qgram" # q-gram distance (fuzzy)
  }else if(method=="jaccard"){
    dist<-"jaccard" # Jaccard distance between q-gram profiles
  }else if(method=="ja-wi"){
    dist<-"jw" # Jaro-Winkler distance (fuzzy)
  }else if(method=="osa"){
    dist<-"osa" # Optimal string alignment
  }else if(method=="dam-le"){
    dist<-"dl" # Full Damerau-Levenshtein distance.
  }else if(method=="hamming" && divLength==T){
    dist<-"hamming" # Hamming distance (a and b must have same nr of characters).
  }else if(method=="hamming" && divLength==F){
    stop("--> Sequences must have same length")
  }else if(method=="lcs"){
    dist<-"lcs" # Longest common substring distance.
  }else{
    stop("--> Method is missing or unknown")
  }
  
  out<-list()
  if(divLength==T){
    seqlength<-unlist(apply(data.frame(sequences),1,function(x){nchar(x)}))
    uniseqlength<-sort(unique(seqlength))
    for(i in 1:length(uniseqlength)){
      seq.utf8<-lapply(sequences[which(seqlength==uniseqlength[i])], utf8ToInt) # translate strings to a list of integer sequences      
      tab<-as.matrix(seq_distmatrix(seq.utf8, method = dist))
      colnames(tab)<-if(length(groups)>0){paste(groups[which(seqlength==uniseqlength[i])],sequences[which(seqlength==uniseqlength[i])],sep=": ")}else{sequences[which(seqlength==uniseqlength[i])]}
      rownames(tab)<-if(length(groups)>0){paste(groups[which(seqlength==uniseqlength[i])],sequences[which(seqlength==uniseqlength[i])],sep=": ")}else{sequences[which(seqlength==uniseqlength[i])]}
      out<-c(out,list(tab))
    }
    names(out)<-paste("Sequence length = ",uniseqlength,sep="")
  }else{
    seq.utf8<-lapply(sequences, utf8ToInt) # translate strings to a list of integer sequences
    tab<-as.matrix(seq_distmatrix(seq.utf8, method = dist))
    colnames(tab)<-if(length(groups)>0){paste(groups,sequences,sep=": ")}else{sequences}
    rownames(tab)<-if(length(groups)>0){paste(groups,sequences,sep=": ")}else{sequences}
    out<-c(out,list(tab))
  }
  return(out)
}
