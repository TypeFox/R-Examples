sfsplit<-function(sform){
if(nchar(sform)==1){stop("S-form is improperly formed, requires > 1 event symbol")}
if(grepl("[[:punct:]]|[[:digit:]]",sform)){stop("S-form is improperly formed, no special characters allowed")}
 sf.split<-unlist(strsplit(sform,""))
 blue<-vector("list",length(sf.split)-1)
 for(i in 2:length(sf.split)){
   sf.red<-combn(sf.split,i)[,1]
   lnA<-length(sf.red)
   blue[[i-1]]<-paste(sf.red[-lnA],collapse="")
   names(blue)[i-1]<-sf.split[i]
 }
# blue<-sapply(rev(sf.split),function(x) strsplit(sform,x),simplify=FALSE)
# blue<-rev(lapply(blue,function(x) x[[1]][[1]]))[-1]
# for(i in 1:length(blue)){
#   if(blue[[i]]==""){blue[[i]]<-names(blue)[[i]]}
#    if(blue[[i]]==""){blue[[i]]<-paste(rep(names(blue)[[i]],i),collapse="")}
# }
 return(blue)
}

