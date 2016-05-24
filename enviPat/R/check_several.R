check_several <-
function(pattern,dmz,ppm=TRUE){

    ############################################################################
    # (1) issue warnings #######################################################
    if(length(pattern[[1]])<2 || colnames(pattern[[1]])[c(1,2)]!=c("m/z","abundance")){stop("WARNING: pattern has invalid entries\n")}
    if(ppm==TRUE & dmz<0){stop("\n WARNING: ppm=TRUE -> dmz must be >0\n");}
    if(ppm!=TRUE & ppm!=FALSE){stop("WARNING: ppm TRUE or FALSE")}
    if(length(pattern)<2){stop("WARNING: nothing to compare")}
    ############################################################################
    # (2) check ################################################################    
    # (2.1) disassemble ######################################################## 
    getit1<-c();
    getit2<-c();
    getit3<-c();
    for(i in 1:length(pattern)){
        getit1<-c(getit1,rep(i,length(pattern[[i]][,1])))
        getit2<-c(getit2,pattern[[i]][,1]);
        getit3<-c(getit3,seq(1,length(pattern[[i]][,1]),1));
    }
    getit1<-getit1[order(getit2)];
    getit3<-getit3[order(getit2)];    
    getit2<-getit2[order(getit2)];
    # (2.2) ####################################################################
    get1<-rep("",length(pattern))
    get2<-rep("",length(pattern))
    get3<-rep("",length(pattern))
    for(i in 2:length(getit1)){
      if(ppm==FALSE){
        if(getit1[i-1]!=getit1[i]){
          if((getit2[i-1]+dmz)>=getit2[i]){
             get1[getit1[i]]<-"TRUE";
             get1[getit1[i-1]]<-"TRUE";
             #get2[getit1[i]]<-paste(get2[getit1[i]],names(pattern)[getit1[i-1]],"/",sep="")
             #get2[getit1[i-1]]<-paste(get2[getit1[i-1]],names(pattern)[getit1[i]],"/",sep="")
             get2[getit1[i]]<-paste(get2[getit1[i]],getit1[i-1],"/",sep="")
             get2[getit1[i-1]]<-paste(get2[getit1[i-1]],getit1[i],"/",sep="")
             get3[getit1[i]]<-paste(get3[getit1[i]],getit3[i-1],"/",sep="")
             get3[getit1[i-1]]<-paste(get3[getit1[i-1]],getit3[i],"/",sep="")        
          }
        }
      }else{ # ppm==TRUE
        if(getit1[i=1]!=getit1[i]){
          if((getit2[i-1]+(getit2[i-1]*dmz/1e6))>=getit2[i]){
             get1[getit1[i]]<-"TRUE";
             get1[getit1[i-1]]<-"TRUE";
             #get2[getit1[i]]<-paste(get2[getit1[i]],names(pattern)[getit1[i-1]],"/",sep="")
             #get2[getit1[i-1]]<-paste(get2[getit1[i-1]],names(pattern)[getit1[i]],"/",sep="")
             get2[getit1[i]]<-paste(get2[getit1[i]],getit1[i-1],"/",sep="")
             get2[getit1[i-1]]<-paste(get2[getit1[i-1]],getit1[i],"/",sep="")
             get3[getit1[i]]<-paste(get3[getit1[i]],getit3[i-1],"/",sep="")
             get3[getit1[i-1]]<-paste(get3[getit1[i-1]],getit3[i],"/",sep="")        
          }
        }
      }
    }
    if(any(get1!="")){cat("\n Overlaps detected!\n\n")}
    checked<-data.frame(names(pattern),get1,get2,get3)
    names(checked)<-c("compound","warning","to?","peak#")
    return(checked)
    ############################################################################

}
