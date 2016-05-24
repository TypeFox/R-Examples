#' make the bins with equal number of observations  
#' or using user-defined breaks
#'
#' @param COV.data numeric vector that need to make bins
#' @param N.covbin the number of bins
#' @param breaks.data user-defined breaks
#' @return information of the binning with summary
#' @export
#' @author Eun-Kyung Lee \email{lee.eunk@@gmail.com}
#' @examples
#' data(origdata)
#' makeCOVbin(origdata$TIME,7)


makeCOVbin<-function(COV.data,N.covbin=NULL,breaks.data=NULL){   
   find.LU<-function(bin.ID){
      LU<-NULL
      for(i in 1:length(bin.ID)){
         LU<-rbind(LU,as.numeric(unlist(strsplit(unlist(strsplit(
                   unlist(strsplit(unlist(strsplit(
                   unlist(strsplit(as.character(bin.ID[i]),",")),"]")),
                                 '\\[')),"\\(")),"\\("))))    
      }
      colnames(LU)<-c("Lower","Upper")
      return(LU) 
   }

   data.temp<-data.frame(COV.data=COV.data)
   if(!is.null(breaks.data)){
      range.temp<-range(COV.data)
      if(min(breaks.data)>range.temp[1]){
         breaks.data[1]<-range.temp[1]-(range.temp[2]-range.temp[1])*0.1
      } else if(max(breaks.data)<range.temp[2]){
         breaks.data[length(breaks.data)]<-range.temp[2]+
                                         (range.temp[2]-range.temp[1])*0.1      
      }
      cut.temp<-cut(COV.data, breaks=breaks.data) 
      tab<-ddply(data.temp,.(cut.temp), summarize,
                       mid.COV=round(mean(COV.data, na.rm=T),2),.drop=FALSE)
        
   } else{
      if(N.covbin<length(table(COV.data))){
         cutpoints<-quantile(COV.data,(0:N.covbin)/N.covbin)
         temp.id<-which(diff(cutpoints)==0)
         if(length(temp.id)!=0)
            cutpoints<-cutpoints[-temp.id]
         cut.temp<-cut(COV.data,cutpoints,include.lowest=TRUE)
         if(sum(table(cut.temp)==0)!=0){
            temp.id<-which(table(cut.temp)==0)
            cutpoints<-cutpoints[-temp.id]
            cut.temp<-cut(COV.data,cutpoints,include.lowest=TRUE)
         }
         tab<-ddply(data.temp,.(cut.temp),summarize, 
                      med.COV=round(median(COV.data, na.rm=T),2),.drop=FALSE)
      } else{
         cut.temp.id<-as.numeric(names(table(COV.data)))
         cut.temp.id[1]<-cut.temp.id[1]-0.1
         cut.temp.id[length(cut.temp.id)]<-cut.temp.id[length(cut.temp.id)]+0.1
         cut.temp<-cut(COV.data,cut.temp.id)     
         tab<-ddply(data.temp,.(cut.temp),summarize,
                   med.COV=round(median(COV.data, na.rm=T),2),.drop=FALSE)      
      }
   } 
   LU.temp<-find.LU(tab[,1])
   colnames(LU.temp)<-c("lower.COV","upper.COV")
   tab<-data.frame(tab,n.bin=c(table(cut.temp)),LU.temp)
   return(list(COV.bin=cut.temp,COV.bin.summary=tab))
}

