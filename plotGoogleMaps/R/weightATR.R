weightATR <-
function(attribute,strokeWeight=1) {
      # attribute=soil.ll@data$ID
        if(length(strokeWeight)==1){ x<-rep(strokeWeight,length(attribute)) 
                                         return(as.numeric(x))}
                
        if(is.factor(attribute)){
        
                        if(length(strokeWeight)!=nlevels(attribute)) {
                               stop("length of strokeWeight should match number of factor levels")   }
                              
                              }else{
                               bre<-quantile(attribute, seq(1,length(strokeWeight))/length(strokeWeight))
                                breakss<-factor(c(min(attribute),bre))
                                break_unique<-as.numeric(levels(breakss))
                                     if(length(strokeWeight)>=length(break_unique)){
                                         strokeWeight<-strokeWeight[1:length(break_unique)] } else{
                                                                     strokeWeight<- strokeWeight[1:length(break_unique)-1] } 
                                
                          if(sum(as.numeric(levels(factor(attribute))))-sum(break_unique) ==0 ){
                             attribute<-factor(attribute)}   
                                                                      
                                                                      }
                
                
if (length(strokeWeight)>1){
   if(!is.numeric(attribute) && !is.character(attribute) ) {

              if(nlevels(attribute)== length(strokeWeight)){
                x<-factor(attribute,labels=strokeWeight)
                return(as.numeric(x))
                }else{
                stop("length of strokeWeight should match number of factor levels")
                }
   }else{

    if( is.numeric(attribute)){
     x<-factor(cut(attribute, break_unique ,include.lowest = TRUE, dig.lab=6),labels=strokeWeight[1:length(break_unique)-1] )
    return(as.numeric(x))
     }
   }

 if(nlevels(factor(attribute))== length(strokeWeight)){
 x<-factor(factor(attribute),labels=strokeWeight)
         return(as.numeric(x))
         }else{
          stop("length of strokeWeight should match number of factor levels")
         }

                          }

                                  
                                  }
