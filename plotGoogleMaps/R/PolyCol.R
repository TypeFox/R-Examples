PolyCol <-
function(attribute,colPalette=NULL,at=NULL) {

pal<-colorRampPalette(c( "green", "orange","brown"), space = "Lab")
          
# create caracter from factor string by replacing one category i.e. key.entries with colors
reclassify = function(data, inCategories, outCategories)  { outCategories[ match(data, inCategories)]  }
          
if(!is.numeric(attribute)){ attribute<-as.factor( attribute)}
          
if(length(colPalette)==1) {
              x<- rep(colPalette,length(attribute))
              col.data<-list(cols=as.character(substr(x,1,7)),col.uniq=colPalette, 
                             att=ifelse(!is.factor(attribute),paste("[",min(attribute)," , ",max(attribute),"]",sep=""), " "),
                             brks=ifelse(!is.factor(attribute),paste("[",min(attribute)," , ",max(attribute),"]",sep=""), " "))
              return(col.data) }
          
if(is.null(colPalette) ){
      colPalette<-pal(min(10,length(attribute)) ) }else{ xx<-colPalette<-as.character(substr(colPalette,1,7)) }

if(is.null(at) & !is.factor(attribute)){
  numcolor = length(colPalette)
  bre<-quantile(attribute, seq(1,numcolor)/numcolor, na.rm=TRUE)
  breakss<-factor(c(min(attribute,na.rm=T),bre))
  break_unique<-as.numeric(levels(breakss))
  break_unique[length(break_unique)]<-max(attribute,na.rm=T)
  break_unique=unique(break_unique)
  at<-break_unique
    if(length(colPalette)>=length(break_unique)){
      print("using original PolyCol")
        colPalette<-colPalette[1:length(break_unique)] 
        } else{ colPalette<- as.character(substr(colPalette[1:length(break_unique)],1,7))
        }
}

if(is.factor(attribute)){
        
    if(length(colPalette)!=nlevels(attribute)) {
         xx<-colPalette<- as.character(substr(pal(nlevels(attribute)),1,7))    }                           
    x<- reclassify(attribute, inCategories=levels(attribute),  outCategories=colPalette)
    col.data<-list(cols=as.character(substr(x,1,7)),col.uniq=colPalette, att=levels(attribute),brks=levels(attribute) )
    return(col.data)
                              
                }else{
                          
                          if(length(colPalette)==length(attribute)){
                            min=min(attribute)
                            max=max(attribute)
                            df=data.frame(attribute,colPalette)
                            df=df[ order(df[,1]), ]
                            attribute=df[,1]
                            colP=df[,2][!is.na(df[,2])] 
                            lut=unique(as.character(colP) )

                            tt=lapply(lut, function(i){
                              x=max(df[df[,2]==i,1]  )
                            })
                            
                            break_unique<-unique(as.vector(do.call(cbind,tt)))
                            break_unique[length( break_unique)]=max(attribute,na.rm=T)
                            
                            break_unique=c(min(attribute,na.rm=T), break_unique)
                            
#                             break_unique=unique(break_unique)
                            
                            
                            col.data<-list(cols=as.character(substr(colPalette,1,7)),col.uniq=as.character(substr(lut,1,7)), 
                                           att=attribute, brks=break_unique )
                            return(col.data)
                            
                          } # end of Lcp=LAtt
                          
#                           if(is.null(at)){
#                                numcolor = length(colPalette) + 1
#                                bre<-quantile(attribute, seq(1,numcolor)/numcolor, na.rm=TRUE)
#                                breakss<-factor(c(min(attribute,na.rm=T),bre))
#                                break_unique<-as.numeric(levels(breakss))
#                                break_unique[length(break_unique)]<-max(attribute,na.rm=T)
#                                break_unique=unique(break_unique)
#                                       
#                                     if(length(colPalette)>=length(break_unique)){
#                                          colPalette<-colPalette[1:length(break_unique)] } else{ 
#                                            colPalette<- as.character(substr(colPalette[1:length(break_unique)-1],1,7))}
#                                                                       
#                                    }else{
                                      at[1]=min(attribute,na.rm=T)
                                      at[length(at)]=max(attribute,na.rm=T)
                                      break_unique=at
                                      colPalette<- as.character(substr(colPalette[1:length(break_unique)-1],1,7))
#                                    }   
                          

                            
                            
                      atr<-cut(attribute, break_unique ,include.lowest = TRUE, dig.lab=6)
#                                     x<-factor(atr,labels=colPalette[1:(length(break_unique)-1)] )
                      x<- reclassify(atr, inCategories=levels(atr),  outCategories=colPalette[1:(length(break_unique)-1)] )
                     col.data<-list(cols=as.character(substr(x,1,7)),col.uniq=colPalette, att=levels(atr), brks=break_unique )
                                    return(col.data)                               
                               
                               }
                                         
                                         
                                         }
