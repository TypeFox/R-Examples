MarkerCat <-
  function(attribute,colPalette=NULL) {
    # attribute=soil.ll@data$ID
    pal<-colorRampPalette(c( "green", "orange","brown"), space = "Lab")
    
    # create caracter from factor string by replacing one category i.e. key.entries with colors
    reclassify = function(data, inCategories, outCategories)  { outCategories[ match(data, inCategories)]  }
    
    if(!is.numeric(attribute)){ attribute<-as.factor( attribute)}
    
    if(length(colPalette)==1) {
      x<- rep(colPalette,length(attribute))
      col.data<-list(cols=x,col.uniq=colPalette, att=ifelse(!is.factor(attribute),paste("[",min(attribute)," , ",max(attribute),"]",sep=""), " "))
      return(col.data) }
    
    if(is.factor(attribute)){
      
      if(length(colPalette)!=nlevels(attribute)) {
        xx<-colPalette   }
      
      x<- reclassify(attribute, inCategories=levels(attribute),  outCategories=colPalette)
      
      
      
      col.data<-list(cols=x,col.uniq=colPalette, att=levels(attribute) )
      return(col.data)
      
    }else{
      bre<-quantile(attribute, seq(1,length(colPalette))/length(colPalette))
      breakss<-factor(c(min(attribute),bre))
      break_unique<-as.numeric(levels(breakss))
      
      if(length(colPalette)>=length(break_unique)){
        colPalette<-colPalette[1:length(break_unique)] } else{
          colPalette<-colPalette[1:length(break_unique)-1] }
      
      atr<-cut(attribute, break_unique ,include.lowest = TRUE, dig.lab=6)                                 
      #                                     x<-factor(atr,labels=colPalette[1:(length(break_unique)-1)] )
      x<- reclassify(atr, inCategories=levels(atr),  outCategories=colPalette[1:(length(break_unique)-1)] )
      col.data<-list(cols=x,col.uniq=colPalette, att=levels(atr) )
      return(col.data)                               
      
    }
    
    
  }
