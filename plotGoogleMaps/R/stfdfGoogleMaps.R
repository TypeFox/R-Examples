stfdfGoogleMaps<-function(stfdf,
                          zcol=1,
                          filename='',
                          layerName="",
                          plotNames=row.names(stfdf@sp),
                          aggregateFUN='mean',
                          round.att=1,
                          plot.height=300,
                          plot.width=300,
                           ...) {
  
  nameOfSP<-sapply(as.list(substitute({stfdf})[-1]), deparse)
  nameOfSP<-gsub("\\s","", nameOfSP)
  nameOfSP<-gsub('[!,",#,$,%,&,(,),*,+,-,.,/,:,;,<,=,>,?,@,_,^,`,|,~]', "x", nameOfSP)
  nameOfSP<-gsub('[[]', "X", nameOfSP)
  nameOfSP<-gsub('[]]', "X", nameOfSP)
  temporary = FALSE 
  if(filename==""){
    filename <- tempfile("map", fileext = c(".html"))
    plotdir=tempdir()
    temporary = TRUE
  }else{
    plotdir=paste('plot',nameOfSP,sep="")
    dir.create(plotdir)
  }
  
  if(layerName==""){
    layerName <- nameOfSP
  }

if(any('data'==slotNames(stfdf)) ){
  attribute=stfdf@data[,zcol] 
  for(i in 1:length(stfdf@data)) {
    if( identical(attribute,stfdf@data[,i])){
      attributeName<-names(stfdf@data)[i]  }
  }
}

if(!is.numeric(stfdf@data[,zcol]) ){
  stop("attribute must be numeric, change zcol in stfdfGoogleMaps ") 
}

SP= as(stfdf[,,zcol],'Spatial' )
SP$plot=rep(NA,length(SP@data[,1]))
wd=getwd()
setwd(plotdir)

plotdir <- ifelse(temporary,'',paste(plotdir,'/',sep=''))


for (i in 1:length(SP@data[,1])) {
  ff= paste(attributeName,i,'.png',sep="")
  png(filename =ff, width = plot.width, height = plot.height,  units = "px", pointsize = 10)
  plot(stfdf[i,,zcol], main=plotNames[i],ylab=attributeName)
  dev.off()
  SP$plot[i]=paste("<img src='",plotdir,ff,"'  />", sep="")
  
}
setwd(wd)

xxx = as(stfdf[,,zcol], "xts")
SP@data[,paste(attributeName,aggregateFUN,sep='') ]=sapply(xxx, FUN = aggregateFUN, na.rm=TRUE)

if(is.numeric(round.att)){
  SP@data[,paste(attributeName,aggregateFUN,sep='') ] <- round(SP@data[,paste(attributeName,aggregateFUN,sep='') ], round.att)
}
SP <- SP [, c(paste(attributeName,aggregateFUN,sep=''), 'plot')]
m <- plotGoogleMaps(SP,
                 zcol=paste(attributeName,aggregateFUN,sep=''),
                  filename=ifelse(temporary,"",filename),
                  layerName=layerName, ...)
return(m)
}
