legendbar <- function(attribute, 
                      colPalette=NULL,
                      legendName="Legend",
                      bgc='white',
                      nticks=11, 
                      title='',
                      temp=FALSE) {
 
  pal<-colorRampPalette(c( "green", "orange","brown"), space = "Lab")
  if(is.null(colPalette) ){
    colPalette<-pal(min(10,length(attribute) ) ) }else{ xx<-colPalette<-as.character(substr(colPalette,1,7)) }


  
if(is.numeric(attribute)){


  min=min(attribute)
  max=max(attribute)
  
  if((length(attribute)-length(colPalette) == 1)) {
    ticks<-signif(attribute, 4)
    if(ticks[length(ticks)]==ticks[length(ticks)-1]){ ticks[length(ticks)]='NA'}
    min=1 #0
    max=length(ticks)
    ticks.at=attribute
    
  }else{
    
    df=data.frame(attribute,colPalette)
    df=df[ order(df[,1]), ]
    attribute=df[,1]
    colPalette=as.character(df[,2] )
    
    ticks<-signif(seq(min, max, len=nticks), 4)
    min=0
    max=length(ticks)
    ticks.at=seq(min,max, length.out =length(ticks))
  }
  
  
  
  }else{
   min=0
   max=length(factor(attribute))
   ticks=as.character(levels(factor(attribute)))
   xxx=seq(min,max, length.out =(length(ticks)+1) )
   ticks.at=sapply(2:(length(xxx)), function(i) mean(c(xxx[i],xxx[i-1])) )
 }
 


  
  lut=as.character(colPalette)
  if(length(lut)>1){
  scale = (length(lut))/(max-min)
  
#   ticks.at[1:(length(ticks.at)-1)]=sapply(1:(length(ticks.at)-1), function(i) (i-1)/scale + min )
  ticks.at=sapply(1:(length(ticks.at)), function(i) (i-1)/scale + min )
  
  if(is.factor(attribute) ) {
    dx=diff(ticks.at)[1]
    ticks.at[1:(length(ticks.at)-1)]=sapply(1:(length(ticks.at)-1), function(i) (i-1)/scale + min +dx*0.5 )  
  }
                      
  }else{scale=1}
  
  if(length(lut)>25){h=length(lut)*15}else{h=360}
  
  png(filename =ifelse(temp,paste(tempdir(),'/',legendName,'.png',sep=""),paste(legendName,'.png',sep="") ), height=h, width=108,pointsize = 8, units = "px", bg="white")
  par(mar=c(0,1,1,4),bg=bgc)
  #dev.new(width=1.75, height=5)
  plot( c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4,ticks.at, labels=ticks,las=1, lwd = 0, lwd.ticks =1) # , las=1
  for (i in 1:(length(lut))) {
    y = (i-1)/scale + min
    rect(0,y,1,y+1/scale, col=lut[i], border=NA)
  }  
  
  graph1 <- dev.cur()
  dev.off(graph1)
}
