plotrd<-function(clustrdOut,what=c("all","all"),obslabel=0,attlabel=0,density=T,fname=0){
  x1=x2=NULL
  ..level..=obslab=itlab=NULL
  group=NULL
  df=data.frame(x1=clustrdOut$obscoord[,1],x2=clustrdOut$obscoord[,2],group=factor(clustrdOut$cluID))
  dfAtt=data.frame(x1=clustrdOut$attcoord[,1],x2=clustrdOut$attcoord[,2])
  dfG=data.frame(x1=clustrdOut$centroid[,1],x2=clustrdOut$centroid[,2])
  
  if(obslabel[1]==0){
    obslabel=paste(" ")
  }
  df$obslab=obslabel
  
  if(attlabel[1]==0){
    attlabel=paste("v.",1:nrow(dfAtt),sep="")
  }
  dfAtt$itlab=attlabel

  if (density==T) {
    dmap=ggplot(df,aes(x=x1,y=x2))+stat_density2d(aes(fill=..level..), geom="polygon") + scale_fill_gradient(low="lightyellow", high="blue")+ theme_bw()
  }
  else {
    dmap=ggplot(df,aes(x=x1,y=x2)) + theme_bw()
  }
  
  
  if (what[1] == "all" & what[2] == "all") {
    dmap=dmap + xlim(min(df$x1,dfAtt$x1,dfG$x1),max(df$x1,dfAtt$x1,dfG$x1))
    dmap=dmap + ylim(min(df$x2,dfAtt$x2,dfG$x2),max(df$x2,dfAtt$x2,dfG$x2))
    dmap=dmap + geom_text(data=df,aes(x=x1,y=x2,label=obslab),colour="black")
    dmap=dmap + geom_text(data=dfAtt,aes(x=x1,y=x2,label=itlab),colour="darkgreen")
    dmap= dmap+ geom_point(data=df,aes(x=x1,y=x2,colour=group),alpha=.25)
    dmap=dmap + geom_point(data=dfG,aes(x=x1,y=x2),colour="red",size=5)
    
  }
  else if (what[1] == "none" & what[2] == "all") {
    dmap=dmap + xlim(min(dfAtt$x1,dfG$x1),max(dfAtt$x1,dfG$x1))
    dmap=dmap + ylim(min(dfAtt$x2,dfG$x2),max(dfAtt$x2,dfG$x2))
    dmap=dmap + geom_text(data=dfAtt,aes(x=x1,y=x2,label=itlab),colour="darkgreen")
    dmap=dmap + geom_point(data=dfG,aes(x=x1,y=x2),colour="red",size=5)
  }
  else if (what[1] == "all" & what[2] == "none") {
    dmap=dmap + xlim(min(df$x1,dfG$x1),max(df$x1,dfG$x1))
    dmap=dmap + ylim(min(df$x2,dfG$x2),max(df$x2,dfG$x2))
    dmap=dmap + geom_text(data=df,aes(x=x1,y=x2,colour=group,label=obslab))
    dmap=dmap + geom_point(data=dfG,aes(x=x1,y=x2),colour="red",size=5)
  }  
  dmap= dmap+ theme(legend.position="none")
  dmap=dmap+xlab("")+ ylab("")
  dmap=dmap+geom_vline(xintercept=c(0,0))+geom_hline(yintercept=c(0,0))
  dmap= dmap+ theme(panel.grid.major = element_line(colour = 'black', linetype = 'dashed',size=.1), panel.grid.minor = element_line(colour = "black", linetype = 'dashed',size=.1)) 
  dmap= dmap+ theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  if(fname==0){
    fname=paste("figure",".png",sep="")
    ggsave(dmap,filename=fname,dpi=300)
  }else{
    fname=paste(fname,".png",sep="")
    ggsave(dmap,filename=fname,dpi=300)
  }
  dmap
}