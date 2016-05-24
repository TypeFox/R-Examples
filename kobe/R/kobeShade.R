setMethod('kobeShade', signature(object='numeric'),
          function(object,breaks=c(-0.1,50,60,70,80,90,100),
                     kobeShades=c("","\\cellcolor{gray90}","\\cellcolor{gray80}","\\cellcolor{gray70}","\\cellcolor{gray60}","\\cellcolor{gray50}"),
                     pct="\\%",...){

  #Kobe II strategy matrices to be prepared by the SCRS should highlight in a similar format as
  #shown in Annex Table 2 a progression of probabilities over 50 % and in the range of 50-59 %, 60-
  #69 %, 70-79 %, 80-89 % and ??? 90 %.
  object=pmin(pmax(object,0),1)*100         
  res   =data.frame("order"=seq(length(object)),object=round(object),"level"=cut(object,breaks))
  gry   =data.frame(level=attributes(unique(res$level))$levels,kobeShades)
  res  =merge(res,gry,all.x=TRUE)
  res  =res[order(res$order),]
  
  res$object=paste(ac(round(res$object)),pct,sep="")
  
  res=with(res,paste(kobeShades,object,sep=" "))
  
  return(res)})

setMethod('kobeShade', signature(object='data.frame'),
          function(object,breaks =c(-0.1,50,60,70,80,90,100),
                   kobeShades=c("","\\cellcolor{gray90}","\\cellcolor{gray80}","\\cellcolor{gray70}","\\cellcolor{gray60}","\\cellcolor{gray50}"),
                   pct="\\%",...){

     as.data.frame(apply(object,2,kobeShade,breaks=breaks,kobeShades=kobeShades,pct=pct))
     })

# setMethod('kobeShade', signature(object='cast_df'),
#           function(object,breaks =c(-0.1,50,60,70,80,90,100),
#                    kobeShades=c("","\\cellcolor{gray90}","\\cellcolor{gray80}","\\cellcolor{gray70}","\\cellcolor{gray60}","\\cellcolor{gray50}"),
#                    pct="\\%",...){
#             
#             as.data.frame(apply(object,2,kobeShade,breaks=breaks,kobeShades=kobeShades,pct=pct))})
setMethod('kobeShade', signature(object='matrix'),
          function(object,breaks =c(-0.1,50,60,70,80,90,100),
                   kobeShades=c("","\\cellcolor{gray90}","\\cellcolor{gray80}","\\cellcolor{gray70}","\\cellcolor{gray60}","\\cellcolor{gray50}"),
                   pct="\\%",...){

     apply(object,2,kobeShade,breaks=breaks,kobeShades=kobeShades,pct=pct)})

setMethod('kobe2sm', signature(object='data.frame'),
          function(object,cex   =1.0,
                         image  =list(levels=seq(0.0,1.0,0.05),
                                      col   =c(colorRampPalette(c("red4","red"))(12),colorRampPalette(c("yellowgreen","darkgreen"))(8))),
                         contour=list(levels=c(.6,.7,1.0,.9),
                                      col   =c("black"))){
                 
    nms   =dimnames(object)[[2]]
    nPlots=length(nms[nms %in% c("overFishing","overFished","green")])
                 
    ops<-par(mfrow=c(nPlots,1), mex=.5,mai=c( 0.5, 0.75 ,.750, 0.1),cex=par()$cex)
    
    res=list()
    if ("overFishing" %in% nms){
       x=transform(object, NotOverFishing=as.numeric(1-object[,"overFishing"]))[,c(nms[1:2],"NotOverFishing")]
       res[["F"]]    =kobe2smFn(x, image=image,contour=contour)
       mtext(expression(plain(P) (F<=F[MSY])),line=3, cex=cex)
       }
    if ("overFished" %in% nms){
       x=transform(object, NotOverFished=as.numeric(1-object[,"overFished"]))[,c(nms[1:2],"NotOverFished")]
       res[["SSB"]]  =kobe2smFn(x, image=image,contour=contour)
       mtext(expression(plain(P) (SSB>=SSB[MSY])),line=3, cex=cex)
       }
    if ("green" %in% nms){
       res[["Joint"]]=kobe2smFn(cbind(object[,c(1:2)],object[,"green"]), image=image,contour=contour)
       mtext(expression(plain(P) (F<=F[MSY]) %*% plain(P)(SSB>=SSB[MSY])),line=3, cex=cex)
       }

    par(mfrow=ops$mfrow,mex=ops$mex,mai=ops$mai,cex=ops$cex)

    invisible(res)})          
