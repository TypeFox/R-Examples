discretVar <-
function(fixBrks,listBrks,pt,var,nclass,
                     style,palCol,diverg,divergBrk,palColPos,
                     palColNeg,NACol,lgdRnd){
  if(fixBrks==TRUE){
    nclass<-length(listBrks)-1
    distr<-classIntervals(pt[is.na(pt[,var])==FALSE,var],nclass,style="fixed",fixedBreaks=listBrks)$brks
  }else{distr<-classIntervals(pt[is.na(pt[,var])==FALSE,var],nclass,style=style)$brks}

  # diverging management
  if (diverg) {
    nb.pos <- sum(distr > divergBrk)
    nb.neg <- sum(distr < divergBrk)
    if (nb.pos > 0) {
      if (nb.pos < 3) 
        palpos <- brewer.pal(3, palColPos)[1:nb.pos]
      else palpos <- brewer.pal(nb.pos, palColPos)
      palette <- palpos
    }


    if (nb.neg > 0) {
      if (nb.neg < 3) 
        palneg <- brewer.pal(3, palColNeg)[1:nb.neg]
      else palneg <- brewer.pal(nb.neg, palColNeg)
      
    }
    if (nb.neg+nb.pos > nclass){
    colours <- c(rev(palneg[2:nb.neg]), "#F5F5F5",palette[2:nb.pos])
    }else {colours <- c(rev(palneg),palette)}
    
  }else{colours<-brewer.pal(nclass,palCol)}

  # color assignement
  pt$col<-colours[(findInterval(pt[,var],distr,all.inside=TRUE))]
  pt[is.na(pt$col)==TRUE,"col"]<-NACol

  if ((length(pt[is.na(pt[,var]),var])>0)==TRUE){
    pdd<-1
  }else{pdd<-0}
  lblLeg<-myLeg(distr,lgdRnd)
return(list(pt,colours,pdd,lblLeg))
}  
