# eloplot 14_10_18

#########################################################################################
# PLOT FUNCTION
eloplot <- function(eloobject, ids="all", interpolate="yes", from="start", to="end", color=TRUE){
  res <- eloobject
  # plotdata handling
  if(interpolate=="yes"){
    plotdata=res$cmat
  }else{
    plotdata=res$lmat
  }
  
  # exclude IDs that had interactions only on one day...
  temp <- rbind(rowSums(table(res$logtable$winner, res$logtable$Date) > 0)[res$allids], rowSums(table(res$logtable$loser, res$logtable$Date) > 0)[res$allids])
  colnames(temp) <- res$allids
  if(1 %in% colSums(temp, na.rm=T)) {
    good <- colnames(temp)[colSums(temp, na.rm=T) > 1]
    bad <- colnames(temp)[colSums(temp, na.rm=T) == 1]
    plotdata <- plotdata[, good]
    
  }
  
  
  # ids handling
  if(ids[1] %in% c("random.20", "first.20", "all")){
    if(ids[1]=="random.20"){
      if(length(colnames(plotdata))>20){
        ids=sample(colnames(plotdata),20)
      }else{
        ids=colnames(plotdata)
      }
    }
    if(ids[1]=="first.20"){
      if(length(colnames(plotdata))>20){
        ids=colnames(plotdata)[1:20]
      }else{
        ids=colnames(plotdata)
      }
    }
    if(ids[1]=="all"){
      ids=colnames(plotdata)
      if(exists("bad")) { 
        warning("IDs for which interactions were observed on only one day were excluded", call. = FALSE)
        rm(bad)
      }
    }
  }else{
    ids=ids
  }
  
  if(exists("bad")) { 
    if(length(intersect(ids, bad)) >= 1) { 
      warning("IDs for which interactions were observed on only one day were excluded", call. = FALSE) 
      ids <- intersect(ids, good)
    }
  }
  
  plotdata=plotdata[,ids]
  
  # data.range handling
  if(from=="start" & to=="end"){
    dates=seq(min(res$truedates),max(res$truedates),"day")
    ids.wo=""
  }
  if(from!="start" | to!="end"){
    if(from!="start" & to=="end"){
      dates=seq(as.Date(from),max(res$truedates),"day")
    }
    if(from=="start" & to!="end"){
      dates=seq(min(res$truedates), as.Date(to),"day")
    }
    if(from!="start" & to!="end"){
      dates=seq(as.Date(from),as.Date(to),"day")
    }
    plotdata=plotdata[which(res$truedates %in% dates), ]
    xx=apply(plotdata, 2, function(x) { sum(is.na(x)) } )
    plotdata=plotdata[, which(xx < nrow(plotdata))]
    ids=colnames(plotdata)
    ids.wo=names(xx)[xx == nrow(plotdata)]
  }
  
  
  
  fst.month=unique(as.Date(as.yearmon(dates)))
  if(dates[1]>fst.month[1]){fst.month=fst.month[-1]}
  if(tail(dates, 1) < tail(fst.month, 1)) { fst.month=fst.month[-(length(fst.month))] }
  labs=c(min(dates), fst.month[length(fst.month)%/%3], fst.month[(length(fst.month)%/%3)*2], max(dates))
  if(labs[1] == labs[2]) { labs <- labs[-1] }
  ats=which(dates%in%labs)
  
  # specify colors
  colo=colors()[c(552,254,652,26,33,259,32,610,51,148,31,47,128,7,8,12,24,53,56,68,547,116,142,30,204,498,22,62,146)][1:length(ids)]
  
  if(color){
    layout(matrix(c(1,2), ncol=2), heights=c(5,5), widths=c(4,1))
    par(mar=c(5,4,4,0))   #b5 l4 t4 r2
    plot(1:length(dates),1:length(dates),ylim=range(plotdata, na.rm=T), ylab="Elo-rating", xlab="date", type="n", xaxt="n", las=1, cex.axis=0.8, cex.lab=0.8)
    axis(1, at = ats, labels = labs, cex.axis=0.8, las=1)
    mtext(c("first day", "last day"), side = 1, line=2, at = c(1,tail(ats,1)), cex=0.8)
    for (i in 1:length(ids)){
      lines(plotdata[,ids[i]], lty=1, type="l", col=colo[i])
    }
    par(mar=c(5,1,3.8,1))   #b5 l4 t4 r2
    plot(1:2,1:2, xaxt="n", yaxt="n", type="n", bty="n", ylab="", xlab="")
    if(length(c(ids,ids.wo))>23){
      legend(1,2.04, c(c(ids,ids.wo)[1:23],"..."), cex=0.8, bty = "n", col=c(colo[1:23],rep("white",length(1:length(c(ids,ids.wo))))), pch=20, pt.cex=1.5)
    }else{
      legend(1,2.04, c(ids,ids.wo), cex=0.8, bty = "n", col=c(colo[1:length(ids)],rep("white",length(ids.wo))), pch=20, pt.cex=1.5)
    }
  }else{
    p.dates=which(dates%in%pretty(dates,15))
    pointsdata= plotdata[p.dates,]
    if(length(c(ids,ids.wo))/25>1){
      p.times=rep(1:25,round(length(c(ids,ids.wo))/25))
    }else{
      p.times=1:25
    }
    layout(matrix(c(1,2), ncol=2), heights=c(5,5), widths=c(4,1))
    par(mar=c(5,4,4,0))   #b5 l4 t4 r2
    plot(1:length(dates),1:length(dates),ylim=range(plotdata, na.rm=T), ylab="Elo-rating", xlab="date", type="n", xaxt="n", las=1, cex.axis=0.8, cex.lab=0.8)
    axis(1, at = ats, labels = labs, cex.axis=0.8, las=1)
    mtext(c("first day", "last day"), side = 1, line=2, at = c(1,tail(ats,1)), cex=0.8)
    for (i in 1:length(ids)){
      lines(plotdata[,ids[i]], lty=i)
      points(p.dates,pointsdata[,ids[i]], pch=p.times[i], cex=0.8, bg="grey")
    }
    par(mar=c(5,1,3.8,1))   #b5 l4 t4 r2
    plot(1:2,1:2, xaxt="n", yaxt="n", type="n", bty="n", ylab="", xlab="")
    if(length(c(ids,ids.wo))>23){
      legend(1,2.04, c(c(ids,ids.wo)[1:23],"..."), cex=0.8, bty = "n", pch=p.times[c(1:23,(length(p.times)+1))], pt.cex=1, pt.bg="grey")
    }else{
      legend(1,2.04, c(ids,ids.wo), cex=0.8, bty = "n", pch=p.times[c(1:length(ids),26)], pt.cex=1, pt.bg="grey")
    }
  }
}
