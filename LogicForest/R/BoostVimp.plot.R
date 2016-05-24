BoostVimp.plot <-
function(fit, num=10, pred, norm=TRUE, titles=TRUE)
{
 if (missing(pred)) {pred<-fit$PredImp}
 if (pred==TRUE) {
    tpred<-sort(fit$Pred.import, decreasing=TRUE)[1:num]
    if (norm==TRUE){
      tpred.norm<-tpred/tpred[1]
      up1<-1.1
      }
    else {tpred.norm<-tpred; up1<-tpred[1]+0.5}
    if (titles==TRUE) {title<-c("Predictor Importance", "Importance Scores")}
    else {title<-c("","")}
    par(mfrow=c(1,1), cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1, ask=TRUE)
    plot(x=tpred.norm, y=c(num:1), xlab=title[2], ylab="", xlim=c(0.025, up1),
            ylim=c(0.5, num+.5), yaxt="n", main=title[1])
    segments(0, num, tpred.norm[1], num, lty=1)
    for (i in 1:(num-1)) {segments(0,(num-i),tpred.norm[i+1],(num-i),lty=1)} 
    axis(2, at=c(num:1), labels=names(tpred))
    plot.new()
    }
 if (fit$PIimp=="AddRemove")
    {
    tPIs<-sort(fit$AddRemove.PIimport, decreasing=T)[1:num]
    if (norm==TRUE) {
      tPIs.norm<-tPIs/tPIs[1]
      up2<-1.1
      if (titles==TRUE) {titles<-c("Normalized Scores", "Add/Remove importance")}
      else {titles<-c("","")}  
      }
    else {
      tPIs.norm<-tPIs
      up2<-tPIs[1]+.5
      if (titles==TRUE) {titles<-c("Importance Scores", "Add/Remove Importance")}
      else {titles<-c("","")}
      }
    par(cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1)
    plot(x<-tPIs.norm, y<-c(num:1), xlab=titles[1], ylab="", xlim=c(0.025,up2),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[2])
    segments(0, num, tPIs.norm[1], num, lty=1)
    for (i in 1:(num-1)) {segments(0,(num-i),tPIs.norm[i+1],(num-i),lty=1)} 
    axis(2, at=c(num:1), labels=names(tPIs))
    }
  if (fit$PIimp=="Permutation")
    {
    tPIs<-sort(fit$Perm.PIimport, decreasing=T)[1:num]
    if (norm==TRUE) {
      tPIs.norm<-tPIs/tPIs[1]
      up2<-1.1
      if (titles==TRUE) {titles<-c("Normalized Scores", "Permutation importance")}
      else {titles<-c("","")}  
      }
    else {
      tPIs.norm<-tPIs
      up2<-tPIs[1]+.5
      if (titles==TRUE) {titles<-c("Importance Scores", "Permutation Importance")}
      else {titles<-c("","")}
      }
    par(cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1)
    plot(x<-tPIs.norm, y<-c(num:1), xlab=titles[1], ylab="", xlim=c(0.025,up2),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[2])
    segments(0, num, tPIs.norm[1], num, lty=1)
    for (i in 1:(num-1)) {segments(0,(num-i),tPIs.norm[i+1],(num-i),lty=1)} 
    axis(2, at=c(num:1), labels=names(tPIs))
    }
  if (fit$PIimp=="Both")
    {
    tPIs1<-sort(fit$Perm.PIimport, decreasing=T)[1:num]
    tPIs2<-sort(fit$AddRemove.PIimport, decreasing=T)[1:num]
    if (norm==TRUE) {
      tPIs.norm1<-tPIs1/tPIs1[1]
      tPIs.norm2<-tPIs2/tPIs2[1]
      up2a<-1.1
      up2b<-1.1
      if (titles==TRUE) {titles<-c("Normalized Scores", "Permutation importance", "Add/Remove importance")}
      else {titles<-c("","","")}  
      }
    else {
      tPIs.norm1<-tPIs1
      up2a<-tPIs1[1]+.5
      tPIs.norm2<-tPIs2
      up2b<-tPIs2[1]+.5
      if (titles==TRUE) {titles<-c("Importance Scores", "Permutation importance", "Add/Remove importance")}
      else {titles<-c("","","")}
      }
    par(mfrow=c(1,2), cex.axis=0.6, cex.main=0.9, mai=c(1,1.4,.5,.1), pch=19, las=1)
    plot(x=tPIs.norm1, y=c(num:1), xlab=titles[1], ylab="", xlim=c(0.025, up2a),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[2])
    segments(0, num, tPIs.norm1[1], num, lty=1)
    for (i in 1:(num-1)) {segments(0,(num-i),tPIs.norm1[i+1],(num-i),lty=1)} 
    axis(2, at=c(num:1), labels=names(tPIs1))
    plot(x<-tPIs.norm2, y<-c(num:1), xlab=titles[1], ylab="", xlim=c(0.025,up2b),
              ylim=c(0.5, num+.5), yaxt="n", main=titles[3])
    segments(0, num, tPIs.norm2[1], num, lty=1)
    for (i in 1:(num-1)) {segments(0,(num-i),tPIs.norm2[i+1],(num-i),lty=1)} 
    axis(2, at=c(num:1), labels=names(tPIs2))
    }
}
