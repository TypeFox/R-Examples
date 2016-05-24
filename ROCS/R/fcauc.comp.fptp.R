fcauc.comp.fptp <-
function(roc.1, roc.2, FDR.cut=0.2, lwd=3, colors=c("blue", "green","cyan", "red"), exp.labels=c("experiment 1", "experiment 2"))
{
    a<-roc.1
    b<-roc.2
    
    o<-order(-a$FP, -a$TP)
    a$FP<-a$FP[o]
    a$TP<-a$TP[o]
    a$TDR<-a$TDR[o]
    
    o<-order(-b$FP, -b$TP)
    b$FP<-b$FP[o]
    b$TP<-b$TP[o]
    b$TDR<-b$TDR[o]
    
    cut.a<-min(which(a$TDR>=1-FDR.cut))
    cut.b<-min(which(b$TDR>=1-FDR.cut))
    plot(a$FP, a$TP, type="n", xlab="FPR", ylab="TPR")
    
    lines(a$FP[1:cut.a], a$TP[1:cut.a], lwd=lwd, col=colors[1])
    lines(a$FP[cut.a:length(a$FP)], a$TP[cut.a:length(a$TP)], lwd=lwd, col=colors[2])
    
    lines(b$FP[1:cut.b], b$TP[1:cut.b], lwd=lwd, col=colors[3])
    lines(b$FP[cut.b:length(b$FP)], b$TP[cut.b:length(b$TP)], lwd=lwd, col=colors[4])
    
    lines(c(0.4, 0.6), c(0.3, 0.3), lwd=3, col=colors[2])
    lines(c(0.6, 0.8), c(0.3, 0.3), lwd=3, col=colors[1])
    text(0.8, 0.3, exp.labels[1],pos=4)
    
    lines(c(0.4, 0.6), c(0.2, 0.2), lwd=3, col=colors[4])
    lines(c(0.6, 0.8), c(0.2, 0.2), lwd=3, col=colors[3])
    text(0.8, 0.2, exp.labels[2],pos=4)
    
    text(0.6, 0.4, "FDR<=0.2",pos=2)
    text(0.6, 0.4, "FDR>0.2",pos=4)
    
}
