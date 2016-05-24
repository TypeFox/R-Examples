bb$winP <- bb$W/bb$G
bb$predWinP <- predict(glm.bb,
    newdata=data.frame(runmargin=bb$runmargin),type='response')
bb$winPdiff <- bb$winP - bb$predWinP
bb[rev(order(abs(bb$winPdiff)))[1:5],c(1,22:24)]
bb.plot1 <- xyplot(winP~predWinP,data=bb,
    panel=function(x,y,...) {
        panel.xyplot(x,y,...)
        panel.abline(0,1)
    })

rm <- seq(-5,5,by=0.1)
wp <- predict(glm.bb, newdata=data.frame(runmargin=rm),type='response')
bb.plot2 <- xyplot(winP~runmargin, data=bb,xlim=c(-2.5,2.5), ylim=c(0,1),
    panel = function(x,y,...){
        panel.xyplot(x,y,...)
        panel.xyplot(rm,wp,type='l',col='gray50')
    })
