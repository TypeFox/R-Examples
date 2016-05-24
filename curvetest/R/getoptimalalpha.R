getoptimalalpha <-
function(formula, data, plotit=F){
 ##To get the optimal alpha value by general cross validations. intensive. 
    require(locfit)
    gg<-gcvplot(formula,data=data,alpha=seq(0.1,0.91,by=0.005)) 
    dd<-data.frame(alpha=as.numeric(gg$alpha),   gcv=gg$values)
    wh<-which.min(gg$values)
    if(plotit) plot(dd)
    return(gg$alpha[wh]) 
 }
