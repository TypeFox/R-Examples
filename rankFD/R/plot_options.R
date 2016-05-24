#################################################################
#Grafik Option
#################################################################
plotting <- function(nf, fac_names, n.hypotheses,
                     Descriptive.Factors, CI.method, ...){


  if(nf ==1){
    Faktor = fac_names[1]}
  
  if(nf > 1){
    print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
    Faktor <- scan("", what = "character")      
  }
  
  Fak.split <- strsplit(Faktor,":")[[1]]
  l.Fak.split <- length(Fak.split)
  
  
  if (!(Faktor %in% fac_names)) {
    stop("Please enter a valid factor name!")
  }
  
  
  for(i in 1:n.hypotheses){
    
    if(names(Descriptive.Factors)[i] == Faktor){
      posP <- which(names(Descriptive.Factors)[i] == Faktor)
      DatenPlot <- data.frame(Descriptive.Factors[[i]])
      
      
      switch(CI.method, Logit={
        upper = DatenPlot$U.Logit 
        lower = DatenPlot$L.Logit},
        Normal={
          upper =DatenPlot$U.Normal
          lower = DatenPlot$L.Normal})
      
      
      if (l.Fak.split==1){
        print(xyplot(pd ~ DatenPlot[,1], group=DatenPlot[,1],data = DatenPlot, 
                     type = 'p',xlab=paste(names(DatenPlot[1])),
                     col = 1, pch = 7, cex = 1.3, ylim = c(0, 1),
                     ylab="",upper = upper,
                     lower = lower,
                     panel = function(x, y, ...){
                     panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                       upper <- upper[subscripts]
                       lower <- lower[subscripts]
                       panel.arrows(x, lower, x, upper,code=4,lwd=4)   
                       panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
                       panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
                     } , ...)
                       panel.xyplot(x, y, ...)}))}
      
      if (l.Fak.split==2){
        print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p',
                     xlab=paste(names(DatenPlot[1])),
                     col = 1, pch = 7, cex = 1.3, ylim = c(0, 1),
                     upper = upper,
                     lower = lower,
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                         panel.arrows(x, lower, x, upper,code=4,lwd=4)   
                         panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
                         panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     }))
        
      }
      
      
      if (l.Fak.split==3){
        print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2]*DatenPlot[,3], 
                     group=DatenPlot[,1],data = DatenPlot, type = 'p',
                     xlab=paste(names(DatenPlot[1])),
                     col = 1, pch = 7, cex = 1.3, ylim = c(0, 1),
                     upper = upper,
                     lower = lower,
                     panel = function(x, y, ...){
                       panel.superpose(x, y,panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
                         upper <- upper[subscripts]
                         lower <- lower[subscripts]
                         panel.arrows(x, lower, x, upper,code=4,lwd=4)   
                         panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
                         panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
                       }, ...)
                       panel.xyplot(x, y, ...)
                     }))
        
      }
      
      if (l.Fak.split>=4){
        stop("4 and higher way interactions cannot be plotted!")
      }
    }}
}