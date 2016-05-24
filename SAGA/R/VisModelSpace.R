VisModelSpace <- function(data, cex.u = 3, cex.mtext=1, cex.leg=.8){
  models <- data[[1]]
  data <- data[[3]]
  colors <- heat.colors(1101)
  colors <- colors[1001:1]
  waic <- (exp(-.5 * data) / (sum(exp(-.5 * data))))
  foobar <- waic/max(waic)
  mod.heat <- round(1000 * foobar) + 1
  foo <- round(sqrt(length(data)) + 1)
  yvals <- rep(1:foo, each=foo)[1:length(data)]
  xvals <- rep(1:foo, times=foo)[1:length(data)]
  par(mar=c(2, 2, 2, 7))
  
  ## add the increasing level of model complexity
  
  plot(x=xvals, y=yvals, col=colors[mod.heat], pch=15,
       xaxt="n", yaxt="n", cex = cex.u, 
       )
  x="Increasing Model Complexity"
  y=""
  foobar <- substitute(x %->% y)
  mtext(text=foobar, line=0.45, side=1, cex=cex.mtext)
  mtext(text=foobar, line=0, side=2, cex=cex.mtext)
  abbi <- 1:length(data)
  #for(i in 1:length(data)){
  #  text(x=xvals[i], y=yvals[i], labels = abbi[i], cex = cex.u/4, col= "gray15")
  #}
  locs <- par("usr")  
  color.legend(locs[2] + (locs[2]*.01),  #xl
               locs[4]- ((locs[4]-locs[3]) * 0.5),  #yb
               locs[2] + (locs[2]*.05),    #xr
               locs[4],        #yt
               legend = c("0.00", round(max(waic), digits=3)),
               rect.col = colors[1:length(colors)],
               cex = cex.leg,
               align = "rb", gradient = "y")
  par(xpd = TRUE)
  text(x=(((locs[2] + 1.5) + (locs[2] + (locs[2]*.05)))/2),
       y=(locs[4]+.2), labels="Akaike Weights", pos=3, cex=cex.leg)
    
}
