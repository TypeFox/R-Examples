dplot.cgb2 <- function(x, shape1, scale, shape2, shape3, pl0, pl, w=rep(1,length(x)) ,decomp="r",
                       xmax = max(x)*(2/3), choicecol=1:3, kernel="epanechnikov", adjust=1, title=NULL, ylim=NULL){
      para <- paste(" a=",round(shape1,2),", b=",round(scale),", p=",round(shape2,2),", q=",round(shape3,2))
      sub0=paste("pl0 = (",round(pl0[1],3)) 
      sub=paste("pl = (",round(pl[1],3)) 
      pl1 <- length(pl0)-1
      if (pl1 >= 2){   
          for (i in 2:pl1) {
              sub0 <- paste(sub0,",", round(pl0[i],3))
              sub <- paste(sub,",", round(pl[i],3))
          }
         }
    sub0 <- paste(sub0,",",round(pl0[pl1+1],3),")")
    sub <- paste(sub,",",round(pl[pl1+1],3),")")
    fcgb2 <- function(x) dcgb2(x, shape1, scale, shape2, shape3, pl0, pl, decomp=decomp)
    fgb2 <- function(x) dgb2(x, shape1, scale, shape2, shape3)
#    maxx <- max(x)*2/3   # change 28.04.2014: put as argument
    if (is.null(ylim)) curve(fcgb2, col=choicecol[2], lwd=2, lty=1, from=0, to=xmax, ylab="Density")
    else curve(fcgb2, col=choicecol[2], lwd=2, lty=1, from=0, to=xmax, ylab="Density", ylim=ylim)
      curve(fgb2,col=choicecol[1], lwd=2, lty=2, add=TRUE)
      if (is.null(title)) title <- "Comparison of densities"
   title(title, sub = paste(para,"; ",sub0,";",sub),
         cex.sub = 0.75, font.sub = 3)
      
     
     ## empirical counterparts
     wk <- w/sum(w)
     densk <- density(x, weights=wk, kernel= kernel, from=0, adjust=adjust)              
     lines(densk, col=choicecol[3], lwd=2, lty=3)
#     print("Please, place the cursor for the legend",quote = FALSE)                  
     legend("topright",c("GB2 ","compound GB2 ","Kernel estimate"),                   
       lwd=2,col=choicecol[c(1,2,3)], lty=c(2,1,3))
}