plotDB<-function(DB)
  {
   stacomp = paste(sep=".", DB$sta, DB$comp)



usta = unique(DB$sta)

     pald = c("black", "darkmagenta", "forestgreen", "blueviolet", 
        "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3", 
        "magenta1", "lightsalmon3", "darkcyan", "darkslateblue", 
        "chocolate4", "goldenrod4", "mediumseagreen")

   
   usc = unique(stacomp)


   y = range( 1:length(usc))
   x = range( c(DB$t1, DB$t2))

  SPT1 = split(DB$t1, stacomp)
   SPT2 = split(DB$t2, stacomp)


   uu  = names(SPT1)
   gs = strsplit(uu, split="\\.")

     stass = unlist(lapply(gs, getmem, 1))

   im = match(  stass, usta)

  ###########     1+((im-1) %% length(pald))
   icols = pald[  1+((im-1) %% length(pald)) ]
   

   
    days = pretty(x)

    days = days[days!=0]

   day.yrs =  EPOCHyear(days, origyr=attr(DB, "origyr")  )


  
   plot(x,y, type="n", axes=FALSE, xlab="", ylab="")
   abline(h=1:length(SPT1), lty=2, col=grey(.9))
   box()

   u = par("usr")

   axis(4, pos=u[1]  , at=1:length(SPT1), labels=names(SPT1), las=1, cex=.8, padj=-0.5)

   

  ##  axis(1, at=days, labels=paste(sep=" ", day.yrs$yr, day.yrs$jd), , las=2)
axis(3, at=days, labels=paste(sep="/", day.yrs$yr, day.yrs$jd))
axis(1, at=days, labels=paste(sep="/", day.yrs$yr, day.yrs$jd))

   
   
   abline(v=days, col=grey(.85), lty=2)

   
   for(i in 1:length(SPT1))
     {
       why = rep(i, length(SPT1[[i]]))
       segments(SPT1[[i]], why, SPT2[[i]], why, col=icols[i], lwd=2.5)


     }

   

  }
