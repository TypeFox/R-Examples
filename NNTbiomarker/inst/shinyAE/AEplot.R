#### Adverse event NNT bars ####
#### For CMFT, from table 5 of fisher1997.
#### For Supplement 2, setting the plot box small (342x 272) seems to work.

AEplot = function(RSinput = 30, makeTitle=FALSE){

  boxcolors = colorRampPalette(c("lightgrey", "red"))(6)

  if(is.null(RSinput))
    RSinput = 48
  par.save = par(mai=c(0,0,0.5,0.5))
  aeProb = c(2.9,15,57,20,5,0.1)
  boxwidths = c(1, (nnt[RSinput] - 1) * aeProb / 100)
  opts = options(warn=-1)
  symbols(x=rep(0, 7), y=1:7, inches=F,
          xlim=c(-ceiling(max(boxwidths)), ceiling(max(boxwidths))) * 0.75,
          rectangles = cbind(boxwidths, 1), bg = c("green", boxcolors),
          axes=F,
          xlab="", ylab="")
  par(par.save)
  options(opts)
  "%except%" <-  function (vector, condition) vector[match(vector, condition, 0) == 0]
  verticalsX = lapply(boxwidths[-1], function(bw)
    if(bw <= 1)  numeric(0)  else  -floor(bw/2):floor(bw/2)
  )
  verticalsY = rep(1:6, times=sapply(verticalsX, length))
  segments(x0= unlist(verticalsX),
           y0 = verticalsY + 1/2, y1 = verticalsY + 3/2
  )
  graphics::text(x = boxwidths/2, y=1:7,
       c("benefitted", "no adverse event", "mild AE", "moderate", "severe AE", "life-threatening AE",
         "fatal toxicity"),
       pos=4 , xpd=NA)
  graphics::text(x = - boxwidths/2, y=1:7, round(boxwidths, 1),
       pos=2 , xpd=NA)
  if(makeTitle)
    title(paste0("Outcomes for ",
                 round(nnt[RSinput]), " patients, \n",
                 "if all are treated, and 1 benefits\n",
                 "RS = ", RSinput, "  NNT = ", round(nnt[RSinput])))
  print(paste0("RS = ", RSinput, "  NNT = ", round(nnt[RSinput])))
  print(boxwidths)
}
#
#for(RSinput in c(OncotypeRScutoffs, TailorXRScutoffs))
#for(RSinput in c(17,30,10,25))
#  AEplot(RSinput)
