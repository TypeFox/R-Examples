
opar=par(no.readonly=TRUE)
dev.new()
par(mfrow=c(1,1))

DAT = randFRY(200)
   plot(DAT$x, DAT$y, asp=1, type='n', ask=FALSE)

for(i in 1:length(DAT$x))
  {
    dcirc = GEOmap::darc(rad = 5/2, ang1 = 0, ang2 = 360, x1 =DAT$x[i],  y1 = DAT$y[i] , n = 10)
    lines(dcirc, col='blue')
  }
points(DAT$x, DAT$y, pch=3)


##################

FF = dofry(DAT$x, DAT$y)
plotfry(FF, dis=30)
title("Fry Plot")


##################

dev.new(width=6, height=10)


RDAT = randFRY(400, LIM=c(0,0, 200, 200) , rlen=5   )
par(mfrow=c(3,2))
    shr = 0.0
simpleshear = matrix(c(1, shr, 0,  1), ncol=2)

Showfry(RDAT, simpleshear, 75)
title(main="no shear")

    shr = 1.2
simpleshear = matrix(c(1, shr, 0,  1), ncol=2)

Showfry(RDAT, simpleshear, 75)
title(main=paste("shear=",1.2))

epsilon1 = 0.4
H = matrix(c(1+epsilon1, 0, 0,  1/(1+epsilon1) ), ncol=2)

Showfry(RDAT, H, 75)
title(main=paste("epsilon=", epsilon1 ))


par(opar)

