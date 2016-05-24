station.symbol <-
function (cx, cy, direction=0, speed=NA, fill=rep(0,length(cx)), temp=NA, press=NA, dewpt=NA, circle=TRUE, cex=1)
{
 ### press is actually height in upper air ###
 ns = length(cx)
 if (length(cy) != ns)	stop("X AND Y COORDINATES SHOULD HAVE SAME LENGTH!")
 msg = "ALL VARIABLES SHOULD HAVE SAME LENGTH AS COORDINATES, OR BE MISSING!!!"
 if (ns > 1) {
   if (length(direction) == 1)
     if (!is.na(direction))     stop(msg)
   if (length(speed) == 1)
     if (!is.na(speed))     stop(msg)
   if (length(fill) == 1)
     if (!is.na(fill))       stop(msg)
   if (length(temp) == 1)
     if (!is.na(temp))     stop(msg)
   if (length(press) == 1)
     if (!is.na(press))   stop(msg)
   if (length(dewpt) == 1)
     if (!is.na(dewpt))   stop(msg)
   if (length(direction) > 1 & length(direction) != ns)	stop(msg)
   if (length(speed) > 1 & length(speed) != ns)	stop(msg)
   if (length(fill) > 1 & length(fill) != ns)	stop(msg)
   if (length(temp) > 1 & length(temp) != ns)	stop(msg)
   if (length(press) > 1 & length(press) != ns)	stop(msg)
   if (length(dewpt) > 1 & length(dewpt) != ns)	stop(msg)
 }
 ### mean/sd in order to identify outliers (> 3sd) ###
 if (ns > 3) {
   mspd = mean(as.numeric(speed), na.rm=TRUE)
   rspd = sd(as.numeric(speed), na.rm=TRUE)
   mt = mean(as.numeric(temp), na.rm=TRUE)
   rt = sd(as.numeric(temp), na.rm=TRUE)
   mps = mean(as.numeric(press), na.rm=TRUE)
   rps = sd(as.numeric(press), na.rm=TRUE)
   mdpt = mean(as.numeric(dewpt), na.rm=TRUE)
   rdpt = sd(as.numeric(dewpt), na.rm=TRUE)
 }

 tpar <- par()
 size <- tpar$csi
 scalex <- (tpar$usr[2] - tpar$usr[1])/(tpar$pin[1])
 scaley <- (tpar$usr[4] - tpar$usr[3])/(tpar$pin[2])
 scalex <- (cex * (scalex * size))/5
 scaley <- (cex * (scaley * size))/5
 for (i in 1:ns) {
   spdcolor = "green"
   tcolor = "blue"
   dptcolor = "blue"
   pscolor = "blue"

   x = cx[i]
   y = cy[i]
   if (is.na(x) | is.na(y))	next
   spd = speed[i]
   t   = temp[i]
   ps  = press[i]
   dpt = dewpt[i]

   if (circle) {
     ts <- seq(0, 2 * pi, , 200)
     RX <- sin(ts) * scalex
     X1 <- RX + x
     RY <- cos(ts) * scaley
     Y1 <- RY + y
     if (!is.na(spd)) {
       if (spd == 0)	{
         if (ns > 3) {
           if (spd < mspd-3*rspd)
             spdcolor = "red"
         }
         lines(RX * 2 + x, RY * 2 + y, col=spdcolor)
       }
     }
     if (fill[i] > 0) {
         lim <- c(51, 101, 151, 200)
         polygon(c(x, X1[1:lim[fill[i]]]), c(y, Y1[1:lim[fill[i]]]), 
             density = -1, col = "green")
     }
     lines(RX + x, RY + y, col=spdcolor)
     if (!is.na(t)) {
       if (ns > 3) {
         if (t > mt+3*rt | t < mt-3*rt)
           tcolor = "red"
       }
       temp.text <- paste(t, sep = "")
       temp.x <- x - 3.75 * scalex
       temp.y <- y + 1.75 * scaley
       text(temp.x, temp.y, labels = temp.text, col=tcolor)
     }
     if (!is.na(ps)) {
       if (ns > 3) {
         if (ps > mps+3*rps | ps < mps-3*rps)
           pscolor = "red"
       }
       press.text <- paste(ps, sep = "")
       press.x <- x + 4.75 * scalex
       press.y <- y + 1.25 * scaley
       text(press.x, press.y, labels = press.text, col=pscolor)
     }
     if (!is.na(dpt)) {
       if (ns > 3) {
         if (dpt > mdpt+3*rdpt | dpt < mdpt-3*rdpt)
           dptcolor = "red"
       }
       dewpt.text <- paste(dpt, sep = "")
       dewpt.x <- x - 3.75 * scalex
       dewpt.y <- y - 2.25 * scaley
       text(dewpt.x, dewpt.y, labels = dewpt.text, col=dptcolor)
     }
   } #end of circle

   if (!is.na(spd)) {
     xs <- if (spd > 0) {
       X1 <- 0
       X2 <- 0
       Y1 <- 0
       Y2 <- 5
       if (spd >= 5 & spd < 10) {
           X1 <- c(X1, 0)
           X2 <- c(X2, 1)
           Y1 <- c(Y1, 5)
           Y2 <- c(Y2, 5)
       }
       if (spd >= 10 & spd < 15) {
           X1 <- c(X1, 0)
           X2 <- c(X2, 2)
           Y1 <- c(Y1, 5)
           Y2 <- c(Y2, 5)
       }
       if (spd >= 15 & spd < 20) {
           X1 <- c(X1, 0, 0)
           X2 <- c(X2, 1, 2)
           Y1 <- c(Y1, 4, 5)
           Y2 <- c(Y2, 4, 5)
       }
       if (spd >= 20 & spd < 25) {
           X1 <- c(X1, 0, 0)
           X2 <- c(X2, 2, 2)
           Y1 <- c(Y1, 4, 5)
           Y2 <- c(Y2, 4, 5)
       }
       if (spd >= 25 & spd < 30) {
           X1 <- c(X1, 0, 0, 0)
           X2 <- c(X2, 1, 2, 2)
           Y1 <- c(Y1, 3, 4, 5)
           Y2 <- c(Y2, 3, 4, 5)
       }
       if (spd >= 30 & spd < 35) {
           X1 <- c(X1, 0, 0, 0)
           X2 <- c(X2, 2, 2, 2)
           Y1 <- c(Y1, 3, 4, 5)
           Y2 <- c(Y2, 3, 4, 5)
       }
       if (spd >= 35 & spd < 40) {
           X1 <- c(X1, 0, 0, 0, 0)
           X2 <- c(X2, 1, 2, 2, 2)
           Y1 <- c(Y1, 2, 3, 4, 5)
           Y2 <- c(Y2, 2, 3, 4, 5)
       }
       if (spd >= 40 & spd < 45) {
           X1 <- c(X1, 0, 0, 0, 0)
           X2 <- c(X2, 2, 2, 2, 2)
           Y1 <- c(Y1, 2, 3, 4, 5)
           Y2 <- c(Y2, 2, 3, 4, 5)
       }
       if (spd >= 45 & spd < 50) {
           X1 <- c(X1, 0, 0, 0, 0, 0)
           X2 <- c(X2, 1, 2, 2, 2, 2)
           Y1 <- c(Y1, 1, 2, 3, 4, 5)
           Y2 <- c(Y2, 1, 2, 3, 4, 5)
       }
       if (spd >= 50 & spd < 55) {
           X1 <- c(X1, 0, 0)
           X2 <- c(X2, 2, 2)
           Y1 <- c(Y1, 4, 5)
           Y2 <- c(Y2, 4.5, 4.5)
       }
       if (spd >= 55 & spd < 60) {
           X1 <- c(X1, 0, 0, 0)
           X2 <- c(X2, 1, 2, 2)
           Y1 <- c(Y1, 3, 4, 5)
           Y2 <- c(Y2, 3, 4.5, 4.5)
       }
       if (spd >= 60 & spd < 65) {
           X1 <- c(X1, 0, 0, 0)
           X2 <- c(X2, 2, 2, 2)
           Y1 <- c(Y1, 3, 4, 5)
           Y2 <- c(Y2, 3, 4.5, 4.5)
       }
       dir <- (direction[i]/360) * 2 * pi
       rot <- cbind(c(cos(dir), -sin(dir)), c(sin(dir), 
           cos(dir)))
       S1 <- rbind(X1, Y1)
       S2 <- rbind(X2, Y2)
       S1 <- rot %*% S1
       S2 <- rot %*% S2
       S1 <- S1 * c(scalex, scaley) + c(x, y)
       S2 <- S2 * c(scalex, scaley) + c(x, y)
     }
     if (spd > 0) {
       if (ns > 3) {
         if (spd > mspd+3*rspd | spd < mspd-3*rspd)
           spdcolor = "red"
       }
       segments(S1[1, ], S1[2, ], S2[1, ], S2[2, ], col=spdcolor, lwd=2)
     }
   } #end of (!is.na(spd))
 } #end of ns
 invisible()
}
