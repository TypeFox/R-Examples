bv.boxplot<-function(X, Y, robust=TRUE, D = 7, xlab="X", ylab="Y", pch = 21, 
pch.out = NULL, bg = "gray", bg.out = NULL, hinge.col = 1, fence.col = 1, 
hinge.lty = 2, fence.lty = 3, xlim = NULL, ylim = NULL, names = 1:length(X), 
ID.out = FALSE, cex.ID.out = 0.7, uni.CI = FALSE, uni.conf = 0.95, 
uni.CI.col = 1, uni.CI.lty = 1, uni.CI.lwd = 2, show.points = TRUE, ...){
if(length(X) != length(Y)) stop("X and Y must have equal length")
biweight<-function(a,const1=9,const2=36,err=0.0001){ ###From Everitt (2004)
              x<-a[,1]
              y<-a[,2]
              n<-length(x)
              mx<-median(x)
              my<-median(y)
              madx<-median(abs(x-mx))
              mady<-median(abs(y-my))
              if(madx != 0) { ux<-(x-mx)/(const1*madx)
                              ux1<-ux[abs(ux)<1]
                              tx<-mx+(sum((x[abs(ux)<1]-mx)*(1-ux1*ux1)^2)/
                                     sum((1-ux1^2)^2))
                              sx<- sqrt(n)*sqrt(sum((x[abs(ux)<1]-mx)^2*
                                   (1-ux1*ux1)^4))/abs(sum((1-ux1*ux1)*
                                    (1-5*ux1*ux1)))
                             }
                  else { tx<-mx
                         sx<-sum(abs(x-mx))/n
                       }
              if(mady != 0) { uy<-(y-my)/(const1*mady)
                              uy1<-uy[abs(uy)<1]
                              ty<-my+(sum((y[abs(uy)<1]-my)*(1-uy1*uy1)^2)/
                                     sum((1-uy1^2)^2))
                              sy<- sqrt(n)*sqrt(sum((y[abs(uy)<1]-my)^2*
                                   (1-uy1*uy1)^4))/abs(sum((1-uy1*uy1)*
                                    (1-5*uy1*uy1)))
                             }
                  else { ty<-my
                         sy<-sum(abs(y-my))/n
                       }
               z1<-(y-ty)/sy+(x-tx)/sx
               z2<-(y-ty)/sy-(x-tx)/sx
              mz1<-median(z1)
              mz2<-median(z2)
              madz1<-median(abs(z1-mz1))
              madz2<-median(abs(z2-mz2))
              if(madz1 != 0) { uz1<-(z1-mz1)/(const1*madz1)
                              uz11<-uz1[abs(uz1)<1]
                           tz1<-mz1+(sum((z1[abs(uz1)<1]-mz1)*(1-uz11*uz11)^2)/
                                     sum((1-uz11^2)^2))
                              sz1<- sqrt(n)*sqrt(sum((z1[abs(uz1)<1]-mz1)^2*
                                   (1-uz11*uz11)^4))/abs(sum((1-uz11*uz11)*
                                    (1-5*uz11*uz11)))
                             }
                  else { tz1<-mz1
                         sz1<-sum(abs(z1-mz1))/n
                       }
              if(mady != 0) { uz2<-(z2-mz2)/(const1*madz2)
                              uz21<-uz2[abs(uz2)<1]
                          tz2<-mz2+(sum((z2[abs(uz2)<1]-mz2)*(1-uz21*uz21)^2)/
                                     sum((1-uz21^2)^2))
                              sz2<- sqrt(n)*sqrt(sum((z2[abs(uz2)<1]-mz2)^2*
                                   (1-uz21*uz21)^4))/abs(sum((1-uz21*uz21)*
                                    (1-5*uz21*uz21)))
                             }
                  else { tz2<-mz2
                         sz2<-sum(abs(z2-mz2))/n
                       }
              esq<-((z1-tz1)/sz1)^2+((z2-tz2)/sz2)^2
              w<-numeric(length=n)
              c2<-const2
              for(i in 1:10) {
              w[esq<const2]<-(1-esq[esq<const2]/const2)^2
              w[esq>=const2]<-0
              l<-length(w[w==0])
              if(l<0.5*n) break
                  else const2<-2*const2
                            }
               tx<-sum(w*x)/sum(w)
               sx<-sqrt(sum(w*(x-tx)^2)/sum(w))
               ty<-sum(w*y)/sum(w)
               sy<-sqrt(sum(w*(y-ty)^2)/sum(w))
               r<-sum(w*(x-tx)*(y-ty))/(sx*sy*sum(w))
               const2<-c2
               wold<-w
            for(i in 1:100) {
                     z1<-((y-ty)/sy+(x-tx)/sx)/sqrt(2*(1+r))
                     z2<-((y-ty)/sy-(x-tx)/sx)/sqrt(2*(1+r))
                     esq<-z1*z1+z2*z2
                     for(j in 1:10) {
                                    w[esq<const2]<-(1-esq[esq<const2]/const2)^2
                                    w[esq>=const2]<-0
                                    l<-length(w[w==0])
                                    if(l<0.5*n) break
                                         else const2<-2*const2
                                     }
               tx<-sum(w*x)/sum(w)
               sx<-sqrt(sum(w*(x-tx)^2)/sum(w))
               ty<-sum(w*y)/sum(w)
               sy<-sqrt(sum(w*(y-ty)^2)/sum(w))
               r<-sum(w*(x-tx)*(y-ty))/(sx*sy*sum(w))
               term<-sum((w-wold)^2)/(sum(w)/n)^2
               if(term<-err) break
                    else {wold<-w
                          const2<-c2
                         }
                             }
            param<-c(tx,ty,sx,sy,r)
            param
}

if(robust == TRUE){
bw <- biweight(cbind(X,Y))
R <- bw[5]
S.X <- bw[3]
S.Y <- bw[4]
X.loc <- bw[1]
Y.loc <- bw[2]
}
if(robust == FALSE){
R <- cor(X, Y)
S.X <- sd(X)
S.Y <- sd(Y)
X.loc <- mean(X)
Y.loc <- mean(Y)
}
X.stan <-(X - X.loc)/S.X
Y.stan <-(Y - Y.loc)/S.Y

E.i <- sqrt((X.stan^2 + Y.stan^2 - 2 * R * X.stan * Y.stan)/(1 - R^2))
E.m <- median(E.i)
E.max <- max(E.i[E.i^2 < D * E.m^2])

#Fence
R1f <- E.max * sqrt((1 + R)/2)
R2f <- E.max * sqrt((1 - R)/2)
theta <- seq(2, 360, by = 2)
theta.rad <-(theta * pi)/180
Theta1 <- R1f * cos(theta.rad)
Theta2 <- R2f * sin(theta.rad)
X.pf <- X.loc + (Theta1 + Theta2) * S.X
Y.pf <- Y.loc + (Theta1 - Theta2) * S.Y

#Hinge
R1 <- E.m * sqrt((1 + R)/2)
R2 <- E.m * sqrt((1 - R)/2)
Theta1 <- R1 * cos(theta.rad)
Theta2 <- R2 * sin(theta.rad)
X.p <- X.loc + (Theta1 + Theta2) * S.X
Y.p <- Y.loc + (Theta1 - Theta2) * S.Y


fmmx <- c(min(X.pf), max(X.pf))
fmmy <- c(min(Y.pf), max(Y.pf))

b1 <- (R * S.Y)/S.X
a1 <- Y.loc - b1 * X.loc
Y1 <- a1 + b1 * min(fmmx)
Y2 <- a1 + b1 * max(fmmx)
b2 <- (R * S.X)/S.Y
a2 <- X.loc - b2 * Y.loc
X1 <- a2 + b2 * min(fmmy)
X2 <- a2 + b2 * max(fmmy)
maxx <- max(c(max(X), max(fmmx), X1, X2))
minx <- min(c(min(X), min(fmmx), X1, X2))
maxy <- max(c(max(Y), max(fmmy), Y1, Y2))
miny <- min(c(min(Y), min(fmmy), Y1, Y2))
if(is.null(xlim)) xlim <- c(minx, maxx) 
if(is.null(ylim)) ylim <- c(miny, maxy)
plot(X, Y,  xlab = xlab, ylab = ylab, type = "n", xlim = xlim, ylim = ylim)
lines(X.pf, Y.pf, lty = fence.lty, col = fence.col)
lines(X.p, Y.p, lty = hinge.lty, col = hinge.col)
segments(min(fmmx), Y1, max(fmmx), Y2, lty = 1)
segments(X1, min(fmmy), X2, max(fmmy), lty = 1)

#------- Outlier ID ---------------#
bool <- ifelse(E.i > E.max, "out", "in")

if(is.null(pch.out)) pch.out <- pch
if(is.null(bg.out)) bg.out <- bg
pch <- ifelse(bool == "in", pch, pch.out)
bg <- ifelse(bool == "in", bg, bg.out)

if(show.points == TRUE) points(X, Y, pch = pch, bg = bg, ...)

outl <- data.frame(X = X, Y = Y)[bool == "out",] 
if(ID.out == TRUE){pos <- 2; if(nrow(outl) > 2) pos <-  thigmophobe(outl$X,outl$Y)
text(outl$X, outl$Y, row.names(outl), pos = pos, cex = cex.ID.out)}

if(uni.CI == TRUE){
cX <- ci.median(X, conf = uni.conf)$ci
cY <- ci.median(Y, conf = uni.conf)$ci
yjust <- abs(par("usr")[1]-par("usr")[2])*.02
xjust <- abs(par("usr")[3]-par("usr")[4])*.02
segments(par("usr")[1]+yjust,cY[2], par("usr")[1]+yjust,cY[3],lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
segments(cX[2],par("usr")[3]+xjust,cX[3], par("usr")[3]+xjust,lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
segments(par("usr")[1]+(.5*yjust),cY[2],par("usr")[1]+(yjust+.5*yjust),cY[2],lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
segments(par("usr")[1]+(.5*yjust),cY[3],par("usr")[1]+(yjust+.5*yjust),cY[3],lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
segments(cX[2],par("usr")[3]+(.5*xjust), cX[2],par("usr")[3]+(xjust+.5*xjust), lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
segments(cX[3],par("usr")[3]+(.5*xjust), cX[3],par("usr")[3]+(xjust+.5*xjust), lwd=uni.CI.lwd, lty = uni.CI.lty, col = uni.CI.col)
}
invisible(list(centroid = c(X.loc, Y.loc), scale = c(S.X, S.Y), correlation = R, E.median = E.m, E.max = E.max, outliers = outl))
}

 