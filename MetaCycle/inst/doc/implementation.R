## ------------------------------------------------------------------------
library(MetaCycle)
head(cycMouseLiverRNA[,1:5])

## ------------------------------------------------------------------------
set.seed(100)
row_index <- sample(1:nrow(cycHumanBloodDesign), 4)
cycHumanBloodDesign[row_index,]

## ------------------------------------------------------------------------
sample_id <- cycHumanBloodDesign[row_index,1]
head(cycHumanBloodData[,c("ID_REF", sample_id)])

## ------------------------------------------------------------------------
group_index <- which(cycHumanBloodDesign[, "group"] == "SleepExtension")
cycHumanBloodDesignSE <- cycHumanBloodDesign[group_index,]
sample_index <- which(cycHumanBloodDesignSE[, "subject"] == "AF0004")
sample_AF0004 <- cycHumanBloodDesignSE[sample_index, "sample_library"]
cycHumanBloodDataSE_AF0004 <- cycHumanBloodData[, c("ID_REF", sample_AF0004)]
head(cycHumanBloodDataSE_AF0004)

## ---- echo=FALSE, warning=FALSE, fig.width=6.65, fig.height=5------------
library(shape)
library(diagram)
par(mar = c(1, 0.5, 1, 0.5))
openplotmat()
#number of elements in each row
num_element <- c(1, 1, 2, 3, 4, 5)
#get position information of each element in the flow chart
elpos <- coordinates (num_element, mx = 0)
#adjust x-position of some elements
elposM <- elpos
elposM[1:2, 1] <- elposM[1:2, 1] - 0.21
elposM[4,1] <- elposM[4,1] - 0.269
elposM[7,1] <- elposM[7,1] - 0.186
elposM[11,1] <- elposM[11,1] - 0.076
#give information of strat and end point of each arrow
fromto <- c( c(1, 2), c(2,12), c(2,4), c(4,13), c(4, 7),
             c(7, 14), c(7,11), c(11, 15), c(11, 16) )
fromtoM <- matrix(ncol = 2, byrow = TRUE, data = fromto)
rownum <- nrow(fromtoM)
#draw arrow and store arrow position informaion in 'arrposM'
arrposM <- matrix(ncol = 2, nrow = rownum)
for (i in 1:rownum)
{
    arrposM[i, ] <- bentarrow (from = elposM[fromtoM[i, 1], ], to = elposM[fromtoM[i, 2], ], 
                               lcol = "blue", lwd = 1, arr.pos = 0.56, arr.length = 0.3, arr.lwd = 0.8)
}
#draw elements of flow chart
textparallel(mid = elposM[1,], radx = 0.089, rady=0.06, lab = "Time-series\ndata", lcol = "blue", 
             lwd=2, shadow.size = 0, cex = 0.72, font=2, theta=80)
textdiamond(mid =  elposM[2,], radx = 0.12, rady=0.066, lab = "With\nnon-integer\nintervals?", 
            lcol = "blue", lwd=2, shadow.size = 0, cex = 0.72, font=2)
diamond_index <- c(4, 7, 11)
diamond_lab <- c("Uneven\nsampling?", "Missing\nvalue?", "With\nreplicates?")
for (i in 1:length(diamond_index))
{
  textdiamond(mid = elposM[diamond_index[i],], radx = 0.08, rady=0.066, lab = diamond_lab[i], 
              lcol = "blue", lwd=2, shadow.size = 0, cex = 0.72, font=2)  
}

round_index <- 12:16
round_lab <- c("LS", "LS", "JTK&LS", "JTK&LS", "ARS&\nJTK&LS")
for (j in 1:length(round_index))
{
    textround(mid = elposM[round_index[j],], radx=0.056, rady=0.05, lab = round_lab[j],
              lcol = "blue", lwd=2, shadow.size = 0, cex =0.76, font=2, rx = 0.02)
}
#add 'Y' and 'N' on the flow chart
midposM <- elposM[c(2, 4, 7, 11),]
xpos <- midposM[,1]
ypos <- midposM[,2]
YposM <- cbind(c( (xpos[1] - 0.145), (xpos[2:3] - 0.1), (xpos[4] - 0.086) ), ypos + 0.02)
NposM <- cbind(c( (xpos[1] + 0.145), (xpos[2:3] + 0.1), (xpos[4] + 0.086) ), ypos + 0.02)
text(YposM[,1], YposM[,2], labels="Y", cex=0.8, font=2)
text(NposM[,1], NposM[,2], labels="N", cex=0.8, font=2)

## ---- warning=FALSE------------------------------------------------------
# given three phases
pha <- c(0.9, 0.6, 23.6)
# their corresponding periods
per <- c(23.5, 24, 24.5)
# mean period length
per_mean <- mean(per)
# covert to polar coordinate
polar <- 2*pi*pha/per
# get averaged ploar coordinate
polar_mean <- atan2(mean(sin(polar)), mean(cos(polar)))
# get averaged phase value
pha_mean <- per_mean*polar_mean/(2*pi)
pha_mean

## ---- echo=FALSE, warning=FALSE, fig.width=6.65, fig.height=5------------
getAMP <- function(expr, per, pha, tim=18:65)
{ 
    trendt <- tim - mean(tim[!is.na(tim) & !is.nan(tim)])
    cost <- cos(2*pi/per*(tim - pha))
    fit <- lm(expr~trendt + cost)
    fitcoef <- fit$coefficients
    basev <- fitcoef[1]
    trendv <- fitcoef[2]
    ampv <- fitcoef[3]
    fitexp <- basev + trendv*trendt + ampv*cost
    outL <- list("base"=basev, "trend"=trendv, "amp"=ampv, "fit"=fitexp)
    return(outL)
}

cirD <- cycVignettesAMP
ampL <- getAMP(expr=as.numeric(cirD[1,24:71]), per=cirD[1, "meta2d_period"], pha=cirD[1, "meta2d_phase"])

lay<-layout(cbind(1, 2), widths=c( lcm(cm(4.5)), lcm(cm(1.5)) ), heights=lcm(cm(4.5)) )
par(mai=c(0.65,0.6,0.4,0.05),mgp=c(2,0.5,0),tck=-0.01)
xrange <- c(18, 65)
yrange <- c(200, 2350)
plot(18:65, cirD[1,24:71], type="b", xlim=xrange, ylim=yrange, xlab="Circadian time(CT)", ylab="Expression value",  main=cirD[1,1], cex.main=1.2)
par(new=T)
plot(18:65, ampL[[4]], type="b", xlim=xrange, ylim=yrange, col="red", xlab="", ylab="", main="")
abline(h=ampL[[1]], lty=3, col="purple", lwd=1.5)
lines(18:65, 500+ampL[[2]]*(18:65-mean(18:65)), lty=4, col="orange", lwd=1.5)
legend("topleft", legend=c("Raw value", "OLS fitted value"), col=c("black", "red"), pch=1, bty="n")
legend("topright", legend=c("Baseline", "Trend"), col=c("purple", "orange"), lty=c(3, 4), lwd=1.5, bty="n" )

par(mai=c(0.5,0.05,0.4,0.1),mgp=c(2,0.3,0),tck=-0.01);
plot(x=NULL,y=NULL,xlim=c(0,10),ylim=c(0,10),type="n", xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="")
text(rep(1,3), c(8, 5, 2), c("Base = ", "Trend = ", "AMP = "), adj=0)
text_value <- unlist(ampL)
text(rep(6,6), c(8, 5, 2), round(text_value[1:3], 1), adj=0)

## ---- echo=FALSE, warning=FALSE, fig.width=6.65, fig.height=5------------
cirD <- cycVignettesAMP
cirM <- as.matrix(cirD[2:3, 24:71])
expmax <- apply(cirM, 1, max)
cirM <- cirM/expmax

lay<-layout(cbind(1, 2), widths=c( lcm(cm(4.5)), lcm(cm(1.5)) ), heights=lcm(cm(4.5)) )
par(mai=c(0.65,0.6,0.4,0.05),mgp=c(2,0.5,0),tck=-0.01)
xrange <- c(18, 65)
yrange <- c(0, 1)
colname <- c("red", "blue")
grey_trans <- rgb(191/255,191/255,191/255,0.65);

par(mai=c(0.65,0.6,0.2,0.05),mgp=c(2,0.5,0),tck=-0.01)
plot(NULL,NULL,xlim=xrange,ylim=yrange,xaxt="n",yaxt="n",xlab="Circadian time(CT)",ylab="Exp/Max", main="");
rect_xL <- c(18, 36, 60)
rect_yL <- c(24, 48, 65)
rect(rect_xL, rep(-0.1, 3), rect_yL, rep(1.1,3), col=grey_trans, border=NA, bty="n")

for (i in 1:2)
{
  loessD <- data.frame(expd=as.numeric(cirM[i,]),tp=18:65);
  exploess <- loess(expd~tp, loessD, span = 0.2);
  expsmooth <- predict(exploess, data.frame(tp=18:65));
  lines(18:65,expsmooth,lwd=1.2,col=colname[i]);
}

xpos <- c(seq(18,60,by=6), 65)
axis(side=1,at=xpos,labels=xpos,mgp=c(0,0.2,0),tck=-0.01,cex.axis=0.8)
ypos <- seq(0, 1, by=0.2)
axis(side=2,at=ypos,labels=ypos,mgp=c(0,0.2,0),tck=-0.01,cex.axis=0.8)

par(mai=c(0.5,0.05,0.2,0.1),mgp=c(2,0.3,0),tck=-0.01)
plot(x=NULL,y=NULL,xlim=c(0,10),ylim=c(0,10),type="n", xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="")
lines(c(0.2,2.3), c(9.5,9.5), col="blue", lwd=1.5)
text(c(5,0,0), c(9.5, 8, 6.5), c("Ugt2b34", "AMP = ", "rAMP = "), col="blue", adj=0)
text(c(5,5), c(8, 6.5), round(as.numeric(cirD[3, 22:23]), 2), col="blue", adj=0)

lines(c(0.2,2.3), c(4.5,4.5), col="red", lwd=1.5)
text(c(5,0,0), c(4.5, 3, 1.5), c("Arntl", "AMP = ", "rAMP = "), col="red", adj=0)
text(c(5,5), c(3,1.5), round(as.numeric(cirD[2, 22:23]), 2), col="red", adj=0)

