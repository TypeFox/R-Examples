pyramid <- function(data, Laxis=NULL, Raxis=NULL, 
 AxisFM="g", AxisBM="", AxisBI=3, Cgap=0.3, Cstep=1, Csize=1, 
 Llab="Males", Rlab="Females", Clab="Ages", GL=TRUE, Cadj=-0.03, 
 Lcol="Cyan", Rcol="Pink", Ldens=-1, Rdens=-1, main="", ...) {
 # A function to draw population pyramid
 # rev 1.0: 5th January 2010
 # rev 1.1: 6th January 2010: Added "Cadj" option with faint modification.
 # rev 1.2: 11th March 2010: Added "Csize", "AxisFM", "AxisBM", and "AxisBI"
 #          options, as suggested by Dr. Philippe Guillet.
 # rev 1.3: 15th March 2012: URL changed and GunmaPop2005 dataframe was included.
 # (C) Minato Nakazawa <minato-nakazawa@umin.net>
 Left <- data[,1]
 Right <- data[,2]
 if (length(data)==2) { Center <- row.names(data) } else { Center <- data[,3] }
 if (is.null(Laxis)) { Laxis <- seq(0,max(c(Left,Right)),len=5) }
 if (is.null(Raxis)) { Raxis <- Laxis }
 # setting x-y axes
 BX <- c(-1-Cgap/2,1+Cgap/2)
 BY <- c(-0.05,1.1)
 plot(BX,BY,type="n",axes=FALSE,xlab="",ylab="",main=main,...)
 # scaling factors
 LL <- max(Laxis)
 LR <- min(Laxis)
 LS <- LL-LR
 LI <- length(Laxis)
 RL <- min(Raxis)
 RR <- max(Raxis)
 RS <- RR-RL
 RI <- length(Raxis)
 # ticks of axis
 segments(-(Laxis-LR)/LS-Cgap/2,-0.01,-(Laxis-LR)/LS-Cgap/2,0.01)
 segments((Raxis-RL)/RS+Cgap/2,-0.01,(Raxis-RL)/RS+Cgap/2,0.01)
 # vertical grid lines
 if (GL) {
  segments(-(Laxis-LR)/LS-Cgap/2,0,-(Laxis-LR)/LS-Cgap/2,1,lty=3,col="blue")
  segments((Raxis-RL)/RS+Cgap/2,0,(Raxis-RL)/RS+Cgap/2,1,lty=3,col="blue")
 }
 # axes
 lines(c(-1-Cgap/2,-Cgap/2),c(0,0),lty=1)
 lines(c(-Cgap/2,-Cgap/2),c(0,1),lty=1)
 lines(c(1+Cgap/2,Cgap/2),c(0,0),lty=1)
 lines(c(Cgap/2,Cgap/2),c(0,1),lty=1)
 # labels
 text(-0.5-Cgap/2,1,Llab,pos=3)
 text(0.5+Cgap/2,1,Rlab,pos=3)
 text(0,1,Clab,pos=3)
 Ci <- length(Center)
 for (i in 0:(Ci-1)) { 
  if ((i%%Cstep)==0) { text(0,i/Ci+Cadj,paste(Center[i+1]),pos=3,cex=Csize) }
 }
 text(-(Laxis-LR)/LS-Cgap/2,rep(0,LI),
  paste(formatC(Laxis,format=AxisFM,big.mark=AxisBM,big.interval=AxisBI)),pos=1)
 text((Raxis-RL)/RS+Cgap/2,rep(0,RI),
  paste(formatC(Raxis,format=AxisFM,big.mark=AxisBM,big.interval=AxisBI)),pos=1)
 # draw rectangles
 VB <- 0:(Ci-1)/Ci
 VT <- 1:Ci/Ci
 LeftP <- -(Left-LR)/LS-Cgap/2
 rect(LeftP,VB,rep(-Cgap/2,Ci),VT,col=Lcol,density=Ldens)
 RightP <- (Right-RL)/RS+Cgap/2
 rect(rep(Cgap/2,Ci),VB,RightP,VT,col=Rcol,density=Rdens)
}

pyramidf <- function(data, Laxis=NULL, Raxis=NULL,
                    frame=c(-1.15, 1.15, -0.05, 1.1),
                    AxisFM="g", AxisBM="", AxisBI=3, Cgap=0.3, Cstep=1, Csize=1, 
                    Llab="Males", Rlab="Females", Clab="Ages", GL=TRUE, Cadj=-0.03, 
                    Lcol="Cyan", Rcol="Pink", Ldens=-1, Rdens=-1, main="", ...) {
# frame version, added since rev. 1.4, 4th September 2014.
# (C) Minato Nakazawa <minato-nakazawa@umin.net>
  Left <- data[,1]
  Right <- data[,2]
  if (length(data)==2) { Center <- row.names(data) } else { Center <- data[,3] }
  if (is.null(Laxis)) { Laxis <- seq(0,max(c(Left,Right)),len=5) }
  if (is.null(Raxis)) { Raxis <- Laxis }
  # setting x-y axes
  BX <- c(-1-Cgap/2,1+Cgap/2)
  BY <- c(-0.05,1.1)
  XC <- function(XB) { (XB-BX[1])*(frame[2]-frame[1])/(2+Cgap)+frame[1] }
  YC <- function(YB) { (YB-BY[1])*(frame[4]-frame[3])/1.15+frame[3] }
  # scaling factors
  LL <- max(Laxis)
  LR <- min(Laxis)
  LS <- LL-LR
  LI <- length(Laxis)
  RL <- min(Raxis)
  RR <- max(Raxis)
  RS <- RR-RL
  RI <- length(Raxis)
  # ticks of axis
  segments(XC(-(Laxis-LR)/LS-Cgap/2),YC(-0.01),XC(-(Laxis-LR)/LS-Cgap/2),YC(0.01))
  segments(XC((Raxis-RL)/RS+Cgap/2),YC(-0.01),XC((Raxis-RL)/RS+Cgap/2),YC(0.01))
  # vertical grid lines
  if (GL) {
    segments(XC(-(Laxis-LR)/LS-Cgap/2),YC(0),XC(-(Laxis-LR)/LS-Cgap/2),YC(1),
             lty=3,col="blue")
    segments(XC((Raxis-RL)/RS+Cgap/2),YC(0),XC((Raxis-RL)/RS+Cgap/2),YC(1),
             lty=3,col="blue")
  }
  # axes
  lines(c(XC(-1-Cgap/2),XC(-Cgap/2)),c(YC(0),YC(0)),lty=1)
  lines(c(XC(-Cgap/2),XC(-Cgap/2)),c(YC(0),YC(1)),lty=1)
  lines(c(XC(1+Cgap/2),XC(Cgap/2)),c(YC(0),YC(0)),lty=1)
  lines(c(XC(Cgap/2),XC(Cgap/2)),c(YC(0),YC(1)),lty=1)
  # labels
  text(XC(-0.5-Cgap/2),YC(1),Llab,pos=3)
  text(XC(0.5+Cgap/2),YC(1),Rlab,pos=3)
  text(XC(0),YC(1),Clab,pos=3)
  Ci <- length(Center)
  for (i in 0:(Ci-1)) { 
    if ((i%%Cstep)==0) { text(XC(0),YC(i/Ci+Cadj),paste(Center[i+1]),pos=3,cex=Csize) }
  }
  text(XC(-(Laxis-LR)/LS-Cgap/2),YC(rep(0,LI)),
       paste(formatC(Laxis,format=AxisFM,big.mark=AxisBM,big.interval=AxisBI)),pos=1)
  text(XC((Raxis-RL)/RS+Cgap/2),YC(rep(0,RI)),
       paste(formatC(Raxis,format=AxisFM,big.mark=AxisBM,big.interval=AxisBI)),pos=1)
  # main text (above the frame)
  if (length(main)>0) { text(XC(0), YC(1.1), main, pos=3) }
  # draw rectangles
  VB <- 0:(Ci-1)/Ci
  VT <- 1:Ci/Ci
  LeftP <- -(Left-LR)/LS-Cgap/2
  rect(XC(LeftP),YC(VB),XC(rep(-Cgap/2,Ci)),YC(VT),col=Lcol,density=Ldens)
  RightP <- (Right-RL)/RS+Cgap/2
  rect(XC(rep(Cgap/2,Ci)),YC(VB),XC(RightP),YC(VT),col=Rcol,density=Rdens)
}

pyramids <- function(Left, Right, Center=NULL, ...) {
# Wrapper funuction for pyramid to use separate two vectors
 if (is.null(Center)) {
  dx <- data.frame(Left, Right, row.names=names(Left))
 } else { dx <- data.frame(Left, Right, Center) }
 pyramid(dx, ...)
}

# The census 2005 for Gunma prefecture data obtained from
# http://www.e-stat.go.jp/SG1/estat/List.do?bid=000001005042&amp;cycode=0
# as "a00411.xls"
GunmaPop2005 <- data.frame(
 Males=c(8872, 9144, 9528, 9812, 9817, 10049, 10234, 10047, 10222, 
 10187, 10319, 10420, 10135, 10473, 10219, 10492, 10943, 11243, 10579, 
 9653, 9941, 10252, 10396, 10649, 11343, 12031, 12420, 13015, 13354, 
 14132, 14624, 15873, 16013, 15809, 15338, 14876, 14568, 14536, 14580, 
 10713, 13844, 12911, 12303, 12158, 12086, 12117, 12620, 12331, 12376, 
 13138, 14144, 13553, 14353, 15194, 16064, 17286, 18561, 18347, 18230, 
 12001, 11938, 14363, 13684, 14122, 13475, 12241, 10279, 10866, 10939, 
 10853, 10093, 9996, 9575, 9334, 9147, 8643, 8249, 7760, 7353, 6812, 
 6170, 5041, 4099, 3206, 2711, 2778, 2113, 1895, 1630, 1431, 1146, 966, 
 675, 559, 384, 263, 212, 133, 84, 38, 24, 15, 11, 6, 1, 2, 1, 0, 0),
Females=c(8323, 8750, 8964, 9359, 9559, 9605, 9511, 9800, 9790, 9848, 
 9939, 9755, 9696, 9884, 9734, 9767, 10293, 10616, 10040, 9383, 9731, 
 9748, 10211, 10204, 10500, 11086, 11758, 12248, 12548, 13524, 13907, 
 14930, 15129, 15057, 14344, 13982, 13942, 13587, 13972, 10359, 13212, 
 12435, 11934, 11673, 11668, 11583, 12100, 11867, 11917, 12746, 13364, 
 13170, 13968, 15318, 16251, 17125, 18253, 18042, 17927, 11981, 11773, 
 14450, 14124, 14438, 13502, 12960, 10729, 11710, 11697, 11884, 11413, 
 11442, 11087, 11035, 11209, 10646, 10482, 9784, 9777, 9491, 8891, 8188, 
 7636, 7034, 6123, 6103, 4577, 4415, 3861, 3426, 3035, 2571, 2064, 1683, 
 1209, 878, 739, 536, 333, 193, 134, 86, 56, 28, 13, 8, 5, 3, 1),
Ages=0:108)

# The census  2010 for Gunma prefecture data obtained from
# http://www.e-stat.go.jp/SG1/estat/GL02020101.do?
# method=csvDownload&fileId=000004983709&releaseCount=1
# (connect above 2 lines as URL) as 00310.csv
GunmaPop2010 <- data.frame(
 Males = c(8135, 8555, 8869, 8607, 8688, 9001, 9354, 9583, 9841, 9749, 
 10073, 10118, 10057, 10211, 10142, 10433, 10520, 10206, 9729, 8505, 
 8675, 8908, 9160, 9569, 9704, 10321, 10638, 10815, 10991, 11463, 11991, 
 12360, 13000, 13299, 14177, 14597, 15768, 16004, 15769, 15140, 14694, 
 14481, 14420, 14444, 10596, 13612, 12718, 12129, 12018, 11875, 11937, 
 12370, 12056, 12072, 12788, 13834, 13087, 13946, 14890, 15589, 16841, 
 17968, 17716, 17610, 11614, 11408, 13654, 12899, 13350, 12546, 11414, 
 9526, 9927, 9888, 9761, 8912, 8657, 8137, 7848, 7528, 6966, 6423, 
 5832, 5414, 4841, 4139, 3285, 2431, 1850, 1443, 1445, 991, 810, 653, 
 511, 367, 294, 175, 101, 57, 49, 35, 14, 7, 3, 1, 1, 0, 0),
 Females = c(7612, 7952, 8176, 8397, 8438, 8458, 8851, 8996, 9377, 9518, 
 9686, 9458, 9748, 9804, 9771, 9961, 9769, 9685, 9182, 8354, 8271, 8673, 
 8853, 8938, 9039, 9560, 9937, 10360, 10250, 10525, 11068, 11745, 12109, 
 12477, 13471, 13945, 15021, 14996, 15080, 14352, 13943, 13878, 13507, 
 13877, 10290, 13046, 12352, 11815, 11599, 11532, 11402, 11969, 11779, 
 11758, 12624, 13173, 13065, 13871, 15122, 16013, 16974, 18054, 17850, 
 17627, 11820, 11598, 14099, 13765, 14021, 13154, 12569, 10288, 11267, 
 11209, 11325, 10843, 10727, 10284, 10252, 10177, 9673, 9339, 8711, 8424, 
 8122, 7342, 6478, 5841, 5183, 4330, 4017, 2847, 2599, 2163, 1715, 1304, 
 1048, 743, 539, 333, 232, 151, 101, 57, 29, 14, 8, 3, 3),
Ages=0:108)
