# bass2.R -- 2012030630 DIM B(3), X(3), O(3,3), Y(80,5),Z(13)
rm(list=ls())
bdat<-read.csv(file="bass2.csv")
bd<-subset(bdat, ymale==1)
#3136 LET Y(M,1)=Z(7): LET Y(M,2)=Z(3): LET Y(M,3)=Z(4): REM S,W, t
#3137 LET Y(M,4)=Z(6): LET Y(M,5)=Z(8): REM B, TM
#3500 LET I3=0: REM BASS2.RES residual
#3510 LET Z1=Y(I,1)/Y(I,2): REM S/W
#3520 LET Z2=Y(I,5)/Y(I,2): REM TM/W
#3525 LET Z4=Y(I,3)/Y(I,4): REM t/B
#3530 LET Z3=B(3)/Z4+1: REM 1+GAMMA/(t/B)
#3540 LET R1=B(1)*Z1-B(2)*Z3*Z3*Z1*Z1-Z2
TMw<-bd$TM/bd$W
# W are weights on (both sides)
Sw<-bd$S/bd$W
tB<-bd$t/bd$B
bassdata<-data.frame(TMw=TMw, Sw=Sw, tB=tB)
resx<-"TMw ~ Sw*(b1 - b2*Sw*(1+b3/tB)^2)"
st1<-c(b1=1, b2=1, b3=1)
require(nlmrt)
a1<-nlsmnq(resx, start=st1, trace=TRUE, data=bassdata)
print(a1)
a1nls<-nls(resx, start=st1, trace=TRUE, data=bassdata)
print(a1nls)
