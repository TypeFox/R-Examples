## ------------------------------------------------------------------------
library("RRTCS")
N=10777
n=710
data(HorvitzDataRealSurvey)
p=0.5
alpha=c(1/12,1/10,20/30,1/10,10/30,1/12)
pi=rep(n/N,n)
cl=0.95
out1=Horvitz(HorvitzDataRealSurvey$copied,p,alpha[1],pi,"mean",cl,N)
out1
out2=Horvitz(HorvitzDataRealSurvey$fought,p,alpha[2],pi,"mean",cl,N)
out2
out3=Horvitz(HorvitzDataRealSurvey$bullied,p,alpha[3],pi,"mean",cl,N)
out3
out4=Horvitz(HorvitzDataRealSurvey$bullying,p,alpha[4],pi,"mean",cl,N)
out4
out5=Horvitz(HorvitzDataRealSurvey$drug,p,alpha[5],pi,"mean",cl,N)
out5
out6=Horvitz(HorvitzDataRealSurvey$sex,p,alpha[6],pi,"mean",cl,N)
out6

