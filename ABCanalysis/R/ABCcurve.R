 ABCcurve=function(Data,p){
# res = ABCcurve(Data,GiniSteigung)
# ABC Curve : cumulative fraction of largest Data in population vs fraction of population
#
# INPUT
# Data(1:n)          data vector,only positive data will be used
# 
# OPTIONAL
#              
# p                  x-werte fuer Spline Interpolation: wenn vorgegeben dann werden diese genommen
# 
# OUTPUT
# ABCx(1:n)           fraction of population in [0,1]
# ABCy(1:n)           cumulative fraction of largest Datas in [0,1]
#
# if nargin  < 3; p = (0:0.01:1)';       end; % x-werte fuer Spline Interpolation       
#
#author: MT 11/2014
# 1.Editor MT 01/2015

cleanData=ABCcleanData(Data)$CleanedData
rows=length(cleanData)
if(missing(p)){
  if(rows<101){ p=seq(from=0,to=1,by=0.01)
  }else{ p=seq(from=0,to=1,by=0.001)}
}   

sorted=sort(cleanData,decreasing=TRUE)
#N=sum(cleanData)
#Anteil=sorted/N
Anteil=sorted
y=cumsum(Anteil)
y=y/tail(y,1)
x=(1:rows)/rows
## Die Kurve muss durch 2 Punkte gehen 0 und 1

 if(head(y,1)>0){
   x=c(0,x)
   y=c(0,y)
 }
 if(tail(x,1)<1){ #Nach matlab Implementation, ueberfluessig?
   x=c(x,1)
   y=c(y,1)
 }

## Spline Interpolation

V=spline(x,y,xout=p)
Effort=V$x
Yield=V$y
#Fehlerabfang der Interpolation
inds=which(Yield>=1)
ind1=min(inds)
if(ind1<length(Yield))
  Yield[c(ind1:length(Yield))]=1

n=length(Effort)
Curvengleichung=splinefun(Effort,Yield) 
ableitung=Curvengleichung(1:n/n,1)

return(list(Curve=cbind(Effort=Effort,Yield=Yield),CleanedData=cleanData,Slope=cbind(p=p,dABC=ableitung)))
}
