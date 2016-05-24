ABCplot=function(Data,LineType=0,LineWidth=3,ShowUniform=TRUE,title,ABCcurvedata,defaultAxes=TRUE){
# res= ABCplot(Data)
# display ABC Curve : cumulative percentage of largest Data (Effort) vs cumlative percentage of sum of largest Data (Yield)
# 
# INPUT
# Data(1:n)          oder [Haeufigkeit(1:n),Data(1:n)]  or
#                    Data = [ABCx,ABCy] iff ABCx(1) ==0 
#
# OPTIONAL
# LineType           for plot default:  LineType=0 for Line, other numbers see documentation about pch
# LineWidth          Breite der ABC kurve
#
# ShowUniform    ==1 (default) bedeutet  die ABC kurve der Uniform verteilung Unifirm[0,beliebig] wird eingezeichnet
# title             string, label for the title  of the plot 
# style             type fancy if you would like to plot in an different style
# ABCcurvedata             Input form ABCcurve
# defaultAxes       FALSE
#  
# OUTPUT
# ABCx                  cumulative population in Percent
# ABCy                  cumulative high Datas in Percent
# author: MT 11/2014
# 1.Editor MT 01/2015
style=TRUE
  #require(Hmisc) #Noch anders zu loesen
  check=F
  if(missing(Data)){ 
  curve = ABCcurvedata
  Data=NULL
  check=T
 }
 if(is.null(Data)&!check){ 
  curve = ABCcurvedata
 }
if(missing(title)){title='ABC plot'}
 if(missing(ABCcurvedata)){      
   curve = ABCcurve(Data)}


Effort=curve$Curve[,'Effort']
Yield=curve$Curve[,'Yield']
par(pty="s")# Plot immer quadratisch

if(missing(style)){
  ylab='fraction of sum of largest data'
  xlab='fraction of data'
  farb.col=c('blue',colors()[452],'green',colors()[175])
  farb.labels <- c(expression(italic("data")),expression(italic("identity")),expression(italic("uniform")),'')

}else{

  ylab='fraction of sum of largest data'
  xlab='fraction of data'
  farb.col=c('blue',colors()[452],'green',colors()[175])

  farb.labels <- c(expression(italic("data")),expression(italic("identity")),expression(italic("uniform")),expression(italic("equilibrium")))

}

if(!ShowUniform){#Dann ist dieser Plot im Vordergrund
      #farb.labels=farb.labels[c(1,3)]
    if(LineType==0){
        plot(Effort,Yield,xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',xlab=xlab,ylab=ylab,type='l',lwd=LineWidth,col=farb.col[1],main=title,axes=defaultAxes)
        }else{
        plot(Effort,Yield,xlim=c(0,1),ylim=c(0,1),asp=1,xaxs='i',yaxs='i',xlab=xlab,ylab=ylab,pch=LineType,lwd=LineWidth,col=farb.col[1],main=title,axes=defaultAxes)
      }
}else{

    #gleichverteilung
    pUnif = seq(from=0,by=0.01,to=1)
    if(!is.null(curve$CleanedData)){
    A = min(curve$CleanedData,na.rm=TRUE) 
    MaxX = max(curve$CleanedData,na.rm=TRUE)
      if(A==MaxX){
        A=0
        MaxX=1
      }
    }else{
      A=0
      MaxX=1
    }

    B = MaxX-A
    ABCuniform = (-0.5*B*pUnif^2+MaxX*pUnif)/(A+0.5*B)
    if(missing(style)){

        plot(pUnif,ABCuniform,type='l',col=farb.col[3],asp=1,xaxs='i',yaxs='i',xlab=xlab,ylab=ylab,axes=defaultAxes,main=title) 
        points(c(0,1),c(1,0),type='l',lty=2,lwd=1,col=farb.col[4],asp=1) #diagonale
    }else{
      plot(pUnif,ABCuniform,type='l',col=farb.col[3],asp=1,lwd=1,xaxs='i',yaxs='i',xlab=xlab,ylab=ylab,axes=defaultAxes,main=title) 
      points(c(0,1),c(1,0),type='l',lty=2,lwd=1,col=farb.col[4],asp=1) #diagonale
      #points(c(0,1),c(1,0),type='l',col=colors()[234],asp=1) #diagonale  
    }    
    if(LineType==0){
      points(Effort,Yield,xlim=c(0,1),ylim=c(0,1),lwd=LineWidth,col=farb.col[1],main=title,type='l')
    }else{
      points(Effort,Yield,xlim=c(0,1),ylim=c(0,1),pch=LineType,lwd=LineWidth,col=farb.col[1],main=title,type='l')
    }
}

ableitung=curve$Slope[,'dABC']
# Suche das Minimum zur Differenz Ableitung vordefiniete Steigung
BreakEvenInds=which.min(abs(ableitung - 1))  #Schraenkt B ein, Ableitung==1

# Bestimme die Werte der Spline Kurve an dem BC Punkt 
BreakEvenInd=max(BreakEvenInds)# falls es mehr als 1 gibt nimm den rechtesten
Kurve=cbind(Effort,Yield)
BreakEvenPoint=Kurve[BreakEvenInd,] #Last
points(BreakEvenPoint[1],BreakEvenPoint[2],pch=8,lwd=1.5,col='green',cex=1.5,asp=1)

if(!is.null(Data)){
if(length(curve$CleanedData)<20){
  
  sorted=sort(curve$CleanedData,decreasing=TRUE)
  Anteil=sorted
  y=cumsum(Anteil)
  y=y/tail(y,1)
  x=(1:length(curve$CleanedData))/length(curve$CleanedData)
  
  points(x,y,pch=1,lwd=1.5,col='blue',cex=1.5,asp=1)
}
}

if(missing(style)){
        points(Effort,Effort,type='l',lwd=1,col=farb.col[2],asp=1) #idetitaet 
}else{
          points(Effort,Effort,type='l',lwd=0.1,col=farb.col[2],asp=1) #identitaet
}
if(defaultAxes){
  axis(1, at=seq(from=0,to=1,by=0.1)) 
  axis(2, at=seq(from=0,to=1,by=0.1)) 

  if(defaultAxes){
        legend("bottomright",bty = "n",legend=farb.labels,text.col=farb.col)
  }  
}
if(!missing(style)){
  #requireRpackage('Hmisc')
        minor.tick(ny=20, nx=20)
        box()
        }else{
          minor.tick(ny=20, nx=20)
          box(col='grey')
        }

invisible(list(ABCx=Effort,ABCy=Yield))
}