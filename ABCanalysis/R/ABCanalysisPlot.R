ABCanalysisPlot=function(Data,LineType=0,LineWidth=3,ShowUniform=TRUE,title='ABC analysis',limits=TRUE,MarkPoints=TRUE,ABCcurvedata,ResetPlotDefaults=TRUE){
# res= ABCanalysisPlot(Data=ABCcleanData(Data)$CleanedData,style='2')
# display ABC Curve : cumulative percentage of largest Data (Effort) vs cumlative percentage of sum of largest Data (Yield)
# 
# INPUT
# Data(1:n)          oder [Haeufigkeit(1:n),Data(1:n)]  or
#                    Data = [ABCx,ABCy] iff ABCx(1) ==0 
#
# OPTIONAL
# LineType           for plot default:  LineType=0 for Line, other numbers see documentation about pch
# LineWidth          Breite der ABC kurve 
# ShowUniform        ==1 (default) bedeutet  die ABC kurve der Uniform verteilung Unifirm[0,beliebig] wird eingezeichnet
# title              string, label for the title  of the plot  
# limits             =TRUE Linien zur Einteilung werden gezeichnet, Default=FALSE
# MarkPoints         MarkPoints=True => MarkPointsOfInterest, Default=TRUE
# ABCcurvedata                Liste V aus ABCcurve()
#
# OUTPUT
# ABC                   Output von ABCplot
# A=c(Ax,Ay)            A Point: Minimum distance to (0,1) = minimal Unrealized potential == min in Effort und min in (1-
# ABCanalysis           Liste V aller Daten aus ABCanalysis()
# 
# author: MT 11/2014
# 1.Editor: MT 01/2015
# Nota: Diese Funktion ist eine "Faulheitsfunktion" => auch Warnemldungen werden unterdrueckt
#suppressWarnings(require(plotrix))
style=TRUE
            #farb.col=c('black','red','blue','green',colors()[452],colors()[57])
     farb.col=c('black','red','blue','green',colors()[452],'red')
     
        #  farb.col=c('paleturquoise3','magenta','blue','palegreen4','palegreen3','plum2)
      #farb.labels <- c(expression(italic("Equilibrium")),expression(italic("set limits")),expression(italic("data")),expression(italic("uniform")),expression(italic("identity")))
      farb.labels <- c('',expression(italic("set limits")),expression(italic("data")),expression(italic("uniform")),expression(italic("identity")))
      
  
if(missing(style)){style=FALSE}
  if(missing(Data)){
    Data=NULL #Faulheitsfunktion, auch Listen?bergabe moeglich
    curve=ABCcurvedata
  }else{
     curve = suppressWarnings(ABCcurve(Data))
  }
   def.par <- par(no.readonly = TRUE) # save default, for resetting...
    if(style==FALSE){

          abc=suppressWarnings(ABCplot(Data,LineType=LineType,LineWidth=LineWidth,ShowUniform=ShowUniform,title=title,defaultAxes=FALSE,ABCcurvedata=ABCcurvedata))

   }else{#Equilibrium, set limits, data, uniform, identity, ABC-Guppen Buchstaben

      #if(missing(title)) title='ABC plot for data grouping'
      abc=suppressWarnings(ABCplot(Data,LineType=LineType,LineWidth=LineWidth,ShowUniform=ShowUniform,title=title,defaultAxes=FALSE,ABCcurvedata=ABCcurvedata))
    }   
    axis(1,xlim=c(0,1),col="black",las=1, at=seq(from=0,to=1,by=0.1)) #x-Achse
    axis(2,ylim=c(0,1),col="black",las=1, at=seq(from=0,to=1,by=0.1)) #y-Achse
 
  abcres = ABCanalysis(Data,ABCcurvedata=curve) #Achtung hier darf PlotIt nicht angegeben werden!
  
  if(MarkPoints){
   if(style==FALSE){
    points(abcres$A[1],abcres$A[2],pch=8,lwd=1.5,col='red',cex=1.5,asp=1)
    points(abcres$B[1],abcres$B[2],pch=8,lwd=1.5,col='green',cex=1.5,asp=1)
    #text(abcres$A[1],abcres$A[2]+0.1,labels=paste0('Ax=',round(abcres$A[1],2)),asp=1)
    points(abcres$C[1],abcres$C[2],pch=8,lwd=1.5,col='blue',cex=1.5)
    #text(abcres$C[1],abcres$C[2]+0.1,labels=paste0('Cx=',round(abcres$C[1],2)),asp=1)
    }
  }
    nA=length(abcres$Aind)
    nB=length(abcres$Bind)
    nC=length(abcres$Cind)
if(style==FALSE){ 
  if(limits){
    if(!MarkPoints){abcres = ABCanalysis(Data,ABCcurvedata=curve)} #Achtung hier darf PlotIt nicht angegeben werden

    #points(c(abcres$A[1],abcres$A[1]),c(0,abcres$C[2]),type='l',col='red') 
   # points(c(0,abcres$C[1]),c(abcres$A[2],abcres$A[2]),type='l',col='red')
    #points(c(abcres$C[1],abcres$C[1]),c(0,1),type='l',col='red')  
   # points(c(0,1),c(abcres$C[2],abcres$C[2]),type='l',col='red')
   linientyp=1 #Linie
     #Linien y-achse set limits  
    points(c(0,abcres$C[1]),c(abcres$C[2],abcres$C[2]),type='l',col=farb.col[2],lty=linientyp) 
    points(c(0,abcres$A[1]),c(abcres$A[2],abcres$A[2]),type='l',col=farb.col[2],lty=linientyp) 
 
    #Linien x-achse set limits   
    points(c(abcres$A[1],abcres$A[1]),c(0,abcres$A[2]),col=farb.col[2],type='l',lty=linientyp)
    points(c(abcres$C[1],abcres$C[1]),c(0,abcres$C[2]),col=farb.col[2],type='l',lty=linientyp)  
    
   
  if(!is.null(Data)){
      if(abs(abcres$A[1]-abcres$C[1])>0.1){
        thigmophobe.labels(x=abcres$A[1]/2,y=abcres$A[2],paste0('A:n=',nA),col='black', cex=1) # Gr??e
        thigmophobe.labels(x=(abcres$C[1]-abcres$A[1])/2+abcres$A[1],y=abcres$C[2],paste0('B:n=',nB),col='black',cex=1)
      
      }else{
        thigmophobe.labels(x=abcres$A[1]-0.05,y=abcres$A[2],paste0('A:n=',nA),col='black', cex=1) # Gr??e
        thigmophobe.labels(x=abcres$C[1]+0.025,y=abcres$C[2]+0.025,paste0('B:n=',nB),col='black',cex=1)
      
      }
     thigmophobe.labels(x=(1-abcres$C[1])/2+abcres$C[1],y=abcres$C[2],paste0('C:n=',nC),col='black',cex=1)
  }  
  }
}else{ 

 box(col='grey')

  #Grenzpunkte
  points(abcres$A[1],abcres$A[2],pch=8,lwd=1.5,col=farb.col[2],cex=1.5,asp=1)
  points(abcres$B[1],abcres$B[2],pch=8,lwd=1.5,col=farb.col[4],cex=1.5,asp=1)
  points(abcres$C[1],abcres$C[2],pch=8,lwd=1.5,col=farb.col[3],cex=1.5)

  #Grenzpunkte Zeichen
  
    if(abs(abcres$A[1]-abcres$C[1])>0.1){
      thigmophobe.labels(x=abcres$A[1],y=abcres$A[2],'A|B',col=farb.col[2], cex=1) # Gr??e
      thigmophobe.labels(x=abcres$C[1],y=abcres$C[2],'B|C',col=farb.col[2],cex=1)
    
    }else{
      thigmophobe.labels(x=abcres$A[1]-0.05,y=abcres$A[2],'A|B',col=farb.col[2], cex=1) # Gr??e
      thigmophobe.labels(x=abcres$C[1]+0.025,y=abcres$C[2]+0.025,'B|C',col=farb.col[2],cex=1)
    
    }
  #ABC gruppen (Buchstaben)
    thigmophobe.labels(x=abcres$A[1]/2,y=abcres$A[2]/4,'A',col=farb.col[6], cex=2.6) # Gr??e 
    thigmophobe.labels(x=(abcres$A[1]+abcres$C[1])/2,y=abcres$A[2]/4,'B',col=farb.col[6],cex=2.1)
    thigmophobe.labels(x=(abcres$A[1]+abcres$C[1])/2+max(abs(abcres$A[1]-abcres$C[1]),0.1),y=abcres$A[2]/4,'C',col=farb.col[6],cex=1.8)
if(!is.null(Data)){
   thigmophobe.labels(x=abcres$A[1]/2,y=abcres$A[2]/4-0.05,paste0('n=',nA),col='black', cex=0.8) # Gr??e
   thigmophobe.labels(x=(abcres$A[1]+abcres$C[1])/2,y=abcres$A[2]/4-0.05,paste0('n=',nB),col='black',cex=0.8)
   thigmophobe.labels(x=(abcres$A[1]+abcres$C[1])/2+max(abs(abcres$A[1]-abcres$C[1]),0.1)+0.02,y=abcres$A[2]/4-0.05,paste0('n=',nC),col='black',cex=0.8)
}    
   # detach(package:plotrix)
  #Markierungen y-Achse
    points(c(0,0.01),c(abcres$C[2],abcres$C[2]),type='l',col=farb.col[2],lwd=2) 
    points(c(0,0.01),c(abcres$A[2],abcres$A[2]),type='l',col=farb.col[2],lwd=2) 
  #linientyp=5 #Gestrichelt
  linientyp=1 #Linie
  #Linien y-achse set limits  
    points(c(0,abcres$C[1]),c(abcres$C[2],abcres$C[2]),type='l',col=farb.col[2],lty=linientyp) 
    points(c(0,abcres$A[1]),c(abcres$A[2],abcres$A[2]),type='l',col=farb.col[2],lty=linientyp) 
   #Markierungen x-Achse   
    points(c(abcres$C[1],abcres$C[1]),c(0,0.01),type='l',col=farb.col[2],lwd=2) 
    points(c(abcres$A[1],abcres$A[1]),c(0,0.01),type='l',col=farb.col[2],lwd=2) 
   #Linien x-achse set limits   
    points(c(abcres$A[1],abcres$A[1]),c(0,abcres$A[2]),col=farb.col[2],type='l',lty=linientyp)
    points(c(abcres$C[1],abcres$C[1]),c(0,abcres$C[2]),col=farb.col[2],type='l',lty=linientyp)  
    
    #diagonale: Skewness bzw Gleichgewichtspunkt
#     pUnif = seq(from=0,by=0.001,to=1)
#     A = min(Data,na.rm=TRUE) 
#     if(!is.null(Data)){
#     A = min(Data,na.rm=TRUE) 
#     MaxX = max(Data,na.rm=TRUE)
#     }else{
#       A=0
#       MaxX=1
#     }
#     Bmax = MaxX-A
#     ABCuniform = (-0.5*Bmax*pUnif^2+MaxX*pUnif)/(A+0.5*Bmax)
#     ind=which.min(abs(pUnif-(1-ABCuniform)))
#    # points(c(0,1),c(1,0),type='c',ljoin='mitre',col=farb.col[1])

}
#damit Buchtabe C und Legende sich nicht ueberschneiden
if((abcres$A[1]+abcres$C[1])/2+max(abs(abcres$A[1]-abcres$C[1]),0.1)+0.02<0.8){
  legend('bottomright',legend=farb.labels,text.col=farb.col,bty = "n",y.intersp=0.8)
}else{
    legend('right',legend=farb.labels,text.col=farb.col,bty = "n",y.intersp=0.8)

}
if(ResetPlotDefaults)
    par(def.par)
   
  invisible(list(ABC=abc,ABCanalysis=abcres))
}