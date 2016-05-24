ABCanalysis=function(Data,PlotIt,ABCcurvedata){
# abcres = ABCanalysis(Data=ABCcleanData(Data)$CleanedData)
# divide the Data in 3 classes A, B,C
# A==Data(Aind) : mit wenig aufwand viel ertrag! 
# B==Data(Bind) : Aufwand und ertrag halten sich die Waage
# C==Data(Cind) : viel Aufwand, wenig Ertrag
#
# Grenzziehung AB: minimaler Abstand zum ideal [0,1]
# Grenzziehung BC: 
# 
# INPUT
# Data(1:n)                 Ungleichverteilung so dass ABC Analyse sinnvoll ist
#
# OPTIONAL
# PlotIt                      if variable is used a plot is made, set with PlotIt=1,PlotIt=TRUE, PlotIt='On', etc
# ABCcurvedata                Liste V aus ABCcurve()
#
# OUTPUT
# Aind,Bind,Cind              so dass:
#                             A==Data(Aind) : mit wenig aufwand viel ertrag! 
#                             B==Data(Bind) : Aufwand und ertrag halten sich die Waage
#                             C==Data(Cind) : viel Aufwand, wenig Ertrag
# smallestAData: Grenzziehung AB: minimaler Abstand zum ideal [0,1]
# smallestBData: Grenzziehung BC: Steigung der ABC kurve == 1  
#
# AlimitIndInInte\rpolation,BlimitIndInInterpolation    indices der ABC genzen in [p,ABC]
# [p,ABC]                                              die interpolationskurve des ABC plots.
# 
# A=c(Ax,Ay)              Pareto point, Minimum distance to (0,1) = minimal Unrealized potential == min in Effort und min in (1-
# B=c(Bx,By)              BreakEven Point: dABC(Bx) == 1 
# C=c(Cx,Cy)             Submarginal Point: Minimum distance to (Bx,1) 
# ABexchanged             bool, TRUEif Point A is the BreakEven and point B is the Pareto Point, 0 otherwise

  
# author: MT 11/2014
#

#Uses:
# ABCcurve, im else-Fall: ABCanalysisPlot
if(missing(Data)){
  if(missing(ABCcurvedata)){stop('argument "Data" and ABCcurvedata are missing')}else{
  Data=NULL #Dann muss ABCcurvedata vorhanden sein
  }
}
  
if(missing(PlotIt)){#Wenn nicht geplottet wird muessen die Daten berechnet werden
  if(missing(ABCcurvedata)){
    # CleanData=1; # bedeutet die kleinen Yieldwerte, die in Summme <0.5% der Gesamtyield 
     ABCcurvedata = ABCcurve(Data)
  } 
  
  Effort=ABCcurvedata$Curve[,'Effort']
  Yield=ABCcurvedata$Curve[,'Yield']

#   Indizies=ABCcurvedata$DataInd
#Distanz zum 0,1 Punkt ueber Euklid berechnet
  curve=cbind(Effort,Yield)
  distPareto=c()
  point=t(as.matrix(c(0,1)))
  for(i in 1:length(Effort)){
    distPareto[i]=sum(abs(point-curve[i,])^2)
  }
  ParetoPointInd=which.min(distPareto) #First Minimun
# Bestimme den AB Punkt auf der Spline Kurve  
  ParetoPoint=curve[ParetoPointInd,]
# Fuer den Punkt BC bestimme die Ableitung der Kurve 

#   n=length(Effort)
#   Curvengleichung=splinefun(Effort,Yield) 
#   ableitung=Curvengleichung(1:n/n,1)

   ableitung=ABCcurvedata$Slope[,'dABC']
# Suche das Minimum zur Differenz Ableitung vordefiniete Steigung
  BreakEvenInds=which.min(abs(ableitung - 1))  #Schraenkt B ein, Ableitung==1

# Bestimme die Werte der Spline Kurve an dem BC Punkt 
BreakEvenInd=max(BreakEvenInds)# falls es mehr als 1 gibt nimm den rechtesten
  BreakEvenPoint=curve[BreakEvenInd,] #Last
  
if(Effort[BreakEvenInd]<Effort[ParetoPointInd]){
    ABexchanged=TRUE
    JurenInd=BreakEvenInd
    Bx=Effort[ParetoPointInd]
   # By=Yield[BreakEvenInd]
	A=BreakEvenPoint
	B=ParetoPoint
  }else{
    JurenInd=ParetoPointInd
    Bx=Effort[BreakEvenInd]
    #By=Yield[ParetoPointInd]
    ABexchanged=FALSE
	A=ParetoPoint
	B=BreakEvenPoint
  }

  distBx=c()
  Juren=t(as.matrix(c(Bx,1)))
  for(i in 1:length(Effort)){
    distBx[i]=sum(abs(Juren-curve[i,])^2)
  }
  bgrenze=which.min(distBx) #First Minimun 
#print(curve[bgrenze[1],])
  
C=curve[bgrenze[1],]
## Datenvektor in 3 Gruppen Teilen
if(!is.null(Data)){ 
#Statt nach in Y-Werten suchen, suchen wir in x-Werten, da diese Eindeutig
#y-Werte koennen dagegen in spezialfaellen mehrfach belegt sein
    Indizies=order(Data,decreasing=TRUE)
    rows=length(Data)
    x=1:rows/rows
# Vergleiche empirische Kurve mit generierten Ma? aus theoretischer Spline Kurve
    Aindvor=which(x<A[1],arr.ind=TRUE)# Suche alle Indizes bis zur ABGrenze
    ABind=which(x<C[1],arr.ind=TRUE) # Suche alle Indizes bis zur BCgrenze
#Setzte in unsortierten Datenvektor die Indizes
    Bind=Indizies[setdiff(ABind,Aindvor)] #Bindizes sind Differenz aus den beiden Zeilen davor
    Aind=Indizies[Aindvor] 
    Cind=Indizies[which(x>C[1],arr.ind=TRUE)]

}else{ #Keine empirischer Datensatz sondern vorgabe theoretischer Kurve
  Bind=NULL
  Cind=NULL
  Aind=NULL
  warning('No Data given: Calculating curve and points by given ABCcurvedata')
}
  return(list(Aind=Aind,Bind=Bind,Cind=Cind,ABexchanged=ABexchanged,A=A,B=B,C=C,smallestAData=Yield[JurenInd],smallestBData=Yield[bgrenze],AlimitIndInInterpolation=JurenInd,BlimitIndInInterpolation=bgrenze,p=Effort,ABC=Yield))
  #Falls Plot erwuenscht
  }else{ #Dann werden die Daten in ueber ABCanalysisPlot berechnet und ABCplot verwendet
    if(missing(Data)|is.null(Data)){
      abc=ABCanalysisPlot(ABCcurvedata=ABCcurvedata)$ABCanalysis  
    }else{
      abc=ABCanalysisPlot(Data)$ABCanalysis
    }
  }  
}