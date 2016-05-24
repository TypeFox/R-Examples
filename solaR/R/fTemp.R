fTemp<-function(sol, BD)
{
    ##sol es un objeto con class='Sol'
    ##BD es un objecto con class='Meteo', cuyo slot 'data' contiene dos columnas llamadas "TempMax" y "TempMin"

    stopifnot(class(sol)=='Sol')
    stopifnot(class(BD)=='Meteo')

    checkIndexD(indexD(sol), indexD(BD))

    indSol<-indexI(sol)	
    ind.rep<-indexRep(sol)
    
    TempMax=coredata(BD@data$TempMax)[ind.rep]
    TempMin=coredata(BD@data$TempMin)[ind.rep]
    ws=coredata(sol@solD$ws)[ind.rep]
    w=coredata(sol@solI$w)

    ## Genera secuencia de temperatura a partir de Maxima y Minima de base de datos

    Tm=(TempMin+TempMax)/2
    Tr=(TempMax-TempMin)/2

    wp=pi/4

    a1=pi*12*(ws-w)/(21*pi+12*ws)
    a2=pi*(3*pi-12*w)/(3*pi-12*ws)
    a3=pi*(24*pi+12*(ws-w))/(21*pi+12*ws)

    T1=Tm-Tr*cos(a1)
    T2=Tm+Tr*cos(a2)
    T3=Tm-Tr*cos(a3)

    Ta=T1*(w<=ws)+T2*(w>ws&w<=wp)+T3*(w>wp)

    ## Resultado
    result<-zoo(Ta, indSol)
}
			
