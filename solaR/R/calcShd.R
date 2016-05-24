calcShd<-function(radEf,##class='Gef'
                  modeTrk='fixed',     #c('two','horiz','fixed')
                  modeShd='prom',      #modeShd=c('area','bt','prom')
                  struct=list(), #list(W=23.11, L=9.8, Nrow=2, Ncol=8), 
                  distances=data.frame() #data.frame(Lew=40, Lns=30, H=0)){
                  )
{
    stopifnot(is.list(struct), is.data.frame(distances))

    ##Por ahora sólo utilizo modeShd='area'
    ##Con diferentes modeShd (por definir) podré calcular Gef de diferente forma
    ##Ver tesis de Macagnan
    prom=("prom"  %in%  modeShd)
    prev <- as.zooI(radEf, complete=TRUE)
    ## Cálculo de sombras
    sol <- prev[,c('AzS', 'AlS')]
    theta <- prev[,c('Beta', 'Alfa', 'cosTheta')]
    AngGen <- CBIND(sol, theta, index=indexI(radEf))
    Shd<-fSombra(AngGen, distances, struct, modeTrk, prom)
    FS=coredata(Shd)
    ## Cálculo de irradiancia
    ## gef0 <- prev[,c('G', 'D', 'B', 'R',
    ##              'FTb', 'FTd', 'FTr',
    ##              'Dief', 'Dcef',
    ##              'Gef', 'Def', 'Bef', 'Ref')]
    gef0 <- radEf@GefI
    Bef0=gef0$Bef
    Dcef0=gef0$Dcef
    Gef0=gef0$Gef
    Dief0=gef0$Dief
    Ref0=gef0$Ref
    ##Cálculos
    Bef=Bef0*(1-FS)
    Dcef=Dcef0*(1-FS)
    Def=Dief0+Dcef
    Gef=Dief0+Ref0+Bef+Dcef               #Incluyendo sombras
    ##Cambio nombres
    nms=c('Gef', 'Def', 'Dcef', 'Bef')
    nmsIndex <- which(names(gef0) %in% nms)
    names(gef0)[nmsIndex]<- paste(names(gef0)[nmsIndex], '0', sep='')
    ##Entrego zoo con resultados, incluyendo previos sin sombras
    GefShd=CBIND(gef0, data.frame(Gef, Def, Bef, Dcef, FS), index=indexI(radEf))


    ## Valores diarios, mensuales y anuales
    DayOfMonth=c(31,28,31,30,31,30,31,31,30,31,30,31) ## OJO
    
    if (radEf@type=='prom') {
        Gefdm=aggregate(GefShd[,c('Gef0', 'Def0', 'Bef0', 'G', 'D', 'B', 'Gef', 'Def', 'Bef')]/1000,
                        by=as.yearmon, FUN=P2E, radEf@sample)       #kWh
        names(Gefdm)=paste(names(Gefdm), 'd', sep='')

        GefD=Gefdm*1000                  #Wh
        index(GefD) <- indexD(radEf)     ##para que sea compatible con G0D
        
        Gefy=zoo(t(colSums(Gefdm*DayOfMonth)),
                 unique(year(index(Gefdm))))
    } else {
        GefD=aggregate(GefShd[,c('Gef0', 'Def0', 'Bef0', 'G', 'D', 'B', 'Gef', 'Def', 'Bef')],
                       by=truncDay, FUN=P2E, radEf@sample)         #Wh
        names(GefD)=paste(names(GefD), 'd', sep='')

        Gefdm=aggregate(GefD/1000, by=as.yearmon, mean, na.rm=1)
        Gefy=aggregate(GefD/1000, by=year, sum, na.rm=1)
    }

    ## Entrego un objecto de clase Gef
    ## modificando los slots 'modeShd', 'GefI', 'GefD', 'Gefdm', y 'Gefy'
    ## del objecto original radEf
    radEf@modeShd=modeShd
    radEf@GefI=GefShd
    radEf@GefD=GefD
    radEf@Gefdm=Gefdm
    radEf@Gefy=Gefy
    return(radEf)
}

  
