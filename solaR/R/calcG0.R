calcG0 <- function(lat,
                   modeRad='prom',
                   dataRad,
                   sample='hour',
                   keep.night=TRUE,
                   sunGeometry='michalsky',
                   corr, f, ...){

    if (missing(lat)) stop('lat missing. You must provide a latitude value.')

    stopifnot(modeRad %in% c('prom', 'aguiar','bd', 'bdI'))


###Datos de Radiacion
    if (missing(corr)){
        corr = switch(modeRad,
        bd = 'CPR', #Correlacion entre Fd y Kt para valores diarios
        aguiar = 'CPR', #Correlacion entre Fd y Kt para valores diarios
        prom = 'Page', #Correlacion entre Fd y Kt para promedios mensuales
        bdI = 'BRL' #Correlación entre fd y kt para valores intradiarios
        )
    }

    BD <- switch(modeRad,
                 bd = {
                     if (is(dataRad, 'Meteo')) {
                         dataRad
                     } else {
                         if (!is.list(dataRad)) dataRad <- list(file=dataRad)
                         switch(class(dataRad$file),
                                character={
                                    bd.default=list(file='', lat=lat, format="%d/%m/%Y",
                                    header=TRUE, fill=TRUE, sep=';',
                                    dec='.', dates.col='dates', source='')
                                    bd=modifyList(bd.default, dataRad)
                                    res <- do.call('readBD', bd)
                                    res
                                },
                                data.frame={
                                    bd.default=list(file='', lat=lat, format="%d/%m/%Y",
                                    dates.col='dates', source='')
                                    bd=modifyList(bd.default, dataRad)
                                    res <- do.call('df2Meteo', bd)
                                    res
                                },
                                zoo={
                                    bd.default=list(file='', lat=lat, source='')
                                    bd=modifyList(bd.default, dataRad)
                                    res <- do.call('zoo2Meteo', bd)
                                    res
                                })
                     }},                #Fin de bd
                 prom = {
                     if (!is.list(dataRad)) dataRad <- list(G0dm=dataRad)
                     prom.default <- list(G0dm=numeric(), Ta=25, lat=lat,
                                          year = as.POSIXlt(Sys.Date())$year+1900,
                                          promDays=c(17,14,15,15,15,10,18,18,18,19,18,13),
                                          source='')
                     prom = modifyList(prom.default, dataRad)
                     res <- do.call('readG0dm', prom)
                 },                     #Fin de prom
                 aguiar = {
                     if (is.list(dataRad)) dataRad <- dataRad$G0dm
                     BTd=fBTd(mode='serie')
                     solD <- fSolD(lat, BTd)
                     G0d <- markovG0(G0dm, solD)
                     res <- zoo2Meteo(G0d, lat=lat, source='aguiar')
                 },                     #Fin de aguiar
                 bdI = {
                     if (is(dataRad, 'Meteo')) {
                         dataRad
                     } else {
                         if (!is.list(dataRad)) dataRad <- list(file=dataRad)
                         switch(class(dataRad$file),
                                character = {
                                    bdI.default <- list(file='', lat=lat,
                                                        format="%d/%m/%Y %H:%M:%S",
                                                        header=TRUE, fill=TRUE, sep=';',
                                                        dec='.', time.col='time', source='')
                                    bdI <- modifyList(bdI.default, dataRad)
                                    res <- do.call('readBDi', bdI)
                                    res
                                },
                                data.frame = {
                                    bdI.default <- list(file='', lat=lat,
                                                        format="%d/%m/%Y %H:%M:%S",
                                                        time.col='time', source='')
                                    bdI <- modifyList(bdI.default, dataRad)
                                    res <- do.call('dfI2Meteo', bdI)
                                    res
                                },
                                zoo = {
                                    bdI.default <- list(file='', lat=lat, source='')
                                    bdI <- modifyList(bdI.default, dataRad)
                                    res <- do.call('zoo2Meteo', bdI)
                                    res
                                },
                                stop('dataRad$file should be a character, a data.frame or a zoo.')
                                )}}     #Fin de bdI
                 )                      #Fin del switch general


### Angulos solares y componentes de irradiancia
    if (modeRad=='bdI') {
        sol <- calcSol(lat=lat, BTi=index(getData(BD)),
                       keep.night=keep.night, method=sunGeometry)
        compI <- fCompI(sol=sol, G0I=BD, corr=corr, f=f, ...)
        compD <- aggregate(compI[,c('G0', 'D0', 'B0')],
                           by=truncDay, FUN=P2E, sol@sample) #Wh
        names(compD) <- c('G0d', 'D0d', 'B0d')
        compD$Fd <- compD$D0d/compD$G0d
        compD$Ktd <- compD$G0d/as.zooD(sol)$Bo0d
    } else { ##modeRad!='bdI'
        sol <- calcSol(lat=lat, BTd=indexD(BD),
                       sample=sample, keep.night=keep.night,
                       method=sunGeometry)
        compD <- fCompD(sol=sol, G0d=BD, corr=corr, f)
        compI <- fCompI(sol=sol, compD=compD, ...)
    }

###Temperatura

    ##Compruebo si tengo información de temperatura a partir de la cual
    ##generar una secuencia de datos. Para eso, debo estar leyendo de www.mapa.es
    ##o de una base de datos que contenga dos variables con información sobre
    ##valores diarios máximos y mínimos de temperatura.

    ind.rep <- indexRep(sol) ##para repetir valores diarios de Ta, si es necesario
    indSol <- indexI(sol)

    Ta=switch(modeRad,
    bd={
        if (all(c("TempMax","TempMin") %in% names(BD@data))) {
            fTemp(sol, BD)
        } else {
            if ("Ta" %in% names(BD@data)) {
                zoo(coredata(BD@data$Ta)[ind.rep], indSol)
            } else {
                warning('No temperature information available!')
            }
        }
    },
    bdI={
        if ("Ta" %in% names(BD@data)) {
            Ta=BD@data$Ta
        } else {
            warning('No temperature information available!')
        }
    },
    prom <- zoo(coredata(BD@data$Ta)[ind.rep], indSol) ##zoo(rep(Ta, length(indSol)), indSol) ##idem
    )

###Medias mensuales y anuales
    DayOfMonth=c(31,28,31,30,31,30,31,31,30,31,30,31) ###OJO

    G0dm <- aggregate(compD[,c('G0d', 'D0d', 'B0d')], by=as.yearmon,
                      FUN=function(x, ...)mean(x, na.rm=1)/1000) ##kWh
    if (modeRad=='prom'){
        G0y <- zoo(t(colSums(G0dm*DayOfMonth)),
                   unique(year(index(G0dm))))
    } else {
        G0y <- aggregate(compD[,c('G0d', 'D0d', 'B0d')], by=year,
                         FUN=function(x, ...)sum(x, na.rm=1)/1000) ##kWh
    }

###Resultado
    result <- new(Class='G0',
                  BD,                     #G0 contains "Meteo"
                  sol,                    #G0 contains 'Sol'
                  G0D=compD,              #resultado de fCompD
                  G0dm=G0dm,          #aggregate, medias mensuales
                  G0y=G0y,            #aggregate, valores anuales
                  G0I=compI,          #resultado de fCompI
                  Ta=Ta               #temperatura ambiente
                  )
    return(result)
}
