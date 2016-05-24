setMethod('show', 'Meteo',
          function(object){
            cat('Object of class ', class(object),'\n\n')
            cat('Source of meteorological information: ')
            cat(paste(object@type, object@source, sep='-'),'\n')
            cat('Latitude of source: ',
                paste(round(getLat(object,'deg'), 1), 'degrees\n\n'))
            cat('Meteorological Data:\n')
            print(summary(getData(object)))
          }
          )

header <-function(object){
  cat('Object of class ', class(object),'\n\n')
  cat('Source of meteorological information: ')
  cat(paste(object@type, object@source, sep='-'),'\n\n')
  cat('Latitude of source: ',
      paste(round(getLat(as(object, 'Meteo'),'deg'), 1), 'degrees\n'))
  cat('Latitude for calculations: ',
      paste(round(getLat(object, 'deg'),1), 'degrees\n\n'))
}

setMethod('show', 'Sol',
          function(object){
            cat('Object of class ', class(object),'\n\n')
            cat('Latitude:',
                paste(round(getLat(object, 'deg'),1), 'degrees\n\n'))
            cat('Daily values:\n')
            print(summary(object@solD))
            cat('\nIntradaily values:\n')
            print(summary(object@solI))
          }
          )

setMethod('show', 'G0',
          function(object){
            header(object)
            cat('Monthly averages:\n')
            print(as.zooM(object))
            cat('\nYearly values:\n')
            print(as.zooY(object))          }
          )

## setMethod('show', 'Gef',
##           function(object){
##             header(object)
##             cat('Monthly averages (kWh/m\u00b2):\n')
##             print(object@Gefdm)
##             cat('\nYearly values (kWh/m\u00b2):\n')
##             print(object@Gefy)
##           }
##           )

setMethod('show', 'Gef',
          function(object){
            callNextMethod()
            cat('-----------------\n')
            cat('Mode of tracking: ', object@modeTrk,'\n')
            if (object@modeTrk=='fixed'){
              cat('    Inclination: ', object@angGen$beta, '\n')
              cat('    Orientation: ', object@angGen$alfa, '\n')
            } else {
              cat('    Inclination limit:', object@angGen$betaLim, '\n')
            }
            ## cat('Monthly averages (kWh/kWp):\n')
            ## print(object@prodDm)
            ## cat('\nYearly values (kWh/kWp):\n')
            ## print(object@prody)
          }
          )

setMethod('show', 'ProdGCPV',
          function(object){
            callNextMethod()
            cat('-----------------\n')
            cat('Generator:\n')
            cat('    Modules in series: ', object@generator$Nms, '\n')
            cat('    Modules in parallel: ', object@generator$Nmp, '\n')
            cat('    Nominal power (kWp): ',
                round(object@generator$Pg/1000, 1), '\n\n')

            ## cat('Monthly averages (kWh/kWp):\n')
            ## print(object@prodDm)
            ## cat('\nYearly values (kWh/kWp):\n')
            ## print(object@prody)
          }
          )

setMethod('show', 'ProdPVPS',
          function(object){
            callNextMethod()
            cat('-----------------\n')
            cat('Pump:\n')
            cat('    Qn: ', object@pump$Qn, '\n')
            cat('    Stages: ', object@pump$stages, '\n')
            cat('Height (m): ', object@H, '\n')
            cat('Generator (Wp): ', object@Pg, '\n')
            ## cat('Monthly averages:\n')
            ## print(object@prodDm)
            ## cat('\nYearly values:\n')
            ## print(object@prody)
          }
          )
