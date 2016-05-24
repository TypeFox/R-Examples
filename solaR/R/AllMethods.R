##################################################################
## GET data
##################################################################

setGeneric('getData', function(object){standardGeneric('getData')})

setMethod('getData',##Solo definido para Meteo, de forma que siempre devuelve valores de partida
          signature=(object='Meteo'),
          definition=function(object){
            result=object@data
            return(result)
          }
          )

setGeneric('getG0', function(object){standardGeneric('getG0')})

setMethod('getG0',##Solo definido para Meteo, de forma que siempre devuelve valores de partida
          signature=(object='Meteo'),
          definition=function(object){
            result=getData(object)
            return(result$G0)
          }
          )

###Latitud
setGeneric('getLat', function(object, units='rad'){standardGeneric('getLat')})

setMethod('getLat',
          signature=(object='Sol'),
          definition=function(object, units='rad'){
            stopifnot(units %in% c('deg', 'rad'))
            res=switch(units,
              rad=d2r(object@lat),
              deg=object@lat)
            return(res)
          }
          )

setMethod('getLat',
          signature=(object='Meteo'),
          definition=function(object, units='rad'){
            stopifnot(units %in% c('deg', 'rad'))
            res=switch(units,
              rad=d2r(object@latData),
              deg=object@latData)
            return(res)
          }
          )
setMethod('getLat',
          signature=(object='G0'),
          definition=function(object, units='rad'){
            getLat(as(object, 'Sol'), units=units)
          }
          )

##################################################################
## INDEX
##################################################################


setGeneric('indexD', function(object){standardGeneric('indexD')})
setMethod('indexD',
          signature=(object='Meteo'),
          definition=function(object){
            return(index(object@data))
          }
          )

setMethod('indexD',
          signature=(object='Sol'),
          definition=function(object){
            return(index(object@solD))
          }
          )

setMethod('indexD',
          signature=(object='G0'),
          definition=function(object){
            indexD(as(object, 'Sol'))
          }
          )


setGeneric('indexI', function(object){standardGeneric('indexI')})
setMethod('indexI',
          signature=(object='Sol'),
          definition=function(object){
            return(index(object@solI))
          }
          )

setGeneric('indexRep', function(object){standardGeneric('indexRep')})
setMethod('indexRep',
          signature=(object='Sol'),
          definition=function(object){
            return(object@match)
          }
          )




