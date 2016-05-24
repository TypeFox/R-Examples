setGeneric('as.data.frame')

###as.data.frameI
setGeneric('as.data.frameI',
           function(object, complete=FALSE, day=FALSE){standardGeneric('as.data.frameI')})

setMethod('as.data.frameI',
          signature=(object='Sol'),
          definition=function(object, complete=FALSE, day=FALSE){
            zoo0=as.zooI(object, complete=complete, day=day)
            data0=as.data.frame(zoo0)
            ind=index(zoo0)
            data0$day=doy(ind)##Incorporo dia, mes y año como columnas del data.frame
            data0$month=month(ind)
            data0$year=year(ind)
            return(data0)
          }
          )

###as.data.frameD
setGeneric('as.data.frameD', function(object, complete=FALSE){standardGeneric('as.data.frameD')})

setMethod('as.data.frameD',
          signature=(object='Sol'),
          definition=function(object, complete=FALSE){
            zoo0=as.zooD(object, complete=complete)
            data0=as.data.frame(zoo0)
            ind=index(zoo0)
            data0$day=doy(ind)##Incorporo dia, mes y año como columnas del data.frame
            data0$month=month(ind)
            data0$year=year(ind)
            return(data0)
          }
          )

###as.data.frameM
setGeneric('as.data.frameM', function(object, complete=FALSE){standardGeneric('as.data.frameM')})

setMethod('as.data.frameM',
          signature=(object='G0'),
          definition=function(object, complete=FALSE){
            zoo0=as.zooM(object, complete=complete)
            data0=as.data.frame(zoo0)
            ind=index(zoo0)
            data0$month=month(ind)
            data0$year=year(ind)
            return(data0)
          }
          )

###as.data.frameY
setGeneric('as.data.frameY', function(object, complete=FALSE){standardGeneric('as.data.frameY')})

setMethod('as.data.frameY',
          signature=(object='G0'),
          definition=function(object, complete=FALSE){
            zoo0=as.zooY(object, complete=complete)
            data0=as.data.frame(zoo0)
            ind=index(zoo0)
            data0$year=ind
            return(data0)
          }
          )
