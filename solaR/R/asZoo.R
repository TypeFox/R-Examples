###as.zooM
setGeneric('as.zooM', function(object, complete=FALSE){standardGeneric('as.zooM')})

setMethod('as.zooM',
          signature=(object='G0'),
          definition=function(object, complete=FALSE){
            return(object@G0dm)
          }
          )

setMethod('as.zooM',
          signature=(object='Gef'),
          definition=function(object, complete=FALSE){
            res0 <- object@Gefdm
            if (complete) {
              res1 <- as.zooM(as(object, 'G0'))
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

setMethod('as.zooM',
          signature=(object='ProdGCPV'),
          definition=function(object, complete=FALSE){
            res0 <- object@prodDm
            if (complete) {
              res1 <- as.zooM(as(object, 'Gef'), complete=TRUE)
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

setMethod('as.zooM',
          signature=(object='ProdPVPS'),
          definition=function(object, complete=FALSE){
            res0 <- object@prodDm
            if (complete) {
              res1 <- as.zooM(as(object, 'Gef'), complete=TRUE)
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

###as.zooY
setGeneric('as.zooY', function(object, complete=FALSE){standardGeneric('as.zooY')})

setMethod('as.zooY',
          signature=(object='G0'),
          definition=function(object, complete=FALSE){
            return(object@G0y)
          }
          )

setMethod('as.zooY',
          signature=(object='Gef'),
          definition=function(object, complete=FALSE){
            res0 <- object@Gefy
            if (complete) {
              res1 <- as.zooY(as(object, 'G0'))
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

setMethod('as.zooY',
          signature=(object='ProdGCPV'),
          definition=function(object, complete=FALSE){
            res0 <- object@prody
            if (complete) {
              res1 <- as.zooY(as(object, 'Gef'), complete=TRUE)
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

setMethod('as.zooY',
          signature=(object='ProdPVPS'),
          definition=function(object, complete=FALSE){
            res0 <- object@prody
            if (complete) {
              res1 <- as.zooY(as(object, 'Gef'), complete=TRUE)
              return(CBIND(res1, res0))
            } else {
              return(res0)
            }
          }
          )

###as.zooD
setGeneric('as.zooD', function(object, complete=FALSE){standardGeneric('as.zooD')})

setMethod('as.zooD',
          signature=(object='Sol'),
          definition=function(object, complete=FALSE){#complete esta por compatibilidad con los otros metodos
            res <- object@solD
            return(res)
          }
          )

setMethod('as.zooD',
          signature=(object='G0'),
          definition=function(object, complete=FALSE){
            res1 <- as.zooD(as(object, 'Sol'))
            res2 <- object@G0D
            if (complete) {
              res1=coredata(res1)
              res2=coredata(res2)
              return(zoo(cbind(res1, res2), indexD(object)))
            } else {
              return(res2[,c('G0d', 'D0d', 'B0d')])}
          }
          )

setMethod('as.zooD',
          signature=(object='Gef'),
          definition=function(object, complete=FALSE){
            res1 <- as.zooD(as(object, 'G0'), complete=TRUE)
            res2 <- object@GefD
            if (complete) {
              res1=coredata(res1)
              res2=coredata(res2)
              return(zoo(cbind(res1, res2), indexD(object)))
            } else {
              return(res2[,c('Gefd', 'Defd', 'Befd')])
            }
          }
          )


setMethod('as.zooD',
          signature=(object='ProdGCPV'),
          definition=function(object, complete=FALSE){
            res1 <- as.zooD(as(object, 'Gef'), complete=TRUE)
            res2 <- object@prodD
            if (complete) {
              res1=coredata(res1)
              res2=coredata(res2)
              return(zoo(cbind(res1, res2), indexD(object)))
            } else {
              return(res2[,c('Eac', 'Edc', 'Yf')])
            }
          }
          )

setMethod('as.zooD',
          signature=(object='ProdPVPS'),
          definition=function(object, complete=FALSE){
            res1 <- as.zooD(as(object, 'Gef'), complete=TRUE)
            res2 <- object@prodD
            if (complete) {
              res1=coredata(res1)
              res2=coredata(res2)
              return(zoo(cbind(res1, res2), indexD(object)))
            } else {
              return(res2[,c('Eac', 'Qd', 'Yf')])
            }
          }
          )

###as.zooI
setGeneric('as.zooI',
           function(object, complete=FALSE, day=FALSE){standardGeneric('as.zooI')})

setMethod('as.zooI',
          signature=(object='Sol'),
          definition=function(object, complete=FALSE, day=FALSE){
            res0 <- object@solI
            if (day) {
              ind <- indexRep(object)
              res2 <- coredata(object@solD)[ind,]
              res0=coredata(res0)
              return(zoo(cbind(res0, res2), indexI(object)))
            } else {return(res0)}
          }
          )

setMethod('as.zooI',
          signature=(object='G0'),
          definition=function(object, complete=FALSE, day=FALSE){
            res0 <- object@G0I
            if (complete) {
              res1 <- coredata(as.zooI(as(object, 'Sol'), day=day))
              res0=coredata(res0)
              Ta <- coredata(object@Ta)
              if (day) { ##complete&day
                ind <- indexRep(object)
                res2 <-coredata(object@G0D)[ind,]
                res <- zoo(cbind(res1, res2, res0, Ta), indexI(object))
              } else { ##complete without day
                res=zoo(cbind(res1, res0, Ta), indexI(object))
              }
              return(res)
            } else { ##neither complete nor day
              return(res0[,c('G0', 'B0', 'D0')])
            }
          }
          )

setMethod('as.zooI',
          signature=(object='Gef'),
          definition=function(object, complete=FALSE, day=FALSE){
            res0 <- object@GefI
            if (complete) {
              res1 <- coredata(as.zooI(as(object, 'G0'),
                                       complete=complete, day=day))
              res2 <- coredata(object@Theta)
              res0=coredata(res0)
              if (day) { ##complete&day
                ind <- indexRep(object)
                res3 <-coredata(object@GefD)[ind,]
                res <- zoo(cbind(res1, res2, res3, res0), indexI(object))
              } else { ##complete without day
                res=zoo(cbind(res1, res2, res0), indexI(object))
              }
              return(res)
            } else { ##neither complete nor day
              return(res0[,c('Gef', 'Bef', 'Def')])
            }
          }
          )

setMethod('as.zooI',
          signature=(object='ProdGCPV'),
          definition=function(object, complete=FALSE, day=FALSE){
            res0 <- object@prodI
            if (complete) {
              res1 <- coredata(as.zooI(as(object, 'Gef'),
                                       complete=complete, day=day))
              res0=coredata(res0)
              if (day) { ##complete&day
                ind <- indexRep(object)
                res2 <-coredata(object@prodD)[ind,]
                res <- zoo(cbind(res1, res2, res0), indexI(object))
              } else { ##complete without day
                res=zoo(cbind(res1, res0), indexI(object))
              }
              return(res)
            } else { ##neither complete nor day
              return(res0[,c('Pac', 'Pdc')])
            }
          }
          )

setMethod('as.zooI',
          signature=(object='ProdPVPS'),
          definition=function(object, complete=FALSE, day=FALSE){
            res0 <- object@prodI
            if (complete) {
              res1 <- coredata(as.zooI(as(object, 'Gef'),
                                       complete=complete, day=day))
              res0=coredata(res0)
              if (day) { ##complete&day
                ind <- indexRep(object)
                res2 <-coredata(object@prodD)[ind,]
                res <- zoo(cbind(res1, res2, res0), indexI(object))
              } else { ##complete without day
                res=zoo(cbind(res1, res0), indexI(object))
              }
              return(res)
            } else { ##neither complete nor day
              return(res0[,c('Pac', 'Q')])
            }
          }
          )
