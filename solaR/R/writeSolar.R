setGeneric('writeSolar', function(object, file,
                                  complete=FALSE, day=FALSE,
                                  timeScales=c('i', 'd', 'm', 'y'), sep=',',
                                  ...){
    standardGeneric('writeSolar')})

setMethod('writeSolar', signature=(object='Sol'),
          definition=function(object, file, complete=FALSE, day=FALSE,
          timeScales=c('i', 'd', 'm', 'y'), sep=',', ...){
              name <- strsplit(file, '\\.')[[1]][1]
              ext <- strsplit(file, '\\.')[[1]][2]
              timeScales <- match.arg(timeScales, several.ok=TRUE)
              if ('i' %in% timeScales) {
                  zI <- as.zooI(object, complete=complete, day=day)
                  write.zoo(zI, file=file, sep=sep, ...)
              }
              if ('d' %in% timeScales) {
                  zD <- as.zooD(object, complete=complete)
                  write.zoo(zD,
                            file=paste(name, 'D', ext, sep='.'),
                            sep=sep, ...)
              }

              if ('m' %in% timeScales) {
                  zM <- as.zooM(object, complete=complete)
                  write.zoo(zM,
                            file=paste(name, 'M', ext, sep='.'),
                            sep=sep, ...)
              }

              if ('y' %in% timeScales) {
                  zY <- as.zooY(object, complete=complete)
                  write.zoo(zY,
                            file=paste(name, 'Y', ext, sep='.'),
                            sep=sep, ...)
              }

          })
