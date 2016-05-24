setMethod("sketchGrid",
          signature(                    
                    layers = "ANY",                                        
                    shapeCross = "ANY"),
          function(i,
                   j,
                   rowLimit,
                   colLimit,
                   layers,
                   shapeCross,
                   excludeCenter,
                   ...)
          {
### used for designed grid
            coordCross <- extendHorVer(i = i,
                                       j = j,
                                       shapeCross = shapeCross,
                                       rowLimit = rowLimit,
                                       colLimit = colLimit)
            ##
            coordCirc <- circularExtension(i = i,
                                           j = j,
                                           layers = layers,
                                           rowLimit = rowLimit,
                                           colLimit = colLimit)
            ##
            if(identical(excludeCenter,FALSE))
              {
                coordGrid <- rbind(cbind(i,j),
                                   coordCross,
                                   coordCirc)
              }
            ##
            else
              {
                coordGrid <- rbind(coordCross,
                                   coordCirc)
              }
            ##
            plot.new()
            plot.window(xlim = c(0,colLimit),
                        ylim = c(0,rowLimit))
            ##
            axis(2,at=1:rowLimit,
                 labels=rowLimit:1,
                 cex.axis=0.5)
            ##
            axis(1,at=1:colLimit,
                 labels=1:colLimit,
                 cex.axis=0.5)
            ##
            ##
            mtext("column",
                  1,
                  line = 3,
                  font = 2)
            ##
            mtext("row",
                  2,
                  line = 3,
                  font = 2)
            ##
            ## empty cells
            ##
            points(x = rep(1:colLimit,rowLimit),
                   y = (rowLimit+1)- rep(1:rowLimit, each = colLimit))
            ##
            ## color cells that are included in the grid
            ##
            points(x = coordGrid[,2],
                   y = (rowLimit + 1)-coordGrid[,1],
                   pch = 21, bg = "black")
            ##
          } ## end if definition
          ) 




