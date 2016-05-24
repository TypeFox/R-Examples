setGeneric("plot") 
setMethod(f = "plot", 
          signature = c(x=".MoveTrackSingle", y="missing"), 
          function(x, y,asp=1, ...){
            plot(coordinates(x), asp=asp,...)
          })

setMethod(f = "plot", 
          signature = c(x=".MoveTrackStack", y="missing"), 
          function(x, y, type='p',asp=1, ...){
            plot(coordinates(x), type='n',asp=asp, ...)
            if(type %in% c('p','o','b'))
              res <- points(x,...)
            if(type %in% c('l','o','b'))
              res <- lines(x,...)        
          })

setMethod(f = "plot", 
          signature = c(x=".MoveTrackSingleBurst", y="missing"), 
          function(x, y, type='p',asp=1, ...){
            plot(coordinates(x),asp=asp, type="n",...)
            if(type %in% c("p", "o", "b"))
              res <- points(x, ...)
            if(type %in% c("l", "o", "b"))
              res <- lines(x, ...)
          }) 

