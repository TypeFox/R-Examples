setMethod("showobj", "SphericalDistribution",
          function(object, className = class(object)[1]){
            txt <- gettextf("Distribution Object of Class: %s\n", className)
            parameter = param(object)
            Names = slotNames(parameter)
            if(length(Names) > 1){
              for(i in Names[Names != "name"]){
                if(is.matrix(slot(parameter, i))){
                   txt <- c(txt, gettextf("%s:\n", i))
                   txt <- c(txt, "        ")
                   for(k in 1:ncol(slot(parameter, i)))
                       txt <- c(txt, gettextf("[,%0d]   ", k))
                   txt <- c(txt,"\n")
                   for(j in 1:nrow(slot(parameter, i))){
                       txt <- c(txt, gettextf("[%0d,]  ", j))
                       for(k in 1:ncol(slot(parameter, i)))
                          txt <- c(txt, gettextf("% 2.2f, ",
                                   slot(parameter, i)[j,k]))
                       txt <- c(txt, "\n")
                   }
                }else{ txt0 <- if(length(slot(parameter,i))>1)
                                 paste("(",paste(slot(parameter, i),
                                         collapse=","), ")",sep="")
                               else
                                 paste(slot(parameter, i))
                 txt <- c(txt,
                          gettextf("%s: %s\n", i, txt0))
            }}}
            txt <- c(txt, "\n Distribution of Lengths:\n",
                          showobj(radDistr(object)))
            return(txt)
          })

setMethod("show", "SphericalDistribution",
          function(object){
            cls <- class(object)[1]
            cat(showobj(object, className = cls))
            ws <- distr:::.IssueWarn(object@.withArith, object@.withSim)
            if(!is.null(ws$msgA)) warning(ws$msgA)
            if(!is.null(ws$msgS)) warning(ws$msgS)
            }
          )
