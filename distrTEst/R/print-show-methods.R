## print Methode
setMethod("print","Evaluation",
          function(x, runs0 = 1:nrow(result(x)), dims0 = 1:ncol(result(x)),
                   inList = FALSE, ...){
            cat(gettextf("An Evaluation Object\n"))
            if(!inList)
               {if(!is.null(name(x)))
                    cat(gettextf("name of Dataobject: %s\n",name(x)))
                if(!is.null(filename(x)))
                    cat(gettextf("name of Datafile: %s\n",filename(x)))
               }
            if(!is.null(call.ev(x)))
                cat(gettextf("estimator: %s\n",
                    as.character(call.ev(x)$estimator)))
            if(!is.null(result(x))) {cat(gettextf("Result: "))
                                     cat(str(result(x)[runs0,dims0]))}
          })

setMethod("print","EvaluationList",
          function(x,
                   runs0 =  1:nrow(result(Elist(x)[[1]])),
                   dims0 =  1:ncol(result(Elist(x)[[1]])),
                   evals0 = 1:length(Elist(x)),...){

            levals0 <- min(getdistrTEstOption("MaxNumberofPrintedEvaluations"),
                           length(evals0))

            if(levals0<length(evals0))
                warning(paste("your evaluation list is too big; only ", levals0,
                              "evaluations are printed"))

            cat(gettextf("An EvaluationList Object\n"))
            if(!is.null(name(x)))
                cat(gettextf("name of Evaluation List: %s\n",name(x)))
            if(!is.null(name(Elist(x)[[1]])))
                cat(gettextf("name of Dataobject: %s\n",name(Elist(x)[[1]])))
            if(!is.null(filename(Elist(x)[[1]])))
                cat(gettextf("name of Datafile: %s\n",filename(Elist(x)[[1]])))

            len <- length(Elist(x)[evals0[1:levals0]])
            for(i in 1:len)
                {cat("----------------------------------\n")
                 print(Elist(x)[[evals0[i]]], runs0 = runs0, dims0 = dims0,
                       inList = TRUE)}
          }
          )

setMethod("show","Evaluation", function(object)print(object))
setMethod("show","EvaluationList", function(object)print(object))
