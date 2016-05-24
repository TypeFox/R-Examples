###summary

setMethod("summary","Evaluation",
          function(object, runs0 = 1:nrow(result(object)),
                   dims0 = 1:ncol(result(object)),
                   inList = FALSE, ...){
            cat(gettextf("name of Evaluation: %s\n", name(object)))
            if(!inList)
               {cat(gettextf("name of Dataobject: %s\n", name(Data(object))))
                cat(gettextf("name of Datafile: %s\n", filename(object)))}
            cat(gettextf("estimator: %s\n",
                as.character(call.ev(object)$estimator)))
            cat(gettext("Result:\n"));
            print(summary(result(object)[runs0,dims0]))
                  ### still do not see why necessary...
          }
          )


setMethod("summary","EvaluationList",
          function(object, runs0 = 1:nrow(result(Elist(object)[[1]])),
                   dims0 = 1:ncol(result(Elist(object)[[1]])),
                   evals0 = 1:length(Elist(object)),...){
            levals0 <- min(length(evals0),
                        getdistrTEstOption("MaxNumberofSummarizedEvaluations"))
            if(levals0<length(evals0))
               warning(gettextf(
  "Your evaluation list is too big; only %i% evaluations are printed",
                        levels ))
            cat(gettextf("name of Evaluation List: %s\n",name(object)))
            if(!is.null(name(Elist(object)[[1]])))
                cat(gettextf("name of Dataobject: %s\n",
                              name(Elist(object)[[1]])))
            if(!is.null(filename(Elist(object)[[1]])))
                cat(gettextf("name of Datafile: %s\n",
                              filename(Elist(object)[[1]])))
            len <- length(Elist(object)[evals0[1:levals0]])
            for(i in 1:len)
                {cat("----------------------------------\n")
                 summary(Elist(object)[[evals0[i]]], runs0 = runs0,
                         dims0 = dims0, inList = TRUE)}
          }
          )
