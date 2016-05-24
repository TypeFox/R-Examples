"getAutoClp0" <-
  function (model) {
    # autoclp0 is list with elements
    # oldRes: the output of the fitModel
    # ind: the index in the result to use. defaults to 1.  
    oldRes <- model@autoclp0[["oldRes"]]
    ind <- if(length(model@autoclp0$ind) < 1) 1 else
      model@autoclp0$ind
    clp <- getSpecList(oldRes$currModel, oldRes$currTheta)[[ind]]
    negclp <- which(clp<0, arr.ind=TRUE)
    clpval <- slot(oldRes$currModel@modellist[[ind]],
                   oldRes$currModel@modellist[[ind]]@clpType)
    getConstr <- function(r) list(high = clpval[ r[1] ], low = clpval[ r[1] ],
                                  comp = r[2])
    oldclp0 <- oldRes$currModel@modellist[[ind]]@clp0
    clp0 <- vector("list", nrow(negclp))
    for(i in 1:nrow(negclp)) {
      clp0[[i]] <- getConstr(negclp[i,])
    }
    clpall <- append(oldclp0, clp0)
    clpall
  }
