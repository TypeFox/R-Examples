## summary gives info about either a SOM unit or an object

summary.wccsom <- function(object,
                           type=c("unit", "object", "smoothness", "quality"),
                           nr, labels, data = object$data,
                           classif = object$unit.classif,
                           wccs = object$wccs,
                           properties = NULL, ...)
{
  type <- match.arg(type)

  codes <- object$codes

  if (type != "smoothness") { # incase output of wccassign is given
    if (is.list(classif) &&
        !is.null(classif$classif) &&
        !is.null(classif$wccs)) {
      classif <- classif$classif
      wccs <- classif$wccs
    }
  }
  
  switch(type,
         unit = {
           objs <- which(classif == nr)

           if (!is.null(properties) &
               (is.data.frame(properties) | is.matrix(properties))) {
             if (missing(labels)) labels <- row.names(properties)[objs]
             else labels <- labels[objs]
             
             summtable <- as.data.frame(cbind(wccs = rep(0, length(labels)),
                                              properties[objs,]))
             row.names(summtable) <- labels
           } else {
             if (missing(labels)) labels <- objs
             else labels <- labels[objs]

             summtable <- data.frame(wccs = rep(0, length(labels)),
                                     row.names = labels)
           }

           if (is.null(object$wccs)) {
             for (i in seq(along=objs))
               summtable[i,1] <- wcc(data[objs[i],], codes[nr,],
                                     object$trwdth)
           } else {
             summtable[,1] <- wccs[objs]
           }
           
           cat("Agreement with codebook vector of ", length(objs),
               " objects mapped to cell ",
               nr,":\n", sep="")
           print(summtable)

           invisible(summtable)
         },
         object = {
           objs <- which(classif == classif[nr])
           objs <- objs[objs != nr]
           
           if (!is.null(properties) &
               (is.data.frame(properties) | is.matrix(properties))) {
             if (missing(labels)) {
               target <- row.names(properties)[nr]
               labels2 <- row.names(properties)[objs]
             } else {
               target <- labels[nr]
               labels2 <- labels[objs]
             }
             
             summtable <- as.data.frame(cbind(wccs = rep(0, length(labels2)),
                                              properties[objs,]))
             row.names(summtable) <- labels2
           } else {
             if (missing(labels)) {
               labels2 <- objs
               target <- nr
             } else {
               labels2 <- labels[objs]
               target <- labels[nr]
             }
             
             summtable <- data.frame(wccs = rep(0, length(labels2)),
                                     row.names = labels2)
           }

           for (i in seq(along=objs))
             summtable[i,1] <- wcc(data[objs[i],], data[nr,],
                                   object$trwdth)
           
           cat("Agreement of object ", target,
               " with ", length(objs),
               " other objects mapped to cell ", classif[nr],
               ":\n", sep="")
           print(summtable)

           invisible(summtable)
         },
         smoothness = {
           nunits <- nrow(object$codes)
           wcctab <- matrix(0, nunits, nunits)
           wghts <- 1 - (0:object$trwdth)/object$trwdth
           for (i in 1:(nunits-1)) {
             for (j in (i+1):nunits) {
               wcctab[i,j] <- wcc(object$codes[i,],
                                  object$codes[j,],
                                  object$trwdth,
                                  wghts = wghts,
                                  acors = object$acors[c(i,j)])
             }
           }
           wcctab <- wcctab + t(wcctab)

           nhbrdist <- as.matrix(unit.distances(object$grid, object$toroidal))
           diag(nhbrdist) <- NA
           
           smthmat <- matrix(0, nunits, 2)
           dimnames(smthmat) <- list(NULL,
                                     c("neighbours", "non-neighbours"))
           
           thresh <- ifelse(object$grid$topo == "hexagonal", 1, sqrt(2))
           for (i in 1:nunits) {
             faraway <- which(nhbrdist[i,]-thresh > 1e-5)
             smthmat[i,1] <- mean(wcctab[i,-faraway])
             smthmat[i,2] <- mean(wcctab[i,faraway])
           }

           cat("Mean neighbour WCC: ", formatC(mean(smthmat[,1]), digits=3),
               ", mean other WCC: ", formatC(mean(smthmat[,2]), digits=3),
               "\nRatio: ",
               formatC(mean(smthmat[,1]) / mean(smthmat[,2]), digits=3),
               "\n", sep="")

           invisible(smthmat)
         },
         quality = {
           within <- nwithin <- between <- nbetween <- 0
           nobj <- nrow(data)
           wghts <- 1 - (0:object$trwdth)/object$trwdth
           if (!is.null(object$data.acors)) {
             data.acors <- object$data.acors
           } else {
             data.acors <- wacmat(data, object$trwdth, wghts)
           }
           
           for (i in 1:(nobj-1))
             for (j in (i+1):nobj) {
               huhn <- wcc(data[i,], data[j,],
                           object$trwdth, wghts=wghts,
                           acors=data.acors[c(i,j)])
               if (classif[i] == classif[j]) {
                 nwithin <- nwithin+1
                 within <- within + huhn
               } else {
                 nbetween <- nbetween+1
                 between <- between + huhn
               }
             }

           huhn <- within*nbetween/(nwithin*between)
           cat("WCC enrichment for neighbouring units:",
               formatC(huhn, digits=4), "\n")
           
           invisible(huhn)
         }
         )
}

