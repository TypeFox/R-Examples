setMethod("as.hclust", signature("AggExResult"),
          function(x, base=0.05)
          {
              if (x@maxNoClusters < 2)
                  stop("cannot create 'hclust' object with less than 2 objects")

              if (base < 0 || base >= 1)
                  stop("'base' must be at least 0 and smaller than 1")

              mini <- min(x@height)
              maxi <- max(x@height)
              auxH <- x@height <- base + (1 - base) * (-x@height + maxi) /
                                                      (maxi - mini)

              to <- list(merge=x@merge, height=auxH, labels=x@labels,
                         order=x@order)
              class(to) <- "hclust"

              to
          })

setMethod("as.hclust", signature("ExClust"),
          function(x, base=0.05, ...)
          {
              if (all(dim(x@sim) <= 1))
                  stop("similarity matrix not included in object")

              as.hclust(aggExCluster(x@sim, x, ...))
          })

setMethod("as.dendrogram", signature("AggExResult"),
          function(object, base=0.05, useNames=TRUE)
          {
              if (object@maxNoClusters < 2)
                  stop("cannot create 'hclust' object with less than 2 objects")

              if (base < 0 || base >= 1)
                  stop("'base' must be at least 0 and smaller than 1")

              obj <- as.hclust(object, base=base)

              z <- list()

              oHgt  <- obj$height
              hMax <- oHgt[object@maxNoClusters]

              topLevel <- object@clusters[[object@maxNoClusters]]

              if (length(names(object@exemplars)) == 0 || !useNames)
                  topLevel <- lapply(object@clusters[[object@maxNoClusters]],
                                     as.character)
              else
                  topLevel <- lapply(object@clusters[[object@maxNoClusters]],
                                     names)

              for (k in 1:length(obj$height))
              {
                  x <- obj$merge[k, ]

                  if (x[1] < 0)
                  {
                      if (length(topLevel[[-x[1]]]) == 1)
                      {
                          leftDend <- topLevel[[-x[1]]]
                          attr(leftDend, "label") <- topLevel[[-x[1]]]
                          attr(leftDend, "members") <- 1
                          attr(leftDend, "midpoint") <- 0
                          attr(leftDend, "height") <- 0
                          attr(leftDend, "leaf") <- TRUE
                      }
                      else
                      {
                          leftDend <- lapply(topLevel[[-x[1]]],
                                             function(elem)
                                             {
                                                 attr(elem, "label") <- elem
                                                 attr(elem, "members") <- 1
                                                 attr(elem, "height") <- 0
                                                 attr(elem, "leaf") <- TRUE
                                                 elem
                                             })
                          attr(leftDend, "members") <- length(topLevel[[-x[1]]])
                          attr(leftDend, "height") <- base / 2
                          attr(leftDend, "midpoint") <-
                              (length(topLevel[[-x[1]]]) - 1) / 2
                      }
                  }
                  else
                      leftDend <- z[[as.character(x[1])]]

                  if (x[2] < 0)
                  {
                      if (length(topLevel[[-x[2]]]) == 1)
                      {
                          rightDend <- topLevel[[-x[2]]]
                          attr(rightDend, "label") <- topLevel[[-x[2]]]
                          attr(rightDend, "members") <- 1
                          attr(rightDend, "midpoint") <- 0
                          attr(rightDend, "height") <- 0
                          attr(rightDend, "leaf") <- TRUE
                      }
                      else
                      {
                          rightDend <- lapply(topLevel[[-x[2]]],
                                              function(elem)
                                              {
                                                  attr(elem, "label") <- elem
                                                  attr(elem, "members") <- 1
                                                  attr(elem, "height") <- 0
                                                  attr(elem, "leaf") <- TRUE
                                                  elem
                                              })
                          attr(rightDend, "members") <-
                              length(topLevel[[-x[2]]])
                          attr(rightDend, "height") <- base / 2
                          attr(rightDend, "midpoint") <-
                              (length(topLevel[[-x[2]]]) - 1) / 2
                      }
                  }
                  else
                      rightDend <- z[[as.character(x[2])]]

                  zk <- list(leftDend, rightDend)
                  attr(zk, "height") <- obj$height[k]
                  attr(zk, "members") <- attr(leftDend, "members") +
                                         attr(rightDend, "members")
                  attr(zk, "midpoint") <- (attr(leftDend, "members") +
                                           attr(leftDend, "midpoint") +
                                           attr(rightDend, "midpoint")) / 2

                  z[[as.character(k)]] <- zk
              }

              z <- z[[as.character(k)]]
              class(z) <- "dendrogram"

              z
          })

setMethod("as.dendrogram", signature("ExClust"),
          function(object, base=0.05, useNames=TRUE, ...)
          {
              if (all(dim(object@sim) <= 1))
                  stop("similarity matrix not included in object")

              as.dendrogram(aggExCluster(object@sim, object, ...), base=base,
                            useNames=useNames)
          })
