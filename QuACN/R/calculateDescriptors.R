calculateDescriptors <- function(graphs, ..., labels=FALSE, log=FALSE) {
  argv <- list(...)

  if (class(graphs) != "list") {
    graphs <- list(graphs)
  }

  result <- lapply(1:length(graphs), function(n) {
    g <- graphs[[n]]
    if (log)
      cat("Calculating selected descriptors for network",
          n, "of", length(graphs), "\n")

    values <- list()
    cache <- new.env()
    i <- 1

    while (i <= length(argv)) {
      funcs <- argv[[i]]
      if (is.numeric(funcs))
        funcs <- .getFunctionsByNumber(funcs)
      i <- i + 1

      if (i <= length(argv) && class(argv[[i]]) == "list") {
        params <- argv[[i]]
        i <- i + 1
      }
      else
        params <- list()

      for (func in funcs) {
        value <- lapply(.callFunction(func, g, params, cache), as.numeric)
        valueName <- sub("^\\.", "", func)
        valueName <- sub("^.*@", "", valueName)
        if (length(value) == 1)
          values[[valueName]] <- unlist(value)
        else {
          for (name in names(value))
            values[[paste(valueName, name, sep=".")]] <- value[[name]]
        }
      }
    }

    values
  })

  cn <- names(result[[1]])
  result <- matrix(unlist(result), byrow=TRUE, ncol=length(cn))
  colnames(result) <- cn
  result <- data.frame(result)

  if (labels)
    colnames(result) <- sapply(colnames(result), getLabels)

  result
}

.callFunction <- function(func, g, extra, cache) {
  if (!(func %in% ls(cache, all.names = TRUE))) {
    # determine dependencies of the function
    funcname <- strsplit(func, split = "@")[[1]][[1]]
    f <- match.fun(funcname)
    params <- as.list(formals(f))
    for (key in names(params)) {
      if (key == "g")
        params[[key]] <- g
      else if (key == "am")
        params[[key]] <- .callFunction("adjacencyMatrix", g, list(), cache)
      else if (key == "dist")
        params[[key]] <- .callFunction("distanceMatrix", g, list(), cache)
      else if (key == "dsc")
        params[[key]] <- .callFunction("distSumConnectMatrix", g, list(), cache)
      else if (key == "deg")
        params[[key]] <- .callFunction(".degreeHelper", g, list(), cache)
      else if (key == "wien")
        params[[key]] <- .callFunction("wiener", g, list(), cache)
      else if (key == "loc")
        params[[key]] <- .callFunction("localClusteringCoeff", g, list(), cache)
      else if (key == "ita")
        params[[key]] <- .callFunction("totalAdjacency", g, list(), cache)
      else if (key == "one.eds")
        params[[key]] <- .callFunction("edgeDeletedSubgraphs@1", g, list(`gs` = list(g)), cache)
      else if (key == "two.eds") {
        one.eds <- .callFunction("edgeDeletedSubgraphs@1", g, list(`gs` = list(g)), cache)
        params[[key]] <- .callFunction("edgeDeletedSubgraphs@2", g, list(`gs` = one.eds), cache)
      }
      else if (key %in% names(extra))
        params[[key]] <- extra[[key]]
    }

    assign(func, do.call(f, params), cache)
  }

  get(func, cache)
}

.getFunctionsByNumber <- function(numbers) {
  sapply(numbers, function(number) {
    if (number %% 1000 == 0)
      .functions[[number %/% 1000]]
    else
      .functions[[number %/% 1000]][[number %% 1000]]
  })
}

.degreeHelper <- function(g) graph::degree(g)

.dobrynin <- function(g, dist=NULL) {
  dobrynin(g, dist)[c("eccentricityGraph", "avgeccOfG", "ecentricGraph",
                      "graphIntegration", "unipolarity", "variation",
                      "centralization", "avgDistance",
                      "meanDistVertexDeviation")]
}

.topologicalInfoContent <- function(g, dist=NULL, deg=NULL) {
  topologicalInfoContent(g, dist, deg)[["entropy"]]
}

.metaInfoTheoreticGCM <- function(infofunct, coeff) {
  name <- paste(".infoTheoreticGCM", infofunct, coeff, sep="_")
  assign(name,
    function(g, dist=NULL, lambda=1000, alpha=0.5, prec=53) {
      result <- infoTheoreticGCM(g, dist=dist, coeff=coeff, infofunct=infofunct,
        lambda=lambda, alpha=alpha, prec=prec)
      result[c("entropy", "distance")]
  }, parent.env(environment()))
  name
}

.metaInfoTheoretic <- function(funct, coeff) {
  if (is.null(coeff))
    name <- paste(".", funct, sep="")
  else
    name <- paste(paste(".", funct, sep=""), coeff, sep="_")
  assign(name,
    function(g, lambda=1000) {
      args <- list(`g` = g, `lambda` = lambda)
      if (!is.null(coeff))
        args[["coeff"]] <- coeff
      result <- do.call(funct, args)
      result[c("entropy", "distance")]
  }, parent.env(environment()))
  name
}

.metaEigenvalueBased <- function(matrix_function) {
  name <- paste(".eigenvalueBased", matrix_function, sep="_")
  assign(name,
    function(g, s=1) eigenvalueBased(g, matrix_function=matrix_function, s=s),
    parent.env(environment()))
  name_2 <- paste(name, "2", sep="_")
  assign(name_2,
    function(g) eigenvalueBased(g, matrix_function=matrix_function, s=2),
    parent.env(environment()))
  c(name, name_2)
}

.functions <- list(
  # group 1000
  c(
    "wiener",                                             # 1001
    "harary",                                             # 1002
    "balabanJ",                                           # 1003
    "meanDistanceDeviation",                              # 1004
    "compactness",                                        # 1005
    "productOfRowSums",                                   # 1006
    "hyperDistancePathIndex",                             # 1007
    ".dobrynin"                                           # 1008
  ),
  # group 2000
  c(
    "totalAdjacency",                                     # 2001
    "zagreb1",                                            # 2002
    "zagreb2",                                            # 2003
    "modifiedZagreb",                                     # 2004
    "augmentedZagreb",                                    # 2005
    "variableZagreb",                                     # 2006
    "randic",                                             # 2007
    "complexityIndexB",                                   # 2008
    "normalizedEdgeComplexity",                           # 2009
    "atomBondConnectivity",                               # 2010
    "geometricArithmetic1",                               # 2011
    "geometricArithmetic2",                               # 2012
    "geometricArithmetic3",                               # 2013
    "narumiKatayama"                                      # 2014
  ),
  # group 3000
  c(
    ".topologicalInfoContent",                            # 3001
    "bonchev1",                                           # 3002
    "bonchev2",                                           # 3003
    "bertz",                                              # 3004
    "radialCentric",                                      # 3005
    "vertexDegree",                                       # 3006
    "balabanlike1",                                       # 3007
    "balabanlike2",                                       # 3008
    "graphVertexComplexity",                              # 3009
    "informationBondIndex",                               # 3010
    "edgeEqualityMIC",                                    # 3011
    "edgeMagnitudeMIC",                                   # 3012
    "symmetryIndex",                                      # 3013
    "bonchev3",                                           # 3014
    "graphDistanceComplexity",                            # 3015
    "distanceDegreeMIC",                                  # 3016
    "distanceDegreeEquality",                             # 3017
    "distanceDegreeCompactness",                          # 3018
    "informationLayerIndex"                               # 3019
  ),
  # group 4000
  c(
    "mediumArticulation",                                 # 4001
    "efficiency",                                         # 4002
    "graphIndexComplexity",                               # 4003
    "offdiagonal",                                        # 4004
    "spanningTreeSensitivity",                            # 4005
    "distanceDegreeCentric",                              # 4006
    "distanceCodeCentric"                                 # 4007
  ),
  # group 5000
  c(
    .metaInfoTheoreticGCM("vertcent", "exp"),             # 5001
    .metaInfoTheoreticGCM("vertcent", "lin"),             # 5002
    .metaInfoTheoreticGCM("sphere", "exp"),               # 5003
    .metaInfoTheoreticGCM("sphere", "lin"),               # 5004
    .metaInfoTheoreticGCM("pathlength", "exp"),           # 5005
    .metaInfoTheoreticGCM("pathlength", "lin"),           # 5006
    .metaInfoTheoreticGCM("degree", "exp"),               # 5007
    .metaInfoTheoreticGCM("degree", "lin"),               # 5008
    .metaInfoTheoretic("infoTheoreticLabeledV1", "exp"),  # 5009
    .metaInfoTheoretic("infoTheoreticLabeledV1", "lin"),  # 5010
    .metaInfoTheoretic("infoTheoreticLabeledV2", NULL),   # 5011
    .metaInfoTheoretic("infoTheoreticLabeledE", "exp"),   # 5012
    .metaInfoTheoretic("infoTheoreticLabeledE", "lin")    # 5013
  ),
  # group 6000
  c(
    .metaEigenvalueBased("adjacencyMatrix"),              # 6001, 6002
    .metaEigenvalueBased("laplaceMatrix"),                # 6003, 6004
    .metaEigenvalueBased("distanceMatrix"),               # 6005, 6006
    .metaEigenvalueBased("distancePathMatrix"),           # 6007, 6008
    .metaEigenvalueBased("augmentedMatrix"),              # 6009, 6010
    .metaEigenvalueBased("extendedAdjacencyMatrix"),      # 6011, 6012
    .metaEigenvalueBased("vertConnectMatrix"),            # 6013, 6014
    .metaEigenvalueBased("randomWalkMatrix"),             # 6015, 6016
    .metaEigenvalueBased("weightStrucFuncMatrix_lin"),    # 6017, 6018
    .metaEigenvalueBased("weightStrucFuncMatrix_exp"),    # 6019, 6020
    "energy",                                             # 6021
    "laplacianEnergy",                                    # 6022
    "estrada",                                            # 6023
    "laplacianEstrada",                                   # 6024
    "spectralRadius"                                      # 6025
  ),
  # group 7000
  c(
    "oneEdgeDeletedSubgraphComplexity",                   # 7001
    "twoEdgesDeletedSubgraphComplexity",                  # 7002
    "globalClusteringCoeff"                               # 7003
    # "localClusteringCoeff"                              # 7004 TODO
  ),
  # group 8000
  c(
    "connectivityID",                                     # 8001
    "minConnectivityID",                                  # 8002
    "primeID",                                            # 8003
    "bondOrderID",                                        # 8004
    "balabanID",                                          # 8005
    "minBalabanID",                                       # 8006
    "weightedID",                                         # 8007
    "huXuID"                                              # 8008
  )
)
