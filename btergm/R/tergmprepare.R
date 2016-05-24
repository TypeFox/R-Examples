

# helper function for adjusting matrix dimensions and creating an offset matrix
tergmprepare <- function(formula, offset = TRUE, blockdiag = FALSE, 
    verbose = TRUE) {
  
  # extract response networks and make sure they are saved in a list
  env <- new.env()
  env$lhs.original <- deparse(formula[[2]])  # for reporting purposes later on
  env$networks <- eval(parse(text = deparse(formula[[2]])))
  if (class(env$networks) == "list" || class(env$networks) == "network.list") {
    # do nothing
  } else {
    env$networks <- list(env$networks)
  }
  
  # convert list elements to matrices if unknown data type
  for (i in 1:length(env$networks)) {
    if (!class(env$networks[[i]]) %in% c("network", "matrix", "list")) {
      tryCatch(
        {
          env$networks[[i]] <- as.matrix(env$networks[[i]])
        }, 
        error = function(cond) {
          stop(paste("Object", i, "could not be converted to a matrix."))
        }
      )
    }
  }
  
  # extract additional information
  env$num.vertices <- max(sapply(env$networks, function(x) 
      network::get.network.attribute(network::network(x), "n")))  # number nodes
  if (is.network(env$networks[[1]])) {
    env$directed <- network::is.directed(env$networks[[1]])  # directed?
    env$bipartite <- network::is.bipartite(env$networks[[1]])  # bipartite?
  } else {
    if (xergm.common::is.mat.directed(as.matrix(env$networks[[1]]))) {
      env$directed <- TRUE
    } else {
      env$directed <- FALSE
    }
    if (xergm.common::is.mat.onemode(as.matrix(env$networks[[1]]))) {
      env$bipartite <- FALSE
    } else {
      env$bipartite <- TRUE
    }
  }
  
  # adjust and disassemble formula
  env$form <- update.formula(formula, networks[[i]] ~ .)
  env$time.steps <- length(env$networks)
  tilde <- deparse(env$form[[1]])
  lhs <- deparse(env$form[[2]])
  rhs <- paste(deparse(env$form[[3]]), collapse = "")  # for long formulae
  rhs <- gsub("\\s+", " ", rhs)
  
  # parse rhs of formula and add indices to edgecov and dyadcov terms
  env$rhs.terms <- strsplit(rhs, "\\s*(\\+|\\*)\\s*")[[1]]
  rhs.indices <- gregexpr("\\+|\\*", rhs)[[1]]
  if (length(rhs.indices) == 1 && rhs.indices < 0) {
    rhs.operators <- character()
  } else {
    rhs.operators <- substring(rhs, rhs.indices, rhs.indices)
  }
  
  # preprocess dyadcov and edgecov terms, memory terms, and timecov terms
  covnames <- character()
  for (k in 1:length(env$rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", env$rhs.terms[k])) {  # edgecov or dyadcov
      
      # split up into components
      if (grepl(",\\s*?((attr)|\\\")", env$rhs.terms[k])) { # with attrib arg.
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)"
      } else { # without attribute argument
        s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)"
      }
      x1 <- sub(s, "\\1", env$rhs.terms[k], perl = TRUE)  # before the covariate
      x2 <- sub(s, "\\5", env$rhs.terms[k], perl = TRUE)  # name of the cov.
      if (grepl("\\[.*\\]", x2)) {
       stop(paste0("Covariate names are not allowed to have indices: ", x2, 
           ". Please prepare a list object before estimation."))
      }
      x3 <- sub(s, "\\6", env$rhs.terms[k], perl = TRUE)  # after the covariate
      x.current <- eval(parse(text = x2))
      type <- class(x.current)
      env$covnames <- c(env$covnames, x2)
      env[[x2]] <- x.current
      if (grepl("\\[i\\]+$", x2)) {
        stop(paste0("Error in the following model term: ", env$rhs.terms[k], 
            ". The index 'i' is used internally by btergm. Please use a ", 
            "different index, for example 'j'."))
      }
      
      # add brackets if necessary, convert to list, and reassemble rhs term
      if (grepl("[^\\]]\\]$", x2)) {
        # time-varying covariate with given indices (e.g., formula[1:5])
        env$rhs.terms[k] <- paste0(x1, x2, x3)
        if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix", 
            "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
          x.current <-list(x.current)
          env[[x2]] <- x.current
        }
        if (length(x.current) != env$time.steps) {
          stop(paste(x2, "has", length(x.current), "elements, but there are", 
              env$time.steps, "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
      } else if (type %in% c("matrix", "network", "dgCMatrix", "dgTMatrix", 
          "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
        # time-independent covariate
        if (!type %in% c("matrix", "network")) {
          x.current <- as.matrix(x.current)
        }
        env[[x2]] <- list()
        for (i in 1:env$time.steps) {
          env[[x2]][[i]] <- x.current
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
        env$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "list" || type == "network.list") {
        # time-varying covariate
        if (length(x.current) != env$time.steps) {
          stop(paste(x2, "has", length(get(x2)), "elements, but there are", 
              env$time.steps, "networks to be modeled."))
        }
        if (blockdiag == TRUE) {
          # do not add brackets
        } else {
          x2 <- paste0(x2, "[[i]]")
        }
        env$rhs.terms[k] <- paste0(x1, x2, x3)
      } else {  # something else --> try to convert to matrix list
        tryCatch(
          {
            env[[x2]] <- list(rep(as.matrix(x.current)), env$time.steps)
          }, 
          error = function(cond) {
            stop(paste0("Object '", x2, 
                "' could not be converted to a matrix."))
          }
        )
      }
    } else if (grepl("memory", env$rhs.terms[k])) {  # memory terms
      
      # extract type argument
      s <- "(?:memory\\((?:.*type\\s*=\\s*)?(?:\"|'))(\\w+)(?:(\"|').*\\))"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        type <- "stability"
      } else {
        type <- sub(s, "\\1", env$rhs.terms[k], perl = TRUE)
      }
      
      # extract lag argument
      s <- "(?:memory\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        lag <- 1
      } else {
        lag <- as.integer(sub(s, "\\1", env$rhs.terms[k], perl = TRUE))
      }
      if (lag > length(env$networks) - 1) {
        stop("The 'lag' argument in the 'memory' term is too large.")
      }
      
      # process dependent list of networks
      mem <- env$networks[-(length(env$networks):(length(env$networks) - lag + 
          1))]
      mem <- lapply(mem, as.matrix)
      memory <- list()
      for (i in 1:length(mem)) {
        if (type == "autoregression") {
          memory[[i]] <- mem[[i]]
        } else if (type == "stability") {
          mem[[i]][mem[[i]] == 0] <- -1
          memory[[i]] <- mem[[i]]
        } else if (type == "innovation") {
          memory[[i]] <- mem[[i]]
          memory[[i]][mem[[i]] == 0] <- 1
          memory[[i]][mem[[i]] == 1] <- 0
        } else if (type == "loss") {
          memory[[i]] <- mem[[i]]
          memory[[i]][mem[[i]] == 0] <- 0
          memory[[i]][mem[[i]] == 1] <- -1
        } else {
          stop("'type' argument in the 'memory' term not recognized.")
        }
      }
      rm(mem)
      
      # re-introduce as edgecov and name of model term including brackets
      env[["memory"]] <- memory
      if (blockdiag == TRUE) {
        env$rhs.terms[k] <- "edgecov(memory)"
      } else {
        env$rhs.terms[k] <- "edgecov(memory[[i]])"
      }
      env$covnames <- c(env$covnames, "memory")
    } else if (grepl("delrecip", env$rhs.terms[k])) {  # delayed reciprocity
      
      # extract mutuality argument
      s <- "(?:delrecip\\((?:.*mutuality\\s*=\\s*)?)((TRUE)|(FALSE)|T|F)(?:.*\\))"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        mutuality <- FALSE
      } else {
        mutuality <- as.logical(sub(s, "\\1", env$rhs.terms[k], perl = TRUE))
      }
      
      # extract lag argument
      s <- "(?:delrecip\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"  # get lag
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        lag <- 1
      } else {
        lag <- as.integer(sub(s, "\\1", env$rhs.terms[k], perl = TRUE))
      }
      if (lag > length(env$networks) - 1) {
        stop("The 'lag' argument in the 'delrecip' term is too large.")
      }
      
      # process dependent list of networks
      dlr <- env$networks[-(length(env$networks):(length(env$networks) - lag + 
          1))]
      dlr <- lapply(dlr, function(x) t(as.matrix(x)))
      delrecip <- list()
      for (i in 1:length(dlr)) {
        delrecip[[i]] <- dlr[[i]]
        if (mutuality == TRUE) {
          delrecip[[i]][dlr[[i]] == 0] <- -1
        }
      }
      rm(dlr)
      
      # re-introduce as edgecov and name of model term including brackets
      env[["delrecip"]] <- delrecip
      if (blockdiag == TRUE) {
        env$rhs.terms[k] <- "edgecov(delrecip)"
      } else {
        env$rhs.terms[k] <- "edgecov(delrecip[[i]])"
      }
      env$covnames <- c(env$covnames, "delrecip")
    } else if (grepl("timecov", env$rhs.terms[k])) {  # time covariate
      
      # extract x argument
      s <- "(?:timecov\\((?:.*x\\s*=\\s*)?)(\\w+)(?:.*\\))"
      if (sub(s, "\\1", env$rhs.terms[k], perl = TRUE) %in% c("minimum", 
          "maximum", "transform", "min", "max", "trans")) {
        s <- "(?:timecov\\(?:.*x\\s*=\\s*)(\\w+)(?:.*\\))"
      }
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        x <- NULL
        suffix <- ""
        label <- "timecov"
      } else {
        x <- sub(s, "\\1", env$rhs.terms[k], perl = TRUE)
        label <- paste0("timecov.", x)
      }
      
      # extract minimum argument
      s <- "(?:timecov\\(.*minimum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        minimum <- 1
      } else {
        minimum <- as.integer(sub(s, "\\1", env$rhs.terms[k], perl = TRUE))
      }
      
      # extract maximum argument
      s <- "(?:timecov\\(.*maximum\\s*=\\s*)(\\d+)(?:.*\\))"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        maximum <- env$time.steps
      } else {
        maximum <- as.integer(sub(s, "\\1", env$rhs.terms[k], perl = TRUE))
      }
      
      # extract transform argument
      s <- "(?:timecov\\(.*transform\\s*=\\s*)(.+?)(?:(?:,|\\)$)]*.*)"
      if (grepl(s, env$rhs.terms[k]) == FALSE) {
        transform <- function(t) 1 + (0 * t) + (0 * t^2)
      } else {
        transform <- eval(parse(text = sub(s, "\\1", env$rhs.terms[k], 
            perl = TRUE)))
      }
      
      # process dependent list of networks
      if (is.null(x)) {
        covariate <- env[["networks"]]
        onlytime <- TRUE
      } else {
        onlytime <- FALSE
        covariate <- get(x)
      }
      tc <- timecov(covariate = covariate, minimum = minimum, 
        maximum = maximum, transform = transform, onlytime = onlytime)
      
      # re-introduce as edgecov and name of model term including brackets
      env[[label]] <- tc
      if (blockdiag == TRUE) {
        env$rhs.terms[k] <- paste0("edgecov(", label, ")")
      } else {
        env$rhs.terms[k] <- paste0("edgecov(", label, "[[i]])")
      }
      env$covnames <- c(env$covnames, label)
    }
  }
  env$covnames <- c("networks", env$covnames)
  
  # fix different lengths of DV/covariate lists due to temporal dependencies
  lengths <- sapply(env$covnames, function(cn) length(env[[cn]]))
  mn <- max(lengths)
  if (length(table(lengths)) > 1) {
    mn <- min(lengths)
    env$time.steps <- mn
    for (i in 1:length(env$covnames)) {
      cn <- env$covnames[[i]]
      l <- env[[cn]]
      difference <- length(l) - mn
      if (difference > 0) {
        env[[cn]] <- l[(difference + 1):length(l)]
      }
    }
  }
  t.end <- max(lengths)
  t.start <- t.end - mn + 1
  
  # determine and report initial dimensions of networks and covariates
  if (verbose == TRUE) {
    if (length(env$covnames) > 1) {
      dimensions <- lapply(lapply(env$covnames, function(x) env[[x]]), 
          function(y) sapply(y, function(z) dim(as.matrix(z))))
      rownames(dimensions[[1]]) <- paste(env$lhs.original, c("(row)", "(col)"))
      for (i in 2:length(dimensions)) {
        rownames(dimensions[[i]]) <- c(paste(env$covnames[i], "(row)"), 
            paste(env$covnames[i], "(col)"))
      }
      dimensions <- do.call(rbind, dimensions)
      colnames(dimensions) <- paste0("t=", t.start:t.end) #1:length(env$networks))
      message("\nInitial dimensions of the network and covariates:")
      print(dimensions)
    } else {
      message("\nNo covariates provided.")
    }
  }
  
  # determine whether covariate dimensions need to be automatically adjusted
  env$auto.adjust <- FALSE
  if (length(env$covnames) > 1) {
    # check number of rows and columns
    nr <- lapply(lapply(env$covnames, function(x) env[[x]]), 
        function(y) sapply(y, function(z) nrow(as.matrix(z))))
    nr <- do.call(rbind, nr)
    nc <- lapply(lapply(env$covnames, function(x) env[[x]]), 
        function(y) sapply(y, function(z) ncol(as.matrix(z))))
    nc <- do.call(rbind, nc)
    for (i in 1:ncol(nr)) {
      if (length(unique(nr[, i])) > 1) {
        env$auto.adjust <- TRUE
      }
    }
    for (i in 1:ncol(nc)) {
      if (length(unique(nc[, i])) > 1) {
        env$auto.adjust <- TRUE
      }
    }
    if (verbose == TRUE && env$auto.adjust == TRUE) {
       message(paste("\nDimensions differ across networks within time steps."))
    }
    # check if labels are present
    if (env$auto.adjust == TRUE) {
      for (i in 1:length(env$covnames)) {
        for (t in 1:env$time.steps) {
          if (is.null(rownames(as.matrix(env[[env$covnames[i]]][[t]]))) || 
              is.null(colnames(as.matrix(env[[env$covnames[i]]][[t]])))) {
            stop(paste0("The dimensions of the covariates differ, but ", 
                "covariate '", env$covnames[i], 
                " does not have node labels at t = ", t, 
                ". Automatic adjustment of dimensions is therefore not ", 
                "possible."))
          }
        }
      }
    }
    # check if there are different labels despite identical dimensions
    if (env$auto.adjust == FALSE) {
      for (t in 1:env$time.steps) {
        rlabels.i <- list()
        clabels.i <- list()
        for (i in 1:length(env$covnames)) {
          rlabels.i[[i]] <- rownames(as.matrix(env[[env$covnames[i]]][[t]]))
          clabels.i[[i]] <- colnames(as.matrix(env[[env$covnames[i]]][[t]]))
        }
        rlabels.i <- do.call(rbind, rlabels.i)
        clabels.i <- do.call(rbind, clabels.i)
        flag <- FALSE
        if (!is.null(rlabels.i)) {
          for (j in 1:ncol(rlabels.i)) {
            if (length(unique(rlabels.i[, j])) > 1) {
              env$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
        if (!is.null(clabels.i)) {
          for (j in 1:ncol(clabels.i)) {
            if (length(unique(clabels.i[, j])) > 1) {
              env$auto.adjust <- TRUE
              flag <- TRUE
              break
            }
          }
        }
      }
      if (verbose == TRUE && flag == TRUE) {
        message(paste("\nSame dimensions but different labels across", 
            "networks within time steps."))
      }
    }
  }
  if (verbose == TRUE && env$auto.adjust == TRUE) {
    message("Trying to auto-adjust the dimensions of the networks. ", 
        "If this fails, provide conformable matrices or network objects.")
  } else if (verbose == TRUE) {
    message("\nAll networks are conformable.")
  }
  
  # do mutual adjustment of networks and covariates at each time step
  structzero.df <- data.frame(label = character(), time = integer(), 
      object = character(), where = character())
  if (length(env$covnames) > 0 && env$auto.adjust == TRUE) {
    for (i in 1:env$time.steps) {
      for (j in 1:length(env$covnames)) {
        for (k in 1:length(env$covnames)) {
          if (j != k) {
            nw.j <- env[[env$covnames[j]]][[i]]
            rn.j <- rownames(as.matrix(nw.j))
            cn.j <- colnames(as.matrix(nw.j))
            nr.j <- nrow(as.matrix(nw.j))
            nc.j <- ncol(as.matrix(nw.j))
            nw.k <- env[[env$covnames[k]]][[i]]
            rn.k <- rownames(as.matrix(nw.k))
            cn.k <- colnames(as.matrix(nw.k))
            nr.k <- nrow(as.matrix(nw.k))
            nc.k <- ncol(as.matrix(nw.k))
            if (is.null(rn.j) || is.null(cn.j)) {
              stop(paste0("Missing row or column labels in object '", 
                  env$covnames[j], "'. Provide row and column ", 
                  "labels for all networks and covariates."))
            } else if (is.null(rn.k) || is.null(cn.k)) {
              stop(paste0("Missing row or column labels in object '", 
                  env$covnames[k], "'. Provide row and column ", 
                  "labels for all networks and covariates."))
            } else {
              if (is.null(rn.j) && !is.null(rn.k) && nr.j == nr.k) {
                if (class(nw.j) == "network") {
                  network::set.vertex.attribute(nw.j, "vertex.names", rn.k)
                } else {
                  rownames(nw.j) <- rn.k
                }
              } else if (is.null(rn.k) && !is.null(rn.j) && nr.j == nr.k) {
                if (class(nw.k) == "network") {
                  network::set.vertex.attribute(nw.k, "vertex.names", rn.j)
                } else {
                  rownames(nw.k) <- rn.j
                }
              } else if ((is.null(rn.k) || is.null(rn.j)) && nr.j != nr.k) {
                stop(paste0("Object '", env$covnames[j], 
                    "' is incompatible with object '", env$covnames[k], 
                    "' at t = ", i, "."))
              }
              # adjust j to k
              nw.j.labels <- adjust(nw.j, nw.k, remove = FALSE, 
                  value = 1, returnlabels = TRUE)
              nw.j <- adjust(nw.j, nw.k, remove = FALSE, value = 1)
              env[[env$covnames[j]]][[i]] <- nw.j
              ro <- nw.j.labels$added.row
              co <- nw.j.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)), 
                    object = rep(env$covnames[j], length(ro)), 
                    where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)), 
                    object = rep(env$covnames[j], length(co)), 
                    where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
              # adjust k back to j
              nw.k.labels <- adjust(nw.k, nw.j, remove = FALSE, 
                  value = 1, returnlabels = TRUE)
              nw.k <- adjust(nw.k, nw.j, remove = FALSE, value = 1)
              env[[env$covnames[k]]][[i]] <- nw.k
              ro <- nw.k.labels$added.row
              co <- nw.k.labels$added.col
              if (length(ro) > 0) {
                ro <- data.frame(label = ro, time = rep(i, length(ro)), 
                    object = rep(env$covnames[j], length(ro)), 
                    where = rep("row", length(ro)))
                structzero.df <- rbind(structzero.df, ro)
              }
              if (length(co) > 0) {
                co <- data.frame(label = co, time = rep(i, length(co)), 
                    object = rep(env$covnames[j], length(co)), 
                    where = rep("col", length(co)))
                structzero.df <- rbind(structzero.df, co)
              }
            }
          }
        }
      }
    }
  }
  
  # check whether all dimensions are cross-sectionally conformable now
  nr.net <- sapply(env$networks, function(x) nrow(as.matrix(x)))
  for (i in 1:length(env$covnames)) {
    nr <- sapply(env[[env$covnames[i]]], function(x) {
      nrow(as.matrix(x))
    })
    for (j in 1:env$time.steps) {
      if (nr[j] != nr.net[j]) {
        stop(paste0("Covariate object '", env$covnames[i], 
            "' does not have the same number of rows as the dependent ", 
            "network at time step ", j, "."))
      }
    }
  }
  nc.net <- sapply(env$networks, function(x) ncol(as.matrix(x)))
  for (i in 1:length(env$covnames)) {
    nc <- sapply(env[[env$covnames[i]]], function(x) {
      ncol(as.matrix(x))
    })
    for (j in 1:env$time.steps) {
      if (nc[j] != nc.net[j]) {
        stop(paste0("Covariate object '", env$covnames[i], 
            "' does not have the same number of columns as the dependent ", 
            "network at time step ", j, "."))
      }
    }
  }
  
  # reporting
  if (verbose == TRUE) {
    if (env$auto.adjust == TRUE) {
      sz.row <- unique(structzero.df[structzero.df$where == "row", -3])
      szrownum <- numeric(length(env$networks))
      for (i in 1:length(env$networks)) {
        szrownum[i] <- nrow(sz.row[sz.row$time == i, ])
      }
      sz.col <- unique(structzero.df[structzero.df$where == "col", -3])
      szcolnum <- numeric(length(env$networks))
      for (i in 1:length(env$networks)) {
        szcolnum[i] <- nrow(sz.col[sz.col$time == i, ])
      }
      totrow <- sapply(env$networks, function(x) nrow(as.matrix(x)))
      totcol <- sapply(env$networks, function(x) ncol(as.matrix(x)))
      if (offset == TRUE) {
        dimensions <- rbind(totrow, totcol, szrownum, szcolnum, 
            totrow - szrownum, totcol - szcolnum)
        rownames(dimensions) <- c("total number of rows", 
            "total number of columns", "row-wise structural zeros", 
            "column-wise structural zeros", "remaining rows", 
            "remaining columns")
      } else {
        dimensions <- rbind(szrownum, szcolnum, totrow - szrownum, 
            totcol - szcolnum)
        rownames(dimensions) <- c("maximum deleted nodes (row)", 
            "maximum deleted nodes (col)", "remaining rows", 
            "remaining columns")
      }
      colnames(dimensions) <- paste0("t=", t.start:t.end)
      if (nrow(structzero.df) > 0) {
        if (offset == TRUE) {
          message("\nNodes affected completely by structural zeros:")
        } else {
          message("\nAbsent nodes:")
        }
        szcopy <- structzero.df
        szcopy$time <- szcopy$time - 1 + t.start  # correct lagged starting time
        print(unique(szcopy))
      } else {
        message("\nAll nodes are retained.")
      }
      
      message("\nNumber of nodes per time step after adjustment:")
      print(dimensions)
      
    }
  }
  
  # create list of offset matrices (required both for offset and node removal)
  env$offsmat <- list()
  for (i in 1:env$time.steps) {
    mat <- matrix(0, nrow = nrow(as.matrix(env$networks[[i]])), 
        ncol = ncol(as.matrix(env$networks[[i]])))
    rownames(mat) <- rownames(as.matrix(env$networks[[i]]))
    colnames(mat) <- colnames(as.matrix(env$networks[[i]]))
    env$offsmat[[i]] <- mat
  }
  if (nrow(structzero.df) > 0) {
    for (i in 1:nrow(structzero.df)) {
      if (structzero.df$where[i] == "row") {
        index <- which(rownames(env$offsmat[[structzero.df$time[i]]]) == 
            structzero.df$label[i])
        env$offsmat[[structzero.df$time[i]]][index, ] <- 1
      } else {
        index <- which(colnames(env$offsmat[[structzero.df$time[i]]]) == 
            structzero.df$label[i])
        env$offsmat[[structzero.df$time[i]]][, index] <- 1
      }
    }
  }
  
  # offset preparation or node removal for MPLE
  if (offset == TRUE) {
    # add offset to formula and reassemble formula
    env$rhs.terms[length(env$rhs.terms) + 1] <- "offset(edgecov(offsmat[[i]]))"
    rhs.operators[length(rhs.operators) + 1] <- "+"
  } else {
    # delete nodes with structural zeros
    if (env$auto.adjust == TRUE) {
      env$offsmat <- suppressMessages(handleMissings(env$offsmat, na = 1, 
          method = "remove"))
      for (j in 1:length(env$covnames)) {
        env[[env$covnames[j]]] <- adjust(env[[env$covnames[j]]], env$offsmat)
      }
    }
  }
  
  # determine and report initial dimensions of networks and covariates
  if (verbose == TRUE && length(env$covnames) > 1) {
    dimensions <- lapply(lapply(env$covnames, function(x) env[[x]]), 
        function(y) sapply(y, function(z) dim(as.matrix(z))))
    rownames(dimensions[[1]]) <- paste(env$lhs.original, c("(row)", "(col)"))
    for (i in 2:length(dimensions)) {
      rownames(dimensions[[i]]) <- c(paste(env$covnames[i], "(row)"), 
          paste(env$covnames[i], "(col)"))
    }
    dimensions <- do.call(rbind, dimensions)
    colnames(dimensions) <- paste0("t=", t.start:t.end) #1:length(env$networks))
    message("\nDimensions of the network and covariates after adjustment:")
    print(dimensions)
  }
  
  # assemble formula
  rhs <- env$rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], env$rhs.terms[i + 1])
    }
  }
  f <- paste(lhs, tilde, rhs)
  env$form <- as.formula(f, env = env)
  
  # for mtergm estimation using MCMC: create block-diagonal matrices
  if (blockdiag == TRUE) {
    if (env$bipartite == TRUE) {
      stop(paste("MCMC estimation is currently only supported for one-mode", 
          "networks. Use the btergm function instead."))
    }
    # also save formula without time indices for ergm estimation
    env$form <- update.formula(env$form, networks ~ .)
    env$form <- paste(deparse(env$form), collapse = "")
    env$form <- paste(env$form, "+ offset(edgecov(offsmat))")
    env$form <- as.formula(env$form, env = env)
    # make covariates block-diagonal
    if (length(env$covnames) > 1) {
      for (j in 2:length(env$covnames)) {
        env[[env$covnames[j]]] <- as.matrix(Matrix::bdiag(lapply(
            env[[env$covnames[j]]], as.matrix)))
      }
    }
    # create block-diagonal offset matrix and merge with existing offsmat term
    env$offsmat <- as.matrix(Matrix::bdiag(env$offsmat))  # make block-diagonal
    bdoffset <- lapply(env$networks, as.matrix)
    for (i in 1:length(bdoffset)) {
      bdoffset[[i]][, ] <- 1
    }
    bdoffset <- as.matrix((Matrix::bdiag(bdoffset) - 1) * -1)  # off-diagonal
    env$offsmat <- env$offsmat + bdoffset
    rm(bdoffset)
    env$offsmat[env$offsmat > 0] <- 1
    # make dependent network block-diagonal
    if (class(env$networks[[1]]) == "network") {  # network
      attrnames <- network::list.vertex.attributes(env$networks[[1]])
      attributes <- list()
      for (i in 1:length(env$networks)) {
        attrib <- list()
        for (j in 1:length(attrnames)) {
          attrib[[j]] <- network::get.vertex.attribute(env$networks[[i]], 
              attrnames[j])
        }
        attributes[[i]] <- attrib
        env$networks[[i]] <- as.matrix(env$networks[[i]])
      }
      env$networks <- network::network(as.matrix(Matrix::bdiag(env$networks)), 
          directed = env$directed, bipartite = env$bipartite)
      for (i in 1:length(attrnames)) {  # collate attributes and merge back in
        attrib <- unlist(lapply(attributes, function(x) x[[i]]))
        network::set.vertex.attribute(env$networks, attrnames[i], attrib)
      }
    } else {  # matrix
      env$networks <- network::network(as.matrix(Matrix::bdiag(env$networks)), 
          directed = env$directed, bipartite = env$bipartite)
    }
    if (verbose == TRUE) {
      cat("\n")  # to get a blank line before the MCMC MLE output starts
    }
  }
  return(env)  # return the environment with all the data
}

