# This file contains preprocessing functions for btergm and tnam.

# check if a matrix is a one-mode matrix
is.mat.onemode <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat)) 
      && any(rownames(mat) != colnames(mat))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# check if a matrix represents a directed network
is.mat.directed <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    return(FALSE)
  } else if (!is.null(rownames(mat)) && !is.null(colnames(mat)) 
      && any(rownames(mat) != colnames(mat), na.rm = TRUE)) {
    return(FALSE)
  } else {
    if (any(as.matrix(mat) != t(as.matrix(mat)), na.rm = TRUE)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


# how many NAs are there per row or column?
numMissing <- function(mat, type = "both", na = NA) {
  numrow <- apply(as.matrix(mat), 1, function(x) sum(x %in% na))
  numcol <- apply(as.matrix(mat), 2, function(x) sum(x %in% na))
  if (type == "both") {
    return(numrow + numcol)
  } else if (type == "row") {
    return(numrow)
  } else if (type == "col") {
    return(numcol)
  } else {
    stop("Unknown 'type' argument in the 'numMissing' function.")
  }
}


# process NA values (= remove nodes with NAs iteratively)
handleMissings <- function(mat, na = NA, method = "remove", logical = FALSE) {
  
  # check and convert arguments
  if (is.null(mat)) {
    stop("The 'mat' argument is not valid.")
  } else if (class(mat) == "list") {
    # OK; do nothing; check later in next step
    initialtype <- "list"
  } else if (class(mat) %in% c("matrix", "network", "data.frame")) {
    # wrap in list
    initialtype <- class(mat)
    mat <- list(mat)
  } else if (length(mat) > 1) {
    # vector --> wrap in list
    initialtype <- class(mat)
    mat <- list(mat)
  } else if (is.function(mat)) {
    stop(paste("The input object is a function. Did you choose a name", 
        "for the input object which already exists as a function in the", 
        "workspace?"))
  } else {
    stop("The 'mat' argument is not valid.")
  }
  
  onemode <- list()  # will indicate whether it is a one- or two-mode network
  directed <- list()  # will indicate whether the network is directed
  attribnames <- list()  # will contain the names of nodal attributes
  attributes <- list()  # will contain the nodal attributes at each time step
  type <- list()  # will indicate the type of data structure at time step i
  for (i in 1:length(mat)) {
    if (class(mat[[i]]) == "matrix") {
      # check manually if onemode and directed
      onemode[[i]] <- is.mat.onemode(mat[[i]])  # helper function
      directed[[i]] <- is.mat.directed(mat[[i]])  # helper function
      type[[i]] <- "matrix"
    } else if (class(mat[[i]]) == "network") {
      # save onemode and directed information; save attributes for later use
      if (is.bipartite(mat[[i]])) {
        onemode[[i]] <- FALSE
      } else {
        onemode[[i]] <- TRUE
      }
      if (is.directed(mat[[i]])) {
        directed[[i]] <- TRUE
      } else {
        directed[[i]] <- FALSE
      }
      attribnames[[i]] <- list.vertex.attributes(mat[[i]])
      attrib <- list()  # list of attributes at time i
      if (network.size(mat[[i]]) == 0) {
        attribnames[[i]] <- character()
        attributes[[i]] <- list(character())
      } else {
        for (j in 1:length(attribnames[[i]])) {
          attrib[[j]] <- get.vertex.attribute(mat[[i]], attribnames[[i]][j])
        }
      }
      attributes[[i]] <- attrib
      mat[[i]] <- as.matrix(mat[[i]])
      type[[i]] <- "network"
    } else if (class(mat[[i]]) == "data.frame") {
      type[[i]] <- "data.frame"
    } else {
      type[[i]] <- class(mat[[i]])
    }
  }
  
  if (is.null(logical) || !is.logical(logical) || length(logical) > 1) {
    stop("The 'logical' argument should be either TRUE or FALSE.")
  }
  if (is.null(method) || class(method) != "character") {
    stop("The 'method' argument should be a character object.")
  }
  if (length(method) > 1) {
    method <- method[1]
  }
  
  na.mat <- list()  # will contain matrices indicating which values are NA
  for (i in 1:length(mat)) {
    na.mat[[i]] <- apply(mat[[i]], 1:2, function(x) x %in% NA)
    if (length(mat) == 1) {  # used for reporting later
      time <- ""
    } else {
      time <- paste0("t = ", i, ": ")
    }
    if (class(mat[[i]]) == "matrix") {
      # matrix objects
      # replace by real NAs, then count NAs
      obs <- length(mat[[i]])
      missing.abs <- length(which(mat[[i]] %in% na))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)
      
      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        nwunique <- unique(as.numeric(mat[[i]][!mat[[i]] %in% na]))
        nwmode <- nwunique[which.max(tabulate(match(mat[[i]][!mat[[i]] %in% 
            na], nwunique)))]
        mat[[i]][mat[[i]] %in% na] <- nwmode
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " ties) were replaced by the mode (", nwmode, 
            ") because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        mat[[i]][mat[[i]] %in% na] <- 0
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " ties) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove rows and columns with NA values iteratively
        rowLabels <- rownames(mat[[i]])
        colLabels <- colnames(mat[[i]])
        if (onemode[[i]] == TRUE) {
          while(sum(numMissing(mat[[i]], na = na)) > 0) {
            indices <- which(numMissing(mat[[i]], na = na) == 
                max(numMissing(mat[[i]], na = na)))
            mat[[i]] <- mat[[i]][-indices, -indices]
            rowLabels <- rowLabels[-indices]
            colLabels <- colLabels[-indices]
            na.mat[[i]][indices, ] <- TRUE
            na.mat[[i]][, indices] <- TRUE
            if (type[[i]] == "network") {
              if (length(attribnames[[i]]) > 0) {
                for (j in 1:length(attribnames[[i]])) {
                  attributes[[i]][[j]] <- attributes[[i]][[j]][-indices]
                }
              }
            }
          }
        } else {
          while(sum(numMissing(mat[[i]], type = "row", na = na)) + 
              sum(numMissing(mat[[i]], type = "col", na = na)) > 0) {
            rowNAs <- numMissing(mat[[i]], type = "row", na = na)
            colNAs <- numMissing(mat[[i]], type = "col", na = na)
            maxNA <- max(c(rowNAs, colNAs))
            if (length(which(rowNAs == maxNA)) > 0) {
              indices <- which(rowNAs == maxNA)
              mat[[i]] <- mat[[i]][-indices, ]
              rowLabels <- rowLabels[-indices]
              na.mat[[i]][indices, ] <- TRUE
              if (type[[i]] == "network") {
                if (length(attribnames[[i]]) > 0) {
                  for (j in 1:length(attribnames[[i]])) {
                    attributes[[i]][[j]] <- attributes[[i]][[j]][-indices]
                  }
                }
              }
            } else if (length(which(colNAs == maxNA)) > 0) {
              indices <- which(colNAs == maxNA)
              mat[[i]] <- mat[[i]][, -indices]
              colLabels <- colLabels[-indices]
              na.mat[[i]][, indices] <- TRUE
              # in bipartite networks, attributes for rows and columns are 
              # saved in a single vector consecutively
              indices.bip <- nrow(mat[[i]]) + indices
              if (type[[i]] == "network") {
                if (length(attribnames[[i]]) > 0) {
                  for (j in 1:length(attribnames[[i]])) {
                    attributes[[i]][[j]] <- attributes[[i]][[j]][-indices.bip]
                  }
                }
              }
            }
          }
        }
        rownames(mat[[i]]) <- rowLabels
        colnames(mat[[i]]) <- colLabels
        removed.abs <- obs - length(mat[[i]])
        removed.perc <- round(100 * removed.abs / obs, digits = 2)
        if (is.nan(removed.perc)) {
          removed.perc <- 0
        }
        if (is.nan(missing.perc)) {
          missing.perc <- 0
        }
        message(paste0("t = ", i, ": ", removed.perc, "% of the data (= ", 
            removed.abs, " ties) were dropped due to ", missing.perc, "% (= ", 
            missing.abs, ") missing ties."))
      } else {
        stop("Method not supported.")
      }
      
      # convert back into network if initial item was a network
      if (type[[i]] == "network") {
        bip <- (onemode[[i]] == FALSE)
        mat[[i]] <- network(mat[[i]], directed = directed[[i]], 
            bipartite = bip)
        if (length(attribnames[[i]]) > 0) {
          for (j in 1:length(attribnames[[i]])) {
            mat[[i]] <- set.vertex.attribute(mat[[i]], attribnames[[i]][j], 
                attributes[[i]][[j]])
          }
        }
      }
    } else if (class(mat[[i]]) == "data.frame") {
      # data.frame objects
      # replace by real NAs, then count NAs
      for (j in 1:nrow(mat[[i]])) {
        for (k in 1:ncol(mat[[i]])) {
          if (mat[[i]][j, k] %in% na) {
            mat[[i]][j, k] <- NA
          }
        }
      }
      obs <- nrow(mat[[i]]) * ncol(mat[[i]])
      missing.abs <- length(which(is.na(mat[[i]])))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)
      
      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        for (j in 1:ncol(mat[[i]])) {
          if (is.numeric(mat[[i]][, j])) {
            nwunique <- unique(as.numeric(mat[[i]][, j]))
            nwmode <- nwunique[which.max(tabulate(match(mat[[i]][, j], 
                nwunique)))]
            mat[[i]][, j][is.na(mat[[i]][, j])] <- nwmode
          }
        }
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " ties) were replaced by the mode in the respective ", 
            "column because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        for (j in 1:ncol(mat[[i]])) {
          mat[[i]][, j][is.na(mat[[i]][, j])] <- 0
        }
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " elements) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove rows with NA values
        before <- nrow(mat[[i]])
        mat[[i]] <- mat[[i]][complete.cases(mat[[i]]), ]
        after <- nrow(mat[[i]])
        removed <- before - after
        rem.perc <- 100 * (1 - after / before)
        message(paste0("t = ", i, ": ", removed, " rows (", rem.perc, 
            "% of all rows) were removed due to missing elements."))
      } else {
        stop("Method not supported.")
      }
    } else if (length(mat[[i]]) > 1) {
      # vectors of arbitrary content
      mat[[i]][mat[[i]] %in% na] <- NA
      obs <- length(mat[[i]])
      missing.abs <- length(which(is.na(mat[[i]])))
      missing.perc <- round(100 * missing.abs / obs, digits = 2)
      
      # do the actual work
      if (method == "fillmode") {
        # fill with modal value (often 0 but not always)
        if (!is.numeric(mat[[i]])) {
          stop("'fillmode' is only compatible with numeric objects.")
        }
        nwunique <- unique(as.numeric(mat[[i]]))
        nwmode <- nwunique[which.max(tabulate(match(mat[[i]], nwunique)))]
        mat[[i]][is.na(mat[[i]])] <- nwmode
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " ties) were replaced by the mode in the respective ", 
            "column because they were NA."))
      } else if (method == "zero") {
        # impute 0 when NA
        mat[[i]][is.na(mat[[i]])] <- 0
        message(paste0("t = ", i, ": ", missing.perc, "% of the data (= ", 
            missing.abs, " ties) were replaced by 0 because they were NA."))
      } else if (method == "remove") {
        # remove NA values
        mat[[i]] <- mat[[i]][!is.na(mat[[i]])]
        message(paste0(time, missing.perc, "% of the data (= ", 
            missing.abs, " elements) were removed because they were NA."))
      } else {
        stop("Method not supported.")
      }
    }
  }
  
  if (logical == TRUE) {
    if (length(na.mat) == 1 && initialtype != "list") {
      return(na.mat[[1]])
    } else {
      return(na.mat)
    }
  } else {
    if (length(mat) == 1 && initialtype != "list") {
      return(mat[[1]])
    } else {
      return(mat)
    }
  }
}


# adjust the dimensions of a source object to the dimensions of a target object
adjust <- function(source, target, remove = TRUE, add = TRUE, value = NA, 
    returnlabels = FALSE) {
  
  # make sure the source is a list
  if (is.null(source)) {
    stop("The 'source' argument was not recognized.")
  } else if (class(source) == "matrix") {
    # wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "matrix"
  } else if (class(source) == "network") {
    # wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "network"
  } else if (class(source) == "list") {
    # rename
    sources <- source
    sources.initialtype <- "list"
  } else if (is.vector(source)) {
    # vector of some type; wrap in list
    sources <- list()
    sources[[1]] <- source
    sources.initialtype <- "vector"
  } else {
    stop(paste("Source data type not supported. Supported types are 'matrix',", 
        "'network', and 'list' objects and vectors."))
  }
  
  # make sure the target is a list
  if (is.null(target)) {
    stop("The 'target' argument was not recognized.")
  } else if (class(target) == "matrix") {
    # wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "matrix"
  } else if (class(target) == "network") {
    # wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "network"
  } else if (class(target) == "list") {
    # rename
    targets <- target
    targets.initialtype <- "list"
  } else if (is.vector(target)) {
    # vector of some type; wrap in list
    targets <- list()
    targets[[1]] <- target
    targets.initialtype <- "vector"
  } else {
    stop(paste("Target data type not supported. Supported types are 'matrix',", 
        "'network', and 'list' objects and vectors."))
  }
  
  # make sure that both lists (sources and targets) have the same length
  if (length(sources) == length(targets)) {
    # OK; do nothing
  } else if (length(sources) == 1) {
    for (i in 2:length(targets)) {
      sources[[i]] <- sources[[1]]
    }
  } else if (length(targets) == 1) {
    for (i in 2:length(sources)) {
      targets[[i]] <- targets[[1]]
    }
  } else {
    stop("Different numbers of sources and targets were provided.")
  }
  
  # convert each item if necessary and save nodal attributes
  sources.attribnames <- list()
  sources.attributes <- list()
  sources.types <- list()
  sources.onemode <- list()
  sources.directed <- list()
  targets.attribnames <- list()
  targets.attributes <- list()
  targets.types <- list()
  targets.onemode <- list()
  targets.directed <- list()
  for (i in 1:length(sources)) {
    sources.types[[i]] <- class(sources[[i]])
    if (class(sources[[i]]) == "network") {
      # save source attributes and other meta information in list
      sources.attribnames[[i]] <- list.vertex.attributes(sources[[i]])
      attributes <- list()
      if (length(sources.attribnames) > 0) {
        for (j in 1:length(sources.attribnames[[i]])) {
          attributes[[j]] <- get.vertex.attribute(sources[[i]], 
              sources.attribnames[[i]][j])
        }
      }
      sources.attributes[[i]] <- attributes
      sources.onemode[[i]] <- !is.bipartite(sources[[i]])
      sources.directed[[i]] <- is.directed(sources[[i]])
      sources[[i]] <- as.matrix(sources[[i]])  # convert to matrix
    } else if (class(sources[[i]]) == "matrix") {
      sources.onemode[[i]] <- is.mat.onemode(sources[[i]])
      sources.directed[[i]] <- is.mat.directed(sources[[i]])
    } else {
      sources[[i]] <- as.matrix(sources[[i]], ncol = 1)
    }
    
    targets.types[[i]] <- class(targets[[i]])
    if (class(targets[[i]]) == "network") {
      # save target attributes and other meta information in list
      targets.attribnames[[i]] <- list.vertex.attributes(targets[[i]])
      attributes <- list()
      if (length(targets.attribnames) > 0) {
        for (j in 1:length(targets.attribnames[[i]])) {
          attributes[[j]] <- get.vertex.attribute(targets[[i]], 
              targets.attribnames[[i]][j])
        }
      }
      targets.attributes[[i]] <- attributes
      targets.onemode[[i]] <- !is.bipartite(targets[[i]])
      targets.directed[[i]] <- is.directed(targets[[i]])
      targets[[i]] <- as.matrix(targets[[i]])  # convert to matrix
    } else if (class(targets[[i]]) == "matrix") {
      targets.onemode[[i]] <- is.mat.onemode(targets[[i]])
      targets.directed[[i]] <- is.mat.directed(targets[[i]])
    } else {
      targets[[i]] <- as.matrix(targets[[i]], ncol = 1)
    }
  }
  
  # impute row or column labels if only one of them is present
  for (i in 1:length(sources)) {
    if (is.null(rownames(sources[[i]])) && !is.null(colnames(sources[[i]])) && 
        nrow(sources[[i]]) == ncol(sources[[i]])) {
      rownames(sources[[i]]) <- colnames(sources[[i]])
    }
    if (is.null(colnames(sources[[i]])) && !is.null(rownames(sources[[i]])) && 
        nrow(sources[[i]]) == ncol(sources[[i]])) {
      colnames(sources[[i]]) <- rownames(sources[[i]])
    }
    if (is.null(rownames(targets[[i]])) && !is.null(colnames(targets[[i]])) && 
        nrow(targets[[i]]) == ncol(targets[[i]])) {
      rownames(targets[[i]]) <- colnames(targets[[i]])
    }
    if (is.null(colnames(targets[[i]])) && !is.null(rownames(targets[[i]])) && 
        nrow(targets[[i]]) == ncol(targets[[i]])) {
      colnames(targets[[i]]) <- rownames(targets[[i]])
    }
  }
  
  # throw error if there are duplicate names (first sources, then targets)
  for (i in 1:length(sources)) {
    if (class(sources[[i]]) %in% c("matrix", "data.frame")) {
      # row names
      if (!is.null(rownames(sources[[i]]))) {
        test.actual <- nrow(sources[[i]])
        test.unique <- length(unique(rownames(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate source row names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate source row name."))
        }
      }
      # column names
      if (!is.null(colnames(sources[[i]]))) {
        test.actual <- ncol(sources[[i]])
        test.unique <- length(unique(colnames(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate source column names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate source column name."))
        }
      }
    } else {
      # vector names
      if (!is.null(names(sources[[i]]))) {
        test.actual <- length(sources[[i]])
        test.unique <- length(unique(names(sources[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate source names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate source name."))
        }
      }
    }
  }
  for (i in 1:length(targets)) {
    if (class(targets[[i]]) %in% c("matrix", "data.frame")) {
      # row names
      if (!is.null(rownames(targets[[i]]))) {
        test.actual <- nrow(targets[[i]])
        test.unique <- length(unique(rownames(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate target row names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate target row name."))
        }
      }
      # column names
      if (!is.null(colnames(targets[[i]]))) {
        test.actual <- ncol(targets[[i]])
        test.unique <- length(unique(colnames(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate target column names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate target column name."))
        }
      }
    } else {
      # vector names
      if (!is.null(names(targets[[i]]))) {
        test.actual <- length(targets[[i]])
        test.unique <- length(unique(names(targets[[i]])))
        dif <- test.actual - test.unique
        if (dif > 1) {
          stop(paste0("At t = ", i, ", there are ", dif, 
              " duplicate target names."))
        } else if (dif == 1) {
          stop(paste0("At t = ", i, ", there is ", dif, 
              " duplicate target name."))
        }
      }
    }
  }
  
  # go through sources and targets and do the actual adjustment
  for (i in 1:length(sources)) {
    if (!is.vector(sources[[i]]) && !class(sources[[i]]) %in% c("matrix", 
        "network")) {
      stop(paste("Source item", i, "is not a matrix, network, or vector."))
    }
    if (!is.vector(targets[[i]]) && !class(targets[[i]]) %in% c("matrix", 
        "network")) {
      stop(paste("Target item", i, "is not a matrix, network, or vector."))
    }
    
    # add
    add.row.labels <- character()
    add.col.labels <- character()
    if (add == TRUE) {
      # compile source and target row and column labels
      nr <- nrow(sources[[i]])  # save for later use
      source.row.labels <- rownames(sources[[i]])
      if (!sources.types[[i]] %in% c("matrix", "network")) {
        source.col.labels <- rownames(sources[[i]])
      } else {
        source.col.labels <- colnames(sources[[i]])
      }
      if (sources.types[[i]] %in% c("matrix", "network")) {
        if (is.null(source.row.labels)) {
          stop(paste0("The source at t = ", i, 
              " does not contain any row labels."))
        }
        if (is.null(source.col.labels)) {
          stop(paste0("The source at t = ", i, 
              " does not contain any column labels."))
        }
      }
      
      target.row.labels <- rownames(targets[[i]])
      if (!targets.types[[i]] %in% c("matrix", "network")) {
        target.col.labels <- rownames(targets[[i]])
      } else {
        target.col.labels <- colnames(targets[[i]])
      }
      if (is.null(target.row.labels)) {
        stop(paste0("The target at t = ", i, 
            " does not contain any row labels."))
      }
      if (targets.types[[i]] %in% c("matrix", "network")) {
        if (is.null(target.col.labels)) {
          stop(paste0("The target at t = ", i, 
              " does not contain any column labels."))
        }
      }
      
      add.row.indices <- which(!target.row.labels %in% source.row.labels)
      add.row.labels <- target.row.labels[add.row.indices]
      add.col.indices <- which(!target.col.labels %in% source.col.labels)
      add.col.labels <- target.col.labels[add.col.indices]
      
      # adjust rows
      if (length(add.row.indices) > 0) {
        for (j in 1:length(add.row.indices)) {
          insert <- rep(value, ncol(sources[[i]]))
          part1 <- sources[[i]][0:(add.row.indices[j] - 1), ]
          if (class(part1) != "matrix") {
            if (sources.types[[i]] == "matrix") {
              part1 <- matrix(part1, nrow = 1)
            } else {
              part1 <- matrix(part1, ncol = 1)
            }
          }
          rownames(part1) <- rownames(sources[[i]])[0:(add.row.indices[j] - 1)]
          if (add.row.indices[j] <= nrow(sources[[i]])) {
            part2 <- sources[[i]][add.row.indices[j]:nrow(sources[[i]]), ]
          } else {
            part2 <- matrix(ncol = ncol(sources[[i]]), nrow = 0)
          }
          if (class(part2) != "matrix") {
            part2 <- matrix(part2, nrow = 1)
          }
          if (nrow(part2) > 0) {
            rownames(part2) <- rownames(sources[[i]])[add.row.indices[j]:
                nrow(sources[[i]])]
            sources[[i]] <- rbind(part1, insert, part2)
          } else {
            sources[[i]] <- rbind(part1, insert)
          }
          rownames(sources[[i]])[add.row.indices[j]] <- add.row.labels[j]
          
          # adjust nodal attributes (in the one-mode case)
          if (sources.types[[i]] == "network" && sources.onemode[[i]] == TRUE) {
            for (k in 1:length(sources.attributes[[i]])) {
              at1 <- sources.attributes[[i]][[k]][0:(add.row.indices[j] - 1)]
              at2 <- sources.attributes[[i]][[k]][add.row.indices[j]:length(
                  sources.attributes[[i]][[k]])]
              if (sources.attribnames[[i]][k] == "vertex.names") {
                sources.attributes[[i]][[k]] <- c(at1, add.row.labels[j], at2)
              } else if (sources.attribnames[[i]][k] == "na") {
                sources.attributes[[i]][[k]] <- c(at1, TRUE, at2)
              } else {
                sources.attributes[[i]][[k]] <- c(at1, value, at2)
              }
            }
          }
        }
      }
      
      # adjust columns
      if (length(add.col.indices) > 0 && sources.types[[i]] %in% c("matrix", 
          "network")) {
        for (j in 1:length(add.col.indices)) {
          insert <- rep(value, nrow(sources[[i]]))
          part1 <- sources[[i]][, 0:(add.col.indices[j] - 1)]
          if (class(part1) != "matrix") {
            part1 <- matrix(part1, ncol = 1)
          }
          colnames(part1) <- colnames(sources[[i]])[0:(add.col.indices[j] - 1)]
          if (add.col.indices[j] <= ncol(sources[[i]])) {
            part2 <- sources[[i]][, add.col.indices[j]:ncol(sources[[i]])]
          } else {  # if last column, add empty column as second part
            part2 <- matrix(nrow = nrow(sources[[i]]), ncol = 0)
          }
          if (class(part2) != "matrix") {
            part2 <- matrix(part2, ncol = 1)
          }
          if (ncol(part2) > 0) {
            colnames(part2) <- colnames(sources[[i]])[add.col.indices[j]:
                ncol(sources[[i]])]
            sources[[i]] <- cbind(part1, insert, part2)
          } else {
            sources[[i]] <- cbind(part1, insert)
          }
          colnames(sources[[i]])[add.col.indices[j]] <- add.col.labels[j]
        }
      }
      
      # adjust nodal attributes for two-mode networks
      if (sources.types[[i]] == "network" && sources.onemode[[i]] == FALSE) {
        add.col.indices <- sapply(add.col.indices, function(x) x + nr)
        combined.indices <- c(add.row.indices, add.col.indices)
        for (j in 1:length(sources.attributes[[i]])) {
          if (length(combined.indices) > 0) {
            for (k in 1:length(combined.indices)) {
              at1 <- sources.attributes[[i]][[j]][0:(combined.indices[k] - 1)]
              at2 <- sources.attributes[[i]][[j]][combined.indices[k]:length(
                  sources.attributes[[i]][[j]])]
              if (sources.attribnames[[i]][j] == "vertex.names") {
                sources.attributes[[i]][[j]] <- c(at1, add.col.labels[j], at2)
              } else if (sources.attribnames[[i]][j] == "na") {
                sources.attributes[[i]][[j]] <- c(at1, TRUE, at2)
              } else {
                sources.attributes[[i]][[j]] <- c(at1, value, at2)
              }
            }
          }
        }
      }
    }
    
    removed.rows <- character()
    removed.columns <- character()
    if (remove == TRUE) {
      # compile source and target row and column labels
      nr <- nrow(sources[[i]])  # save for later use
      source.row.labels <- rownames(sources[[i]])
      if (!sources.types[[i]] %in% c("matrix", "network")) {
        source.col.labels <- rownames(sources[[i]])
      } else {
        source.col.labels <- colnames(sources[[i]])
      }
      if (sources.types[[i]] %in% c("matrix", "network")) {
        if (nr == 0) {
          stop(paste0("The source at t = ", i, " has no rows."))
        }
        if (is.null(source.row.labels)) {
          stop(paste0("The source at t = ", i, 
              " does not contain any row labels."))
        }
        if (is.null(source.col.labels)) {
          stop(paste0("The source at t = ", i, 
              " does not contain any column labels."))
        }
      }
      
      target.row.labels <- rownames(targets[[i]])
      if (!targets.types[[i]] %in% c("matrix", "network")) {
        target.col.labels <- rownames(targets[[i]])
      } else {
        target.col.labels <- colnames(targets[[i]])
      }
      if (targets.types[[i]] %in% c("matrix", "network")) {
        if (is.null(target.row.labels)) {
          stop(paste0("The target at t = ", i, 
              " does not contain any row labels."))
        }
        if (is.null(target.col.labels)) {
          stop(paste0("The target at t = ", i, 
              " does not contain any column labels."))
        }
      }
      
      # remove
      source.row.labels <- rownames(sources[[i]])
      source.col.labels <- colnames(sources[[i]])
      target.row.labels <- rownames(targets[[i]])
      target.col.labels <- colnames(targets[[i]])
      keep.row.indices <- which(source.row.labels %in% target.row.labels)
      if (sources.types[[i]] %in% c("matrix", "network") && 
          targets.types[[i]] %in% c("matrix", "network")) {
        keep.col.indices <- which(source.col.labels %in% target.col.labels)
      } else if (sources.types[[i]] %in% c("matrix", "network") 
          && !targets.types[[i]] %in% c("matrix", "network")) { 
        # target is a vector -> keep all columns of source if not onemode
        if (sources.onemode[[i]] == TRUE) {  # columns same as rows
          keep.col.indices <- keep.row.indices
        } else {
          keep.col.indices <- 1:ncol(sources[[i]])
        }
      } else {
        keep.col.indices <- 1
      }
      removed.rows <- which(!1:nrow(as.matrix(sources[[i]])) %in% 
          keep.row.indices)
      removed.columns <- which(!1:ncol(as.matrix(sources[[i]])) %in% 
          keep.col.indices)
      
      sources[[i]] <- as.matrix(sources[[i]][keep.row.indices, 
          keep.col.indices])
      if (sources.types[[i]] == "network") {
        if (sources.onemode[[i]] == TRUE) {
          for (j in 1:length(sources.attributes[[i]])) {
            sources.attributes[[i]][[j]] <- sources.attributes[[i]][[j]][
                keep.row.indices]
          }
        } else {
          keep.col.indices <- sapply(keep.col.indices, function(x) x + nr)
          combined.indices <- c(keep.row.indices, keep.col.indices)
          for (j in 1:length(sources.attributes[[i]])) {
            sources.attributes[[i]][[j]] <- sources.attributes[[i]][[j]][
                combined.indices]
          }
        }
      }
    }
    
    # sort source (and attributes) according to row and column names of target
#    if (length(sources.attributes[[i]]) > 0) {
#      for (j in 1:length(sources.attributes[[i]])) {
#        if (!is.null(sources.attributes[[i]][[j]]) && 
#            length(sources.attributes[[i]][[j]]) > 0) {
#          if (sources.onemode[[i]] == TRUE) {
#            names(sources.attributes[[i]][[j]]) <- rownames(sources[[i]])
#            sources.attributes[[i]][[j]] <- 
#                sources.attributes[[i]][[j]][rownames(sources[[i]])]
#          } else {
#            names(sources.attributes[[i]][[j]]) <- c(rownames(sources[[i]]), 
#                rownames(sources[[i]]))
#            sources.attributes[[i]][[j]] <- 
#                c(sources.attributes[[i]][[j]][rownames(sources[[i]])], 
#                sources.attributes[[i]][[j]][colnames(sources[[i]])])
#          }
#        }
#      }
#    }
#    
    if (sources.types[[i]] %in% c("matrix", "network") && 
        targets.types[[i]] %in% c("matrix", "network") && 
        nrow(sources[[i]]) == nrow(targets[[i]]) && 
        ncol(sources[[i]]) == ncol(targets[[i]])) {
      sources[[i]] <- sources[[i]][rownames(targets[[i]]), 
          colnames(targets[[i]])]
    } else if (sources.types[[i]] %in% c("matrix", "network") && 
        !targets.types[[i]] %in% c("matrix", "network") && 
        nrow(sources[[i]]) == nrow(targets[[i]])) {
      sources[[i]] <- sources[[i]][rownames(targets[[i]]), 
          rownames(targets[[i]])]
    } else if (length(sources[[i]]) == nrow(targets[[i]])) {
      # source is a vector, irrespective of the target
      sources[[i]] <- sources[[i]][rownames(targets[[i]]), ]
    } else if (add == FALSE && (nrow(sources[[i]]) < nrow(targets[[i]]) || 
        any(rownames(sources[[i]]) != rownames(targets[[i]])))) {
    }
    
    # convert back into network
    if (sources.types[[i]] == "network") {
      sources[[i]] <- network(sources[[i]], directed = sources.directed[[i]], 
          bipartite = !sources.onemode[[i]])
      for (j in 1:length(sources.attribnames[[i]])) {
        sources[[i]] <- set.vertex.attribute(sources[[i]], 
            sources.attribnames[[i]][j], sources.attributes[[i]][[j]])
      }
    }
    
    # convert vectors back from one-column matrices to vectors
    if (!sources.types[[i]] %in% c("matrix", "network") && 
        class(sources[[i]]) == "matrix" && ncol(sources[[i]]) == 1) {
      sources[[i]] <- sources[[i]][, 1]
    }
    
    # return added and removed labels instead of actual objects
    if (returnlabels == TRUE) {
      sources[[i]] <- list()
      sources[[i]]$removed.row <- removed.rows
      sources[[i]]$removed.col <- removed.columns
      sources[[i]]$added.row <- add.row.labels
      sources[[i]]$added.col <- add.col.labels
    }
  }
  
  if (sources.initialtype == "list") {
    return(sources)
  } else {
    return(sources[[1]])
  }
}


# handle missing data and absent nodes for multiple time points with lags
preprocess <- function(object, ..., lag = FALSE, covariate = FALSE, 
    memory = c("no", "autoregression", "stability", "innovation", "loss"), 
    na = NA, na.method = "fillmode", structzero = -9, 
    structzero.method = "remove", verbose = FALSE) {
  
  # save objects in a list
  l <- as.list(match.call())[-1]  # names of all objects including 'object'
  arg <- c("lag", "covariate", "memory", "na", "na.method", "structzero", 
      "structzero.method", "verbose")  # these args are not included as objects
  for (i in 1:length(arg)) {
    if (arg[i] %in% names(l)) {
      l <- l[-which(names(l) == arg[i])]
    }
  }
  for (i in 1:length(l)) {  # save the objects in list, not just the names
    l[[i]] <- eval(l[[i]], envir = .GlobalEnv)
    if (class(l[[i]]) != "list") {
      l[[i]] <- list(l[[i]])
    }
  }
  
  # check number of time steps; what is the length of the longest list in l?
  t <- 1
  for (i in 1:length(l)) {
    if (class(l[[i]]) == "list" && length(l[[i]]) > t) {
      t <- length(l[[i]])
    }
  }
  
  # make sure that all lists have the same length
  for (i in 1:length(l)) {
    if (length(l[[i]]) == 1 && t > 1) {
      for (j in 2:t) {
        l[[i]][[j]] <- l[[i]][[1]]
      }
    }
    if (length(l[[i]]) != t) {
      stop(paste("Object", i, "does not contain the right number of elements."))
    }
  }
  
  # infer labels if length corresponds to largest item otherwise in list
  largest.nr <- 0
  largest.row.labels <- character()
  largest.nc <- 0
  largest.col.labels <- character()
  for (i in 1:length(l)) {  # find largest items and save their labels
    for (j in 1:length(l[[i]])) {
      if (class(l[[i]][[j]]) == "matrix" || class(l[[i]][[j]]) == "network") {
        if (nrow(as.matrix(l[[i]][[j]])) > largest.nr && 
            !is.null(rownames(as.matrix(l[[i]][[j]])))) {
          largest.nr <- nrow(as.matrix(l[[i]][[j]]))
          largest.row.labels <- rownames(as.matrix(l[[i]][[j]]))
        }
        if (ncol(as.matrix(l[[i]][[j]])) > largest.nc && 
            !is.null(colnames(as.matrix(l[[i]][[j]])))) {
          largest.nc <- ncol(as.matrix(l[[i]][[j]]))
          largest.col.labels <- colnames(as.matrix(l[[i]][[j]]))
        }
      } else if (class(l[[i]][[j]]) == "data.frame") {
        if (nrow((l[[i]][[j]]) > largest.nr && 
            !is.null(rownames(l[[i]][[j]])))) {
          largest.nr <- nrow(l[[i]][[j]])
          largest.row.labels <- rownames(l[[i]][[j]])
        }
      } else if (is.vector(l[[i]][[j]]) && length(l[[i]][[j]]) > largest.nr && 
          !is.null(names(l[[i]][[j]]))) {
        largest.nr <- length(l[[i]][[j]])
        largest.row.labels <- names(l[[i]][[j]])
      }
    }
  }
  if (largest.nc == 0) {  # if col labels are never given, assume same as row l.
    largest.nc <- largest.nr
    largest.col.labels <- largest.row.labels
  }
  
  for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
      if (class(l[[i]][[j]]) == "matrix") {
        if (nrow(l[[i]][[j]]) == largest.nr && is.null(rownames(l[[i]][[j]]))) {
          rownames(l[[i]][[j]]) <- largest.row.labels
        }
        if (ncol(l[[i]][[j]]) == largest.nc && is.null(colnames(l[[i]][[j]]))) {
          colnames(l[[i]][[j]]) <- largest.col.labels
        }
      } else if (class(l[[i]][[j]]) == "network") {
        if (nrow(as.matrix(l[[i]][[j]])) == largest.nr && is.null(rownames(
            as.matrix(l[[i]][[j]])))) {
          if (!is.bipartite(l[[i]][[j]])) {
            l[[i]][[j]] <- set.vertex.attribute(l[[i]][[j]], "vertex.names", 
                largest.row.labels)
          } else {
            l[[i]][[j]] <- set.vertex.attribute(l[[i]][[j]], "vertex.names", 
                c(largest.row.labels, get.vertex.attribute(l[[i]][[j]], 
                "vertex.names")[(nrow(as.matrix(l[[i]][[j]]) + 1):ncol(
                as.matrix(l[[i]][[j]])))]))
          }
        }
        if (ncol(as.matrix(l[[i]][[j]])) == largest.nc && is.null(colnames(
            as.matrix(l[[i]][[j]]))) && is.bipartite(l[[i]][[j]])) {
          l[[i]][[j]] <- set.vertex.attribute(l[[i]][[j]], "vertex.names", 
          c(get.vertex.attribute(l[[i]][[j]], 
              "vertex.names")[1:nrow(as.matrix(l[[i]][[j]]))], 
              largest.col.labels))
        }
      } else if (class(l[[i]][[j]]) == "data.frame" && is.null(
          rownames(l[[i]][[j]])) && nrow(l[[i]][[j]]) == largest.nr) {
        rownames(l[[i]][[j]]) <- largest.row.labels
      } else if (is.vector(l[[i]][[j]]) && length(l[[i]][[j]]) == largest.nr && 
          is.null(names(l[[i]][[j]]))) {
        names(l[[i]][[j]]) <- largest.row.labels
      }
    }
  }
  
  # remove or replace NA values
  for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
      for (k in 1:length(na)) {
        if (verbose == TRUE) {
          l[[i]][[j]] <- handleMissings(l[[i]][[j]], na = na[k], 
              method = na.method)
        } else {
          l[[i]][[j]] <- suppressMessages(handleMissings(l[[i]][[j]], 
              na = na[k], method = na.method))
        }
      }
      for (k in 1:length(structzero)) {
        if (verbose == TRUE) {
          l[[i]][[j]] <- handleMissings(l[[i]][[j]], na = structzero[k], 
              method = structzero.method)
        } else {
          l[[i]][[j]] <- suppressMessages(handleMissings(l[[i]][[j]], 
              na = structzero[k], method = structzero.method))
        }
      }
    }
  }
  
  # adjust matrix dimensions to each other cross-sectionally (across list items)
  for (i in 1:t) {
    for (j in 1:length(l)) {
      for (k in 1:length(l)) {
        tryCatch({
          l[[j]][[i]] <- adjust(l[[j]][[i]], l[[k]][[i]], add = FALSE)
        }, error = function(e) {
          stop(paste0("The following error was encountered when executing ", 
              "the adjust function at time step ", i, ": ", e))
        })
      }
    }
  }
  
  # forward adjustment for lagged covariates
  forward <- l
  if (lag == TRUE && t > 1) {
    for (i in 1:length(forward)) {
      for (j in 1:(t - 1)) {
        forward[[i]][[j]] <- adjust(forward[[i]][[j]], forward[[i]][[j + 1]], 
            remove = TRUE, add = FALSE)
      }
      if (covariate == TRUE) {  # remove last time step
        forward[[i]] <- forward[[i]][-t]
      }  # result if covariate = TRUE and lag = TRUE
    }
  }
  # backward adjustment for dependent networks when lagged cov are present
  if (lag == TRUE && t > 1) {
    backward <- l
    for (i in 1:length(backward)) {
      for (j in 2:t) {
        backward[[i]][[j]] <- adjust(backward[[i]][[j]], forward[[i]][[j - 1]], 
            add = FALSE)
      }
      backward[[i]] <- backward[[i]][-1]  # remove first time step
    }  # result when covariate = FALSE and lag = TRUE
  }
  
  # create memory term
  if (memory[1] == "no") {
    # no memory term
  } else if (memory[1] %in% c("stability", "autoregression", "innovation", 
      "loss")) {
    if (lag == FALSE || covariate == FALSE) {
      stop(paste("Memory terms can only be created in conjunction with", 
          "'lag = TRUE' and 'covariate = TRUE'."))
    }
    if (t < 2) {
      stop("Only one time step. Memory term cannot be created.")
    }
    
    # compute memory terms by comparing forward with backward matrices
    memlist <- list()
    for (i in 1:length(forward[[1]])) {
      if (memory[1] == "autoregression") {
        memlist[[i]] <- as.matrix(forward[[1]][[i]])
      } else if (memory[1] == "stability") {
        memlist[[i]] <- as.matrix(forward[[1]][[i]])
        memlist[[i]][memlist[[i]] == 0] <- -1
      } else if (memory[1] == "innovation") {
        memlist[[i]] <- (as.matrix(forward[[1]][[i]]) * -1 + 1) * as.matrix(
            backward[[1]][[i]])
      } else if (memory[1] == "loss") {
        memlist[[i]] <- as.matrix(forward[[1]][[i]]) * (1 - as.matrix(
            backward[[1]][[i]]))
      } else {
        stop("'memory' argument not recognized.")
      }
    }
    
    # assign names
    names(memlist) <- paste0("t", 1:length(forward[[1]]), "-t", 
        2:(length(forward[[1]]) + 1))
  } else {
    stop("The 'memory' argument was not recognized.")
  }
  
  # determine return object
  if (covariate == FALSE && lag == TRUE) {
    obj <- backward[[1]]
  } else if (covariate == TRUE && lag == TRUE) {
    if (memory[1] == "no") {
      obj <- forward[[1]]
    } else {
      obj <- memlist
    }
  } else {
    obj <- l[[1]]
  }
  if (class(obj) == "list" && length(obj) == 1) {
    obj <- obj[[1]]
  }
  return(obj)
}

