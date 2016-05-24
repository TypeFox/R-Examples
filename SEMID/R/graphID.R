#' A helper function to validate input matrices.
#'
#' This helper function validates that the two input matrices, L and O, are of
#' the appropriate form to be interpreted by the other functions. In particular
#' they should be square matrices of 1's and 0's with all 0's along their
#' diagonals. We do not require O to be symmetric here.
#'
#' @param L See above description.
#' @param O See above description.
#'
#' @return This function has no return value.
validateMatrices = function(L, O) {
  if (!is.matrix(L) || !is.matrix(O)) {
    stop("L and O must be matrices.")
  } else if (length(unique(c(dim(L), dim(O)))) != 1) {
    stop("L and O must both be square matrices of the same dimensions.")
  }
  t1 = all(L %in% c(0, 1))
  t2 = all(O %in% c(0, 1))
  if (!t1 || !t2) {
    stop("L and O must contain only 1's and 0's.")
  } else if (any(diag(L) != 0) || any(diag(O) != 0)) {
    stop("L and O must have 0's along their diagonals.")
  }
}

#' Identifiability of linear structural equation models.
#'
#' This function checks global and generic identifiability of linear
#' structural equation models. For generic identifiability the function
#' checks a sufficient criterion as well as a necessary criterion but this
#' check may be inconclusive.
#'
#' @export
#'
#' @param L Adjacency matrix for the directed part of the path
#' diagram/mixed graph; an edge pointing from i to j is encoded as L[i,j]=1 and
#' the lack of an edge between i and j is encoded as L[i,j]=0. There should be
#' no directed self loops, i.e. no i such that L[i,i]=1.
#' @param O Adjacency matrix for the bidirected part of the path diagram/mixed
#' graph. Edges are encoded as for the L parameter. Again there should be no
#' self loops. Also this matrix will be coerced to be symmetric so it is only
#' necessary to specify an edge once, i.e. if O[i,j]=1 you may, but are not
#' required to, also have O[j,i]=1.
#' @param output.type A character string indicating whether output is
#' printed ('matrix'), saved to a file ('file'), or returned as a list
#' ('list') for further processing in R.
#' @param file.name A character string naming the output file.
#' @param decomp.if.acyclic A logical value indicating whether an input graph
#' that is acyclic is to be decomposed before applying identifiability criteria.
#' @param test.globalID A logical value indicating whether or not global
#' identifiability is checked.
#' @param test.genericID A logical value indicating whether or not a sufficient
#' condition for generic identifiability is checked.
#' @param test.nonID A logical value indicating whether or not a condition
#' implying generic non-identifiability is checked.
#'
#' @return
#'   A list or printed matrix indicating the identifiability status of the
#'   linear SEM given by the input graph.  Optionally the graph's
#'   components are listed.
#'
#'   With output.type = 'list', the function returns a list of components
#'   for the graph.  Each list entry is again a list that indicates first
#'   which nodes form the component and second whether the component forms
#'   a mixed graph that is acyclic.  The next entries in the list show
#'   HTC-identifiable nodes, meaning nodes v for which the coefficients for
#'   all the directed edges pointing to v can be identified using the
#'   methods from Foygel et al. (2012).  The HTC-identifiable nodes are
#'   listed in the order in which they are found by the recursive
#'   identification algorithm.  The last three list entries are
#'   logical values that indicate whether or not the graph component is
#'   generically identifiable, globally identifiable or not identifiable;
#'   compare Drton et al. (2011) and Foygel et al. (2012).  In the latter
#'   case the Jacobian of the parametrization does not have full rank.
#'
#'   With output.type = 'matrix', a summary of the above
#'   information is printed.
#'
#' @references Drton, M., Foygel, R., and Sullivant, S.  (2011) Global
#' identifiability of linear structural equation models. \emph{Ann. Statist.}
#' 39(2): 865-886.
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
#'
#' @examples
#' L = t(matrix(
#'   c(0, 1, 0, 0, 0,
#'     0, 0, 1, 0, 0,
#'     0, 0, 0, 1, 0,
#'     0, 0, 0, 0, 1,
#'     0, 0, 0, 0, 0), 5, 5))
#' O = t(matrix(
#'   c(0, 0, 1, 1, 0,
#'     0, 0, 0, 1, 1,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 0), 5, 5))
#' O=O+t(O)
#' graphID(L,O)
#'
#' ## Examples from Foygel, Draisma & Drton (2012)
#' demo(SEMID)
graphID <- function(L, O, output.type = 'matrix', file.name = NULL,
                    decomp.if.acyclic = TRUE, test.globalID = TRUE,
                    test.genericID = TRUE, test.nonID = TRUE) {
  if (!is.matrix(L) || !is.matrix(O)) {
    print('L and O must be two square matrices of the same size')
    return(NULL)
  } else {
    m = unique(c(dim(L),dim(O)))
    if (length(m) > 1) {
      print('L and O must be two square matrices of the same size')
      return(NULL)
    } else {
      if (!is.character(file.name) && output.type == 'file') {
        print('need to input a file name for output.type = \'file\'')
        return(NULL)
      } else {
        L <- (L != 0) ; diag(L) <- 0
        O <- O + t(O) ; O <- (O != 0) ; diag(O) <- 0
        Graph.output <- graphID.decompose(L, O, decomp.if.acyclic,
                                          test.globalID, test.genericID,
                                          test.nonID)
        if (all(c(test.globalID, test.genericID, test.nonID))) {
          graphID.output.alltests(Graph.output$Components, Graph.output$Decomp,
                                  output.type, file.name)
        } else {
          graphID.output.notalltests(Graph.output$Components,
                                     Graph.output$Decomp, output.type,
                                     file.name, test.globalID, test.genericID,
                                     test.nonID)
        }
      }
    }
  }
}

###
# Function to format the output when all tests have been run.
###
graphID.output.alltests <- function(Graph.components, Decomp, output.type,
                                    file.name) {
  if (output.type == 'list') {
    return(Graph.components)
  } else {
    Graph.output.matrix <- NULL
    for (comp in 1:length(Graph.components)) {
      if (Graph.components[[comp]]$acyclic) {
        if (Graph.components[[comp]]$GlobalID) {
          Comp.status <- 'Globally ID\'able'
        } else {
          if (Graph.components[[comp]]$GenericID) {
            Comp.status <- 'Generically ID\'able'
          } else {
            if (Graph.components[[comp]]$NonID) {
              Comp.status <- 'Generically non-ID\'able'
            } else {
              Comp.status <- 'HTC-inconclusive'
            }
          }
        }
      } else {
        if (Graph.components[[comp]]$GenericID) {
          Comp.status <- 'Generically ID\'able'
        } else {
          if (Graph.components[[comp]]$NonID) {
            Comp.status <- 'Generically non-ID\'able'
          } else {
            Comp.status <- 'HTC-inconclusive'
          }
        }
      }

      if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0) {
        Graph.components[[comp]]$HTC.ID.nodes <- 'none'
      }

      Comp.row <- c(comp,
                    paste(Graph.components[[comp]]$Nodes,collapse = ','),
                    paste(Graph.components[[comp]]$HTC.ID.nodes,collapse = ','),
                    Comp.status)
      Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
    }

    Graph.output.matrix <- rbind(c('Component', 'Nodes', 'HTC-ID\'able nodes',
                                   'Identifiability'), Graph.output.matrix)
    colnames(Graph.output.matrix) <- rep('',4)
    rownames(Graph.output.matrix) <- rep('',length(Graph.components) + 1)

    if (!Decomp) {
      Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
    }

    Max.string <- max(nchar(Graph.output.matrix))
    for (i in 1:dim(Graph.output.matrix)[1]) {
      for (j in 1:dim(Graph.output.matrix)[2]) {
        Graph.output.matrix[i,j] <-
          paste(c(Graph.output.matrix[i,j],
                  rep(' ', 5 + Max.string - nchar(Graph.output.matrix[i,j]))),
                collapse = '')
      }
    }

    if (output.type == 'file') {
      write.table(Graph.output.matrix, file = file.name, quote = FALSE)
    } else {
      print(Graph.output.matrix, quote = FALSE)
    }
  }
}

##
# Format output when not all tests have been run.
##
graphID.output.notalltests <- function(Graph.components, Decomp, output.type,
                                       file.name, test.globalID, test.genericID,
                                       test.nonID) {
  if (output.type == 'list') {
    return(Graph.components)
  } else {
    Graph.output.matrix <- NULL
    for (comp in 1:length(Graph.components)) {
      if (Graph.components[[comp]]$acyclic && test.globalID) {
        if (Graph.components[[comp]]$GlobalID) {
          Comp.status.GlobalID <- 'yes'
        } else {
          Comp.status.GlobalID <- 'no'
        }
      } else {
        Comp.status.GlobalID <- 'not tested'
      }
      if (test.genericID) {
        if (Graph.components[[comp]]$GenericID) {
          Comp.status.GenericID <- 'yes'
        } else {
          Comp.status.GenericID <- 'no'
        }
      } else {
        Comp.status.GenericID <- 'not tested'
      }
      if (test.nonID) {
        if (Graph.components[[comp]]$NonID) {
          Comp.status.NonID <- 'yes'
        } else {
          Comp.status.NonID <- 'no'
        }
      } else {
        Comp.status.NonID <- 'not tested'
      }

      if (length(Graph.components[[comp]]$HTC.ID.nodes) == 0) {
        if (test.genericID) {
          Graph.components[[comp]]$HTC.ID.nodes <- 'none'
        } else {
          Graph.components[[comp]]$HTC.ID.nodes <- 'not tested'
        }
      }

      Comp.row <- c(comp,
                    paste(Graph.components[[comp]]$Nodes,collapse = ','),
                    paste(Graph.components[[comp]]$HTC.ID.nodes,collapse = ','),
                    Comp.status.GlobalID, Comp.status.GenericID,
                    Comp.status.NonID)
      Graph.output.matrix <- rbind(Graph.output.matrix, Comp.row)
    }

    Graph.output.matrix <- rbind(c('Component', 'Nodes', 'HTC-ID\'able nodes',
                                   'Globally ID\'able?',
                                   'Generically ID\'able?',
                                   'Generically non-ID\'able?'),
                                 Graph.output.matrix)
    colnames(Graph.output.matrix) <- rep('', 6)
    rownames(Graph.output.matrix) <- rep('', length(Graph.components) + 1)

    if (!Decomp) {
      Graph.output.matrix <- Graph.output.matrix[, -(1:2)]
    }

    Max.string <- max(nchar(Graph.output.matrix))
    for (i in 1:dim(Graph.output.matrix)[1]) {
      for (j in 1:dim(Graph.output.matrix)[2]) {
        Graph.output.matrix[i,j] <-
          paste(c(Graph.output.matrix[i,j],
                  rep(' ', 3 + Max.string - nchar(Graph.output.matrix[i,j]))),
                collapse = '')
      }
    }


    if (output.type == 'file') {
      write.table(Graph.output.matrix, file = file.name, quote = FALSE)
    } else {
      print(Graph.output.matrix, quote = FALSE)
    }
  }
}


#' Determine generic identifiability by Tian Decomposition and HTC
#'
#' Split a graph into mixed Tiancomponents and solve each separately
#' using the HTC.
#'
#' @inheritParams graphID
#'
#' @return A list with two named components:
#'
#'   1. Components - a list of lists. Each list represents one mixed Tian component
#'       of the graph. Each list contains named components corresponding to which
#'       nodes are in the component and results of various tests of
#'       identifiability on the component (see the parameter descriptions).
#'
#'  2. Decomp - true if a decomposition occured, false if not.
graphID.decompose <- function(L, O, decomp.if.acyclic = TRUE,
                              test.globalID = TRUE, test.genericID = TRUE,
                              test.nonID = TRUE) {
  m <- nrow(L)
  L <- (L != 0)
  diag(L) <- 0
  O <- O + t(O)
  O <- (O != 0)
  diag(O) <- 0

  Decomp <- FALSE
  if (decomp.if.acyclic) {
    test.acyclic <- TRUE
    nodes.acyclic <- 1:m
    while (test.acyclic && !Decomp) {
      if (any(colSums(L[nodes.acyclic, nodes.acyclic, drop = FALSE]) == 0)) {
        nodes.acyclic <- nodes.acyclic[-which(colSums(L[nodes.acyclic,
                                                        nodes.acyclic,
                                                        drop = FALSE]) == 0)]
        if (length(nodes.acyclic) == 0) {
          Decomp <- TRUE
        }
      } else {
        test.acyclic <- FALSE
      }
    }
  }

  if (Decomp) {
    Bidir.comp <- diag(m) + O
    for (i in 1:(m - 1)) {
      Bidir.comp <- (Bidir.comp %*% (diag(m) + O) > 0) + 0
    }
    Components <- list()
    V <- 1:m
    num.comp = 0
    while (length(V) > 0) {
      num.comp <- num.comp + 1
      i <- min(V)
      Components[[num.comp]] <- which(Bidir.comp[i,] == 1)
      V <- setdiff(V,Components[[num.comp]])
    }
  } else {
    Components <- list()
    Components[[1]] <- 1:m
    num.comp <- 1
  }

  for (comp in 1:num.comp) {
    Component <- Components[[comp]]
    Component.parents <-
      sort(setdiff(which(rowSums(L[,Component,drop = FALSE])>0),Component))
    if (length(Components[[comp]]) == 1) {
      Components[[comp]] <- list()
      Components[[comp]]$Nodes <- Component
      Components[[comp]]$acyclic <- TRUE
      if (test.genericID) {
        Components[[comp]]$HTC.ID.nodes <- Component
        Components[[comp]]$GenericID <- TRUE
      }
      if (test.globalID) {
        Components[[comp]]$GlobalID <- TRUE
      }
      if (test.nonID) {
        Components[[comp]]$NonID <- FALSE
      }
    } else {
      m1 <- length(Component) ; m2 <- length(Component.parents)
      L.Component <- cbind( L[c(Component,Component.parents),Component],
                            matrix(0, m1 + m2, m2) )
      O.Component <- cbind( rbind( O[Component,Component], matrix(0, m2, m1) ),
                            matrix(0, m1 + m2, m2) )
      Components[[comp]] <- c(list(Nodes = Component),
                              graphID.main(L.Component, O.Component,
                                           test.globalID, test.genericID,
                                           test.nonID))
      Components[[comp]]$HTC.ID.nodes <-
        Component[intersect(Components[[comp]]$HTC.ID.nodes,
                            1:length(Component))]
    }
  }
  Graph.output = list()
  Graph.output$Components <- Components
  Graph.output$Decomp <- Decomp
  return(Graph.output)
}

#' Helper function to handle a graph component.
#'
#' Calls the other functions that determine identifiability status.
#'
#' @inheritParams graphID
#'
#' @return A list containing named components of the results of various tests
#' desired based on the input parameters.
graphID.main <- function(L, O, test.globalID = TRUE, test.genericID = TRUE,
                         test.nonID = TRUE) {
  m <- nrow(L)
  L <- (L != 0)
  diag(L) <- 0
  O <- O + t(O)
  O <- (O != 0)
  diag(O) <- 0

  ILinv <- diag(m)
  for (i in 1:m) {
    ILinv <- 0 + (diag(m) + ILinv %*% L > 0)
  }

  Output <- list()

  Output$acyclic <- (max(ILinv + t(ILinv) - diag(m)) == 1)

  if (test.genericID) {
    Output$HTC.ID.nodes <- graphID.htcID(L,O)

    if (length(Output$HTC.ID.nodes) < m) {
      Output$GenericID <- FALSE
      if (Output$acyclic && test.globalID) {
        Output$GlobalID <- FALSE
      }
      if (test.nonID) {
        Output$NonID <- graphID.nonHtcID(L,O)
      }
    } else {
      Output$GenericID <- TRUE
      if (Output$acyclic && test.globalID) {
        Output$GlobalID <- graphID.globalID(L,O)
      }
      if (test.nonID) {
        Output$NonID <- FALSE
      }
    }
  } else {
    if (test.nonID) {
      Output$NonID <- graphID.nonHtcID(L,O)
      if (!Output$NonID) {
        if (Output$acyclic && test.globalID) {
          Output$GlobalID <- graphID.globalID(L,O)
        }
      }
    } else {
      if (Output$acyclic && test.globalID) {
        Output$GlobalID <- graphID.globalID(L,O)
      }
    }
  }
  return(Output)
}

#' Check for global identifiability of a mixed graph.
#'
#' Checks for the global identifiability of a mixed graph using techniques
#' presented in Drton, Foygel, Sullivant (2011).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return TRUE if the graph was globally identifiable, FALSE otherwise.
#'
#' @references
#' Drton, Mathias; Foygel, Rina; Sullivant, Seth. Global identifiability of
#' linear structural equation models. \emph{Ann. Statist.}  39 (2011), no. 2,
#' 865--886.
graphID.globalID <- function(L, O) {
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  if (!is.dag(graph.adjacency(L))) {
    return(F)
  }

  Global.ID = TRUE
  i <- 0
  while (Global.ID == 1 && i < m) {
    i <- i + 1
    S <- 1:m
    change <- 1
    while (change == 1) {
      change <- 0
      S.old <- S
      for (s in setdiff(S.old, i)) {
        if (graph.maxflow(graph.adjacency(L), source = s,
                          target = i)$value == 0) {
          S <- setdiff(S, s)
          change <- 1
        }
      }
      S.old <- S
      for (s in setdiff(S.old, i)) {
        if (graph.maxflow(graph.adjacency(O, mode = 'undirected'),
                          source = s, target = i)$value == 0) {
          S <- setdiff(S, s)
          change <- 1
        }
      }
    }
    if (length(S) > 1) {
      Global.ID <- FALSE
    }
  }

  return(Global.ID)
}

#' Check for generic infinite-to-one via the half-trek criterion.
#'
#' Checks if a mixed graph is infinite-to-one using the half-trek criterion
#' presented by Foygel, Draisma, and Drton (2012).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return TRUE if the graph could be determined to be generically
#' non-identifiable, FALSE if this test was inconclusive.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
graphID.nonHtcID <- function(L, O) {
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  # 1 & 2 = source & target
  # 2 + {1,...,N} = L{i_n,j_n} for the n-th nonsibling pair, n=1,...,N
  # (2+N) + 2*m^2 = 2 copies of R_i(j) -- in copy & out copy
  #         where (2+N) + (i-1)*m + j    = R_i(j) in copy
  #         & (2+N+m^2) + (i-1)*m + j    = R_i(j) out copy

  nonsibs <- NULL
  N <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      if (O[i,j] == 0) {
        N <- N + 1
        nonsibs <- rbind(nonsibs, c(i, j))
      }
    }
  }

  Cap.matrix <- matrix(0, 2*m^2 + N + 2, 2*m^2 + N + 2)
  if (N != 0) {
    Cap.matrix[1, 2 + (1:N)] <- 1 # edges from source to L{i,j} for each
                                  # nonsibling pair
    for (n in 1:N) {  #{i,j} = nonsibs[n,1:2]
      # edge from L{i,j} to R_i(j), and to R_i(k)-in for all siblings k of
      # node j
      Cap.matrix[2 + n, 2 + N +
                   (nonsibs[n,1] - 1)*m +
                   c(nonsibs[n,2], which(O[nonsibs[n,2], ] == 1))] <- 1
      # edge from L{i,j} to R_j(i), and to R_j(i)-in for all siblings k of
      # node i
      Cap.matrix[2 + n, 2 + N +
                   (nonsibs[n,2] - 1)*m +
                   c(nonsibs[n,1], which(O[nonsibs[n,1], ] == 1))] <- 1
    }
  }
  for (i in 1:m) {
    # edge from R_i(j)-out to target when j is a parent of i
    Cap.matrix[2 + N + m^2 + (i - 1)*m + which(L[, i] == 1), 2] <- 1
    for (j in 1:m) {
      # edge from R_i(j)-in to R_i(j)-out
      Cap.matrix[2 + N + (i - 1)*m + j,
                 2 + N + m^2 + (i - 1)*m + j] <- 1
      # edge from R_i(j)-out to R_i(k)-in where j->k is a directed edge
      Cap.matrix[2 + N + m^2 + (i - 1)*m + j,
                 2 + N + (i - 1)*m + which(L[j, ] == 1)] <- 1
    }
  }

  HTC.nonID <-
    graph.maxflow(graph.adjacency(Cap.matrix), source = 1, target = 2)$value
  return(HTC.nonID < sum(L))
}

#' Determine generic identifiability of a mixed graph.
#'
#' If directed part of input graph is cyclic then will check for generic
#' identifiability using the half-trek criterion. Otherwise will use the a
#' slightly stronger version of the half-trek criterion using ancestor
#' decompositions.
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of nodes that could be determined to be generically
#' identifiable.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
#'
#' @references
#' {Drton}, M. and {Weihs}, L. (2015) Generic Identifiability of Linear
#' Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992
graphID.genericID <- function(L, O) {
  if (is.dag(graph.adjacency(L, mode = "directed"))) {
    return(graphID.ancestralID(L, O))
  } else {
    return(graphID.htcID(L, O))
  }
}

#' Determines if a mixed graph is HTC-identifiable.
#'
#' Uses the half-trek criterion of Foygel, Draisma, and Drton (2015) to check
#' if an input mixed graph is generically identifiable.
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of HTC-identifiable nodes.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
graphID.htcID <- function(L, O) {
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  # 1 & 2 = source & target
  # 2 + {1,...,m} = L(i) for i=1,...,m
  # 2+m + {1,...,m} = R(i)-in for i=1,...,m
  # 2+2*m + {1,...,m} = R(i)-out for i=1,...,m

  Cap.matrix.init <- matrix(0, 2 + 3*m, 2 + 3*m)
  for (i in 1:m) {
    # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix.init[2 + i, 2 + m + c(i, which(O[i,] == 1))] <- 1
    # edge from R(i)-in to R(i)-out
    Cap.matrix.init[2 + m + i, 2 + 2*m + i] <- 1
    # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix.init[2 + 2*m + i, 2 + m + which(L[i,] == 1)] <- 1
  }

  # when testing if a set A satisfies the HTC with respect to a node i,
  #    need to add (1) edge from source to L(j) for all j in A
  #            and (2) edge from R(j)-out to target for all parents j of i

  Dependence.matrix <- O + diag(m)
  for (i in 1:m) {
    Dependence.matrix <- (Dependence.matrix + Dependence.matrix %*% L > 0)
  }

  Solved.nodes <- rep(0,m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1

  while (change == 1) {
    change <- 0

    for (i in which(Solved.nodes == 0)) {
      A <- setdiff(c(which(Solved.nodes > 0),
                     which(Dependence.matrix[i, ] == 0)),
                   c(i, which(O[i,] == 1)))

      Cap.matrix <- Cap.matrix.init

      Cap.matrix[1, 2 + A] <- 1
      Cap.matrix[2 + 2*m + which(L[,i] == 1), 2] <- 1

      flow <-
        graph.maxflow(graph.adjacency(Cap.matrix), source = 1, target = 2)$value

      if (flow == sum(L[,i])) {
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
      }
    }
  }

  if (all(Solved.nodes == 0)) {
    Solved.nodes <- NULL
  } else {
    Solved.nodes <- order(Solved.nodes)[(1 + sum(Solved.nodes == 0)):m]
  }
  return(Solved.nodes)
}

#' Get ancestors of nodes in a graph.
#'
#' Get the ancestors of a collection of nodes in a graph g, the ancestors DO
#' include the the nodes themselves.
#'
#' @param g the graph (as an igraph).
#' @param nodes the nodes in the graph of which to get the ancestors.
#'
#' @return a sorted vector of all ancestor nodes.
ancestors <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "in", unreachable = F)$order,
                  na.last = NA))
  #sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode="in"))))
}

#' Get parents of nodes in a graph.
#'
#' Get the parents of a collection of nodes in a graph g, the parents DO include
#' the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the parents.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all parent nodes.
parents <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "in"))))
}

#' Get siblings of nodes in a graph.
#'
#' Get the siblings of a collection of nodes in a graph g, the siblings DO
#' include the input nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the siblings.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all siblings of nodes.
siblings <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  sort(unique(unlist(neighborhood(g, 1, nodes = nodes, mode = "all"))))
}

#' Getdescendants of nodes in a graph.
#'
#' Gets the descendants of a collection of nodes in a graph (all nodes that can
#' be reached by following directed edges from those nodes). Descendants DO
#' include the nodes themselves.
#'
#' @param nodes the nodes in the graph of which to get the descendants.
#' @inheritParams ancestors
#'
#' @return a sorted vector of all descendants of nodes.
descendants <- function(g, nodes) {
  if (vcount(g) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  as.numeric(sort(graph.bfs(g, nodes, neimode = "out", unreachable = F)$order,
                  na.last = NA))
  #sort(unique(unlist(neighborhood(g, vcount(g), nodes=nodes, mode="out"))))
}

#' Get all HTR nodes from a set of nodes in a graph.
#'
#' Gets all vertices in a graph that are half-trek reachable from a set of
#' nodes.
#' WARNING: Often the half-trek reachable nodes from a vertex v are defined to
#' not include the vertex v or its siblings. We DO NOT follow this convention,
#' the returned set will include input nodes and their siblings.
#'
#' @inheritParams getMixedCompForNode
#' @param nodes the nodes in the graph of which to get the HTR nodes.
#'
#' @return a sorted list of all half-trek reachable nodes.
htr <- function(dG, bG, nodes) {
  if (!is.directed(dG) || is.directed(bG) || vcount(dG) != vcount(bG)) {
    stop("dG is undirected or bG is directed.")
  }
  if (vcount(dG) == 0 || length(nodes) == 0) {
    return(numeric(0))
  }
  return(descendants(dG, siblings(bG, nodes)))
}

#' Get the mixed component of a node in a mixed subgraph.
#'
#' For an input mixed graph H and set of nodes A, let GA be the subgraph of
#' H on the nodes A. This function returns the mixed component of GA containing
#' a specified node.
#'
#' @param dG a directed graph representing the directed part of the mixed graph.
#' @param bG an undirected graph representing the undirected part of the mixed
#'        graph.
#' @param subNodes an ancestral set of nodes in the mixed graph, this set should
#'        include the node for which the mixed component sould be found.
#' @param node the node for which the mixed component is found.
#' @return a list with two named elements:
#'          biNodes - the nodes of the mixed graph in the biDirected component
#'                    containing nodeName w.r.t the ancestral set of nodes
#'          inNodes - the nodes in the graph which are not part of biNodes
#'                    but which are a parent of some node in biNodes.
getMixedCompForNode <- function(dG, bG, subNodes, node) {
  VdG = V(dG)
  VbG = V(bG)
  m = vcount(dG)
  if (is.null(VdG$names) || is.null(VbG$names) ||
      any(VdG$names != 1:m) || any(VbG$names != 1:m)) {
    stop(paste("Input graphs to getMixedCompForNode must have vertices named",
               "1:m in order."))
  }

  bidirectedComp =
    as.numeric(sort(graph.bfs(bG, root = node, restricted = subNodes - 1,
                              neimode = "total", unreachable = F)$order))
  incomingNodes =
    intersect(setdiff(parents(dG, bidirectedComp), bidirectedComp),
              subNodes)
  return(list(biNodes = bidirectedComp, inNodes = incomingNodes))
}


#' Size of largest HT system Y satisfying the HTC for a node v except perhaps
#' having |parents(v)| < |Y|.
#'
#' For an input mixed graph H, constructs the Gflow graph as described in Foygel
#' et al. (2012) for a subgraph G of H. A max flow algorithm is then run on
#' Gflow to determine the largest half-trek system in G to a particular node's
#' parents given a set of allowed nodes. Here G should consist of a bidirected
#' part and nodes which are not in the bidirected part but are a parent of some
#' node in the bidirected part. G should contain the node for which to compute
#' the max flow.
#'
#' @param allowedNodes the set of allowed nodes.
#' @param biNodes the set of nodes in the subgraph G which are part of the
#'        bidirected part.
#' @param inNodes the nodes of the subgraph G which are not in the bidirected
#'        part but are a parent of some node in the bidirected component.
#' @param node the node (as an integer) for which the maxflow the largest half
#' trek system
#' @inheritParams graphID
#'
#' @return See title.
#'
#' @references
#' Foygel, R., Draisma, J., and Drton, M.  (2012) Half-trek criterion for
#' generic identifiability of linear structural equation models.
#' \emph{Ann. Statist.} 40(3): 1682-1713.
getMaxFlow <- function(L, O, allowedNodes, biNodes, inNodes, node) {
  if (!(node %in% biNodes) || any(O[biNodes, inNodes] != 0)) {
    stop(paste("When getting max flow either some in-nodes were connected to",
               "some bi-nodes or node was not in biNodes."))
  }
  if (length(intersect(allowedNodes, c(node, which(O[node,] == 1)))) != 0) {
    stop("Allowed nodes contained siblings of input node or the node itself.")
  }
  allowedNodes = union(allowedNodes, inNodes)
  m = length(biNodes) + length(inNodes)

  if (m == 1) {
    return(0)
  }

  oldNumToNewNum = numeric(m)
  oldNumToNewNum[biNodes] = 1:length(biNodes)
  oldNumToNewNum[inNodes] = (length(biNodes) + 1):m

  nodesToUse = c(biNodes, inNodes)
  allowedNodes = intersect(allowedNodes, nodesToUse)
  L[c(inNodes, biNodes), inNodes] = 0
  O[inNodes, inNodes] = 0
  L = L[nodesToUse, nodesToUse]
  O = O[nodesToUse, nodesToUse]

  # 1 & 2 = source & target
  # 2 + {1,...,m} = L(i) for i=1,...,m
  # 2 + m + {1,...,m} = R(i)-in for i=1,...,m
  # 2 + 2*m + {1,...,m} = R(i)-out for i=1,...,m

  Cap.matrix <- matrix(0, 2 + 3*m, 2 + 3*m)

  for (i in 1:m) {
    # edge from L(i) to R(i)-in, and to R(j)-in for all siblings j of i
    Cap.matrix[2 + i, 2 + m + c(i, which(O[i,] == 1))] <- 1
    # edge from R(i)-in to R(i)-out
    Cap.matrix[2 + m + i, 2 + 2*m + i] <- 1
    # edge from R(i)-out to R(j)-in for all directed edges i->j
    Cap.matrix[2 + 2*m + i, 2 + m + which(L[i,] == 1)] <- 1
  }

  allowedNodes = oldNumToNewNum[allowedNodes]
  node = oldNumToNewNum[node]
  Cap.matrix[1, 2 + allowedNodes] = 1
  Cap.matrix[2 + 2*m + which(L[,node] == 1), 2] = 1

  return(graph.maxflow(graph.adjacency(Cap.matrix),
                       source = 1, target = 2)$value)
}

#' Determine generic identifiability of an acyclic mixed graph using ancestral
#' decomposition.
#'
#' For an input, acyclic, mixed graph attempts to determine if the graph is
#' generically identifiable using decomposition by ancestral subsets. See
#' algorithm 1 of Drton and Weihs (2015).
#'
#' @export
#'
#' @inheritParams graphID
#'
#' @return The vector of nodes that could be determined to be generically
#' identifiable using the above algorithm.
#'
#' @references
#' {Drton}, M. and {Weihs}, L. (2015) Generic Identifiability of Linear
#' Structural Equation Models by Ancestor Decomposition. arXiv 1504.02992
graphID.ancestralID <- function(L, O) {
  m <- nrow(L)
  validateMatrices(L, O)
  O <- 1 * ((O + t(O)) != 0)

  dG = graph.adjacency(L)
  newOrder = as.numeric(topological.sort(dG))
  L = L[newOrder, newOrder]
  O = O[newOrder, newOrder]

  dG = graph.adjacency(L)
  bG = graph.adjacency(O)
  V(dG)$names = 1:m
  V(bG)$names = 1:m

  # Generates a list where, for each node v, we have a vector
  # corresponding to all the nodes that could ever be in a
  # half-trek system for v
  halfTrekSources = vector("list", length = m)
  for (i in 1:m) {
    halfTrekSources[[i]] = siblings(bG, ancestors(dG, i))
  }

  # A matrix determining which nodes are half-trek reachable from each node
  Dependence.matrix <- O + diag(m)
  for (i in 1:m) {
    Dependence.matrix <- ((Dependence.matrix + Dependence.matrix %*% L) > 0)
  }

  Solved.nodes <- rep(0, m)
  Solved.nodes[which(colSums(L) == 0)] <- 1 # nodes with no parents
  change <- 1
  count <- 1
  while (change == 1) {
    change <- 0

    Unsolved.nodes <- which(Solved.nodes == 0)
    for (i in Unsolved.nodes) {
      # A <- (Solved Nodes 'union' nodes not htr from i) \ ({i} 'union' sibs(i))
      A = setdiff(c(which(Solved.nodes > 0), which(Dependence.matrix[i,] == 0)),
                  c(i, which(O[i,] == 1)))
      # A <- A intersect (nodes that can ever be in a HT system for i)
      A = intersect(A, halfTrekSources[[i]])

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, c(i,A)), i)
      flow = getMaxFlow(L, O, A,
                        mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])) {
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }

      mixedCompList = getMixedCompForNode(dG, bG, ancestors(dG, i), i)
      A = intersect(A, unlist(mixedCompList))
      flow = getMaxFlow(L, O, A,
                        mixedCompList$biNodes, mixedCompList$inNodes, i)
      if (flow == sum(L[,i])) {
        change <- 1
        count <- count + 1
        Solved.nodes[i] <- count
        next
      }
    }
  }

  if (all(Solved.nodes == 0)) {
    Solved.nodes <- NULL
  } else {
    Solved.nodes <- order(Solved.nodes)[(1 + sum(Solved.nodes == 0)):m]
  }

  return(newOrder[Solved.nodes])
}
