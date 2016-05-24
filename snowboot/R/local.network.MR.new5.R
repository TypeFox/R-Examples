#' Construct Artificial Networks
#'
#' This function constructs an artificial network from a given distribution.
#' Only 11 distributions are available.
#'
#' @param n The number of nodes in the desired network.
#' @param distrib An atomic character representing the desired degree
#' distribution. User may choose from 11 available distributions: "fixed",
#' "pois", "ztpois", "geom", "nbinom", "ztgeom", "poly.log", "logarithmic",
#' "power.law", "full" (fully connected), or "none" (no element connected).
#' @param param The distribution parameters. If the function is "fixed",
#' \code{param} is a vector of degrees.
#' @param degree An optional vector of degrees that must be of length \code{n}.
#' The default is \code{degree = NULL}.
#' @param take.p A number between 0 and 1 representing the proportion to take
#'  for elimiation with each iteration.
#' @return A list consisting of
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{degree.left}{A vector of length \code{n} that should be all zeroes.}
#'    \item{n}{The network order. The order for every network is 2000.}
#' @export
#' @examples
#' a <- local.network.MR.new5(1000,"poly.log",c(2,13))

###################### NETWORK CONSTRUCTION MAIN FUNCTION ######################
local.network.MR.new5 <- function(n, distrib, param = NULL, degree = NULL,
                                  take.p = 0.05) {
  # param is not a list We develop the network construction by sampling in two stages.
  # We use the frequency tables of the vertices' degree
  # and arch weights (that are the product of the vertices's degree they connect)

  # network creates the network based in the edge distribution
  # and the number of nodes n number of
  # individuals (susceptible and infective).
  #  distrib is the degree distribution.
  # param is the distribution parameter (if the function is "fixed" it is
  # a vector of degrees)

  # In this interative algorithm we eliminate the archs that
  # cannot be longer be selected to construct a simple network In
  # order to optimize the algorithm we select more than one arch at each step

  if (!is.null(distrib) && is.null(degree)) {
    degree <- sdegree(n, distrib, param)
  } else {
    if (length(degree) != n)
      stop("degree has to have length ", n, "\n")
  }
  id <- 1:n
  edges <- matrix(NA, ceiling(sum(degree)/2), 2)  # I want to reserve the memory
  # for this variable using the
  # maximum number of edges
  # edges1<-edges2<-rep(NA,ceiling(sum(degree)/2))
  # of edges
  flag <- TRUE
  degree.left <- degree
  edge.row <- 1
  count <- 1

  while (sum(degree.left > 0) >= 2 & flag) {
    ## *********************************************************************##
    continue = TRUE
    eff.nodes <- (1:n)[degree.left > 0]
    # eff.degree<-degree.left[degree.left>0] it is the same as degree.left[eff.nodes]
    # as noted next. Now, if the network has more than 1e4,
    # I want to consider only cutnet nodes that I select by weight (degree)
    # I cut the number of nodes here
    m <- length(eff.nodes)
    # pl<-m*(m-1)/2 #number of elements in products ****************************
    # a<-Sys.time()
    tab.degree <- table(degree.left[degree.left > 0])  # the frequency of degree.left
    # that are greater than zero
    vals <- as.numeric(names(tab.degree))
    tab.prod.deg <- unlist(sapply(X = 1:(length(tab.degree)), FUN = table.mult.degree,
                                  vector = as.numeric(tab.degree),
                                  m = length(tab.degree), freq = TRUE))
    vals.prod.deg <- unlist(sapply(X = 1:(length(tab.degree)),
                                   FUN = table.mult.degree, vector = vals,
                                   m = length(tab.degree), freq = FALSE))
    take <- max(ceiling(length(tab.prod.deg) * take.p), 1)

    first.sample <- resample(1:length(tab.prod.deg), take,
                             prob = (tab.prod.deg * vals.prod.deg), rep = TRUE)
    # Here I sample the arch group
    # I want to allow to sample from the same group, but if in the second step
    # I repeat a vertex, I remove the rest of sampled groups
    # and repeated vertices ahora tengo que poder traducir el sitio muestreado
    # y el grupo de arcos con uno y otro grados
    degree1 <- vals[id1 <- cut(first.sample, breaks = br <- cumsum(c(0, length(vals):1)),
                               label = FALSE, include.lowest = TRUE)]
    degree2 <- vals[id2 <- id1 - 1 + first.sample - br[id1]]
    sam <- rep(NA, 2 * length(degree1))
    dd <- c(degree1, degree2)
    tdd <- table(dd)
    vtdd <- as.numeric(names(tdd))
    for (ss in 1:length(tdd)) {
      sam[dd == vtdd[ss]] <- resample(id[degree.left == vtdd[ss]], tdd[ss], rep = TRUE)
    }
    tt <- which(
      duplicated(as.vector(rbind(sam1 <- sam[1:length(degree1)],
                                 sam2 <- sam[(length(degree1) + 1):(2 * length(degree1))]))))
    if (length(tt) > 0) {
      tt <- min(tt)
    } else {
      tt <- NULL
    }
    # now I remove all the edges from tt and up
    if (!is.null(tt)) {
      sam1 <- sam1[-c(ceiling(tt/2):length(degree1))]
      sam2 <- sam2[-c(ceiling(tt/2):length(degree1))]
    }
    if (length(sam1) > 0) {
      new.node1 <- new.node2 <- rep(NA, length(sam1))
      new.node1[sam1 <= sam2] <- sam1[sam1 <= sam2]
      new.node1[sam1 > sam2] <- sam2[sam1 > sam2]
      new.node2[sam1 <= sam2] <- sam2[sam1 <= sam2]
      new.node2[sam1 > sam2] <- sam1[sam1 > sam2]
      # Now new.node1 is always smaller than new.node2
    } else {
      continue = FALSE
      cat("entra \n")
    }

    if (continue)
    {
      # if at least one edge candidate
      new.edge <- cbind(new.node1, new.node2)
      # *******************## avoid to repeat new edges when whe have more to
      # compare #edge.row ==1 for the first pair of nodes connected and >1 for the rest
      if (edge.row > 1 & length(new.node1) > 0)
      {
        a <- which(is.element(edges[1:edge.row, 1], new.edge[, 1]))
        # since new.edges sorted, a's who interest us are in the first column
        if (length(a) > 0) {
          b <- is.element(paste(new.edge[, 1], new.edge[, 2]), paste(edges[a, 1],
                                                                     edges[a, 2]))  #always in the second column
          if (any(b)) {
            # at least one repeated at leat one not repeated (we want to keep that one)
            if (any(!b)) {
              new.edge <- new.edge[!b, ]  #I eliminate the repeated
              if (is.vector(new.edge))
                new.edge <- t(new.edge)
            } else {
              new.edge <- NULL
            }  #all are repeated
          }
        }
      }
      ## ***************************************************************##
      if (!is.null(new.edge)) {
        edges[edge.row:(edge.row + dim(new.edge)[1] - 1), ] <- new.edge
        edge.row <- edge.row + dim(new.edge)[1]
        degree.left[new.edge[, 1]] <- degree.left[new.edge[, 1]] - 1
        degree.left[new.edge[, 2]] <- degree.left[new.edge[, 2]] - 1  #I can decrease only 1 in both cases because no nodes are repeated
      } else {
        if (m < 10)
          count <- count + 1
      }  #few nodes left to connect but unable to find edges that are not repeated
      if (count > 10)
        flag <- FALSE
      if (any(degree.left < 0)) {
        cat("raro\n")
        browser()
      }
    }  #continue
  }  #while

  edges <- edges[!is.na(edges[, 1]), ]  #remove NA's I put at the beginning
  edges <- order.edges(edges, ord.col = FALSE)
  list(edges = edges, degree = degree, degree.left = degree.left, n = n)
}  #Main function end
