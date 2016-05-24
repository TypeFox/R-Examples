# S3method for "tocher"
tocher <-
function(d, algorithm = c("original", "sequential")) UseMethod("tocher")

# ---------------------------------
# method for object of class dist
tocher.dist <-
function(d, algorithm = c("original", "sequential"))
{
   if (!inherits(d, "dist"))
      stop("'d' must be an object of class 'dist'!")
   algorithm <- match.arg(algorithm)
   d <- as.matrix(d)
   n <- nrow(d)

   # object labels
   if (is.null(dimnames(d))) {
       lab <- as.character(1:n)
       } else {
       lab <- colnames(d)
   }
   dimnames(d) <- list(lab, lab)

   # aux function to find the two closest objects
   fun.min <- function(mat)
   {
      n <- ncol(mat)
      v1 <- v2 <- NULL
      aux <- data.frame(v1 = rep(colnames(mat), each = n),
         v2 = rep(colnames(mat), times = n),
         val = as.vector(mat))
      aux2 <- subset(aux, v1 != v2)
      ind <- which.min(aux2[, "val"])
      mi <- aux2[ind, c("v1", "v2")]
      return(c(as.matrix(mi)))
   }

   # initial definitions (cluster 1)
   min1 <- fun.min(d)
   g <- list()
   ig <- 1
   g[[ig]] <- min1

   # (original) clustering criterion
   d. <- d
   diag(d.) <- NA
   theta <- max(apply(d., 2, min, na.rm = TRUE))
   criterion <- c()

   # clustering
   repeat {
      criterion[ig] <- theta
      newlab <- lab[-charmatch(unlist(g), lab)]
      n <- length(newlab)
      if (n < 1) break()
      m <- matrix(0, n + 1, n + 1)
      colnames(m) <- rownames(m) <- c("G", newlab)
      m[newlab, newlab] <- d[newlab, newlab]
      for(j in 1:n) {
          m["G", newlab[j]] <- m[newlab[j], "G"] <-
             mean(d[g[[ig]], newlab[j]])
      }
      comp <- newlab[which.min(m["G", newlab])]
      if (m["G", comp] <= theta) {
         g[[ig]] <- c(g[[ig]], comp)
      # -------------------------------------------
      # forming a new cluster
      } else {
         ig <- ig + 1
         if (n > 1) {
            # theta according to the algorithm
            theta <- ifelse(algorithm == "original", theta,
               max(apply(d.[newlab, newlab], 2, min, na.rm = TRUE)) )
            newcomp <- fun.min(d[newlab, newlab])
            if (d[newcomp[1], newcomp[2]] <= theta) {
               g[[ig]] <- newcomp
            } else {
               for(i in 1:n) g[[ig + i - 1]] <- newlab[i]
            }
         } else {
            g[[ig]] <- newlab
         }
      }
   }

   # output
   ng <- length(g)
   names(g) <- paste("cluster", 1:ng)
   class <- NULL
   for(k in 1:ng) {
      g[[k]] <- noquote(g[[k]])
      for(i in 1:ncol(d)) {
         if (any(lab[i] == g[[k]])) class[i] <- k
      }
   }
   nopc <- sapply(g, length)
   dc <- distClust(as.dist(d), nopc, unlist(g))
   out <- list(call = match.call(), 
      algorithm = algorithm,
      clusters = g, 
      class = class, 
      criterion = criterion, 
      distClust = dc, 
      d = as.dist(d))
   class(out) <- "tocher"
   return(out)
}

# ---------------------------------
# print method
print.tocher <- 
function (x, ...) 
{
   cat("\n          Tocher's Clustering \n\n")
   n <- attr(x$d, "Size")
   cat("Call: ")
      print(x$call)
   cat("\nCluster algorithm:", x$algorithm,
      "\nNumber of objects:", n,
      "\nNumber of clusters:", length(x$clusters), "\n")
   if (n < 25) {
      id <- maxmat(x$distClust)
      cat("Most contrasting clusters: ", id[1], " and ", id[2],
          ", with \n   average intercluster distance: ", 
             x$distClust[id[1], id[2]], "\n\n", sep = "")
      print(x$clusters)
   } else {
      cat("Most contrasting clusters (first 4 objects):\n")
      id <- maxmat(x$distClust)
      nopc <- sapply(x$clusters[id], length)
      rest <- max(nopc) - nopc
      dd <- data.frame(1:max(nopc))
      for(i in 1:length(nopc)) {
         dd[, i] <- c(x$clusters[[id[i]]], rep("", rest[id[i]]))
      }
      colnames(dd) <- paste(names(x$clusters[id]), 
         " (size ", nopc, ")", sep = "")
      print(head(dd, 4))
      cat("... with average intercluster distance:", 
         x$distClust[id[1], id[2]], "\n")
   }
   invisible(x)
}

# -------------------------------------------
# cophenetic method
coph.tocher <- function(x) {
   .Deprecated("cophenetic", package = "biotools")
   cophenetic.tocher(x=x)
}

cophenetic.tocher <-
function(x)
{
   dc <- x$distClust
   d <- x$d
   ord <- attr(d, "Labels")
   lab <- unlist(x$clusters)
   n <- length(lab)
   nc <- length(x$clusters)
   nopc <- sapply(x$clusters, length)
   cl <- rep(1:nc, nopc)
   coph <- matrix(0, n, n)
   dimnames(coph) <- list(lab, lab)
   for(j in 1:n) {
      for(i in 1:n) {
         if (i != j) {
            coph[lab[i], lab[j]] <- dc[cl[lab == lab[i]],
               cl[lab == lab[j]]]
         }
      }
   }
   return(as.dist(coph[ord, ord]))
}
