# Sequential Bifurcations (Bettonvil, Kleijnen 1997)
#
# Gilles Pujol 2007

# TO DO: store computed effects in the object, hard-coded loop, return
# warning when hypotheses violated.


# sbgroups : returns the table the effects per groups

sbgroups <- function(x) {

  # sort the data on i

  y <- x$y
  y[rank(x$i)] <- y
  if (! is.null(x$ym)) {
    ym <- x$ym
    ym[rank(x$i)] <- ym
  }
  i <- sort(x$i)
  
  # groups

  n <- length(i)
  lower <- i[1 : (n - 1)] + 1
  upper <- i[2 : n]

  # effects
  
  if (is.null(x$ym)) {
    effect = diff(y)
  } else {
    effect = (diff(y) - diff(ym)) / 2
  }

  return(data.frame(lower = lower, upper = upper, effect = effect))
}


sb <- function(p, sign = rep("+", p), interaction = FALSE) {
  x <- list(p = p, sign = sign, interaction = interaction, i = NULL)
  # (remark: i = NULL is necessary, else x$i would match x$interaction...)
  class(x) <- "sb"
  return(x)
}


ask.sb <- function(x, i = NULL, ...) {
  out <- list()

  if (! is.null(i)) {

    # 1. if i is given, return x(i), where i can be negative, or a vector
    
    for (j in i) {
      if (j >= 0) {
        tmp <- c(rep(1, j), rep(-1, x$p - j))
      } else {
        tmp <- c(rep(-1, -j), rep(1, x$p + j))
      }
      tmp[x$sign == "-"] <- - tmp[x$sign == "-"]
      out[paste(j)] <- list(tmp)
    }
  } else {

    # 2. if i is not given, the function proposes the 'best' i and
    #    the correspondind x(i)
    
    if (length(x$i) == 0) {
      out["0"] <- Recall(x, 0)
      out[paste(x$p)] <- Recall(x, x$p)
    } else {
      
      # table of groups that can be subdivied (i.e non reduced to one element)
      
      T <- sbgroups(x)
      T <- T[T[,1] != T[,2],]

      # groups to subdivide

      g <- T[which.max(T[,3]), c(1,2)]

      # bifurcation

      j <- floor((g[1] + g[2])/ 2) # dichotomy
      out[paste(j)] <- Recall(x, j)
      if (x$interaction) {     # if interactions, mirror observation
        out[paste(-j)] <- Recall(x, -j)
      }
    }
  }

  return(out)
}


tell.sb <- function(x, y, ...) {
  id <- deparse(substitute(x))

  # (y is a list where the labels are "i" or "-i")
  
  # if interactions, the observations x(0) dans x(p)
  # are duplicated to provide x(-p) and x(-0)

  if (x$interaction) {
    if ("0" %in% names(y)) {
      y[[paste(-x$p)]] <- y[["0"]]
    }
    if (paste(x$p) %in% names(y)) {
      y[["-0"]] <- y[[paste(x$p)]]
    }
  }
  
  for (i in names(y)) {
    j <- as.numeric(i) # warning: 0 and -0 will be both converted into "0"

    # if there is not a line corresponding to i, it is created
    if (! (abs(j) %in% x$i)) {
      x$i <- c(x$i, abs(j))
      x$y <- c(x$y, NA)
      if (x$interaction) {
        x$ym <- c(x$ym, NA)
      }
    }

    # number of the line corresponding to i
    ligne <- (1 : length(x$i))[x$i == abs(j)]

    # y is set  this line, at the right place...
    if (paste(abs(j)) == i) {
      # ... if "i": 2nd column
      x$y[ligne] <- y[[i]]
    } else {
      # ... if "-i": 3rd column
      x$ym[ligne] <- y[[i]]
    }
  }

  assign(id, x, envir = parent.frame())
}


print.sb <- function(x, ...) {
  if (! is.null(x$i)) {
##     cat("\nObservations:\n")
##     if (! x$interaction) {
##       print(data.frame(i = x$i, observation = x$y))
##     } else {
##       print(data.frame(i = x$i, observation = x$y, mirror = x$ym))
##     }
    
    cat("\nGroups:\n")
    T <- sbgroups(x)
    groupnames <- paste(T[, 1], T[, 2], sep = "-")
    groupnames[T[, 1] == T[, 2]] <- paste(T[T[, 1] == T[, 2], 1])
    print(data.frame(group = groupnames, effect = T[, 3]))
  } else {
    cat("(empty)\n")
  }
}


plot.sb <- function(x, ...) {
  if (! is.null(x$i)) {
    T <- sbgroups(x)
    n <- nrow(T)

    xmax <- max(T[, c(1, 2)])
    ymax <- max(T[, 3])
    plot(0, type = "n", xlim = c(1, xmax), ylim = c(0, ymax), axes = FALSE)
    axis(1, at = 1 : xmax)
    axis(2)
    box()
    lines(c(1, xmax), c(0, 0))
    
    for (i in 1 : n){
      if (T[i, 1] == T[i, 2]){
        lines(c(T[i, 1], T[i, 1]), c(0, T[i, 3]))
      } else {
        lines(c(T[i, 1], T[i, 1], T[i, 2], T[i, 2]),
              c(0, T[i, 3], T[i, 3], 0))
      }
    }
    
  }
}
