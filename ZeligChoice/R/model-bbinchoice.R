#' Bivariate Binary Choice object for inheritance across models in ZeligChoice
#'
#' @import methods
#' @export Zelig-bbinchoice
#' @exportClass Zelig-bbinchoice



zbbinchoice <- setRefClass("Zelig-bbinchoice",
                          contains = "Zelig",
                          field = list(family = "ANY",
                                       linkinv = "function"
                          ))

zbbinchoice$methods(
  initialize = function() {
    callSuper()
    .self$fn <- quote(VGAM::vglm)
    .self$authors <- "Kosuke Imai, Gary King, Olivia Lau"
    .self$packageauthors <- "Thomas W. Yee"
    .self$year <- 2007
    .self$category <- "dichotomous"
  }
)

zbbinchoice$methods(
  zelig = function(formula, data, ..., weights = NULL, by = NULL) {
    .self$zelig.call <- match.call(expand.dots = TRUE)
    .self$model.call <- match.call(expand.dots = TRUE)
    .self$model.call$family <- .self$family
    # Zelig 4 compatibility layer
#     if (is.list(formula)) {
#       formula[[1]] <- formula(f1)
#       formula[[2]] <- formula(f2)
#       f1 <- f[[1]]
#       f2 <- f[[2]]
#       formula <- update(f1, cbind(f1, f2) ~ . )
#     }
    callSuper(formula = formula, data = data, ..., weights = NULL, by = by)
  }
)

zbbinchoice$methods(
  param = function(z.out) {
    return(mvrnorm(.self$num, coef(z.out), vcov(z.out)))
  }
)

zbbinchoice$methods(
  # From Zelig 4
  qi = function(simparam, mm) {
    .pp <- function(object, constr, all.coef, x) {
      xm <- list()
      xm <- rep(list(NULL), 3)
      sim.eta <- NULL
      for (i in 1:length(constr))
        for (j in 1:3)
          if (sum(constr[[i]][j,]) == 1)
            xm[[j]] <- c(xm[[j]], x[,names(constr)[i]])
      sim.eta <- cbind(
        all.coef[[1]] %*% as.matrix( xm[[1]] ),
        all.coef[[2]] %*% as.matrix( xm[[2]] ),
        all.coef[[3]] %*% as.matrix( xm[[3]] )
      )
      # compute inverse (theta)
      ev <- .self$linkinv(sim.eta)
      # assign correct column names
      colnames(ev) <- c("Pr(Y1=0, Y2=0)",
                        "Pr(Y1=0, Y2=1)",
                        "Pr(Y1=1, Y2=0)",
                        "Pr(Y1=1, Y2=1)"
      )
      return(ev)
    }
    
    .pr <- function(ev) {
      mpr <- cbind(ev[, 3] + ev[, 4], ev[, 2] + ev[, 4])
      index <- matrix(NA, ncol=2, nrow=nrow(mpr))
      index[, 1] <- rbinom(n=nrow(ev), size=1, prob=mpr[, 1])
      index[, 2] <- rbinom(n=nrow(ev), size=1, prob=mpr[, 2])
      pr <- matrix(NA, nrow=nrow(ev), ncol=4)
      pr[, 1] <- as.integer(index[, 1] == 0 & index[, 2] == 0)
      pr[, 2] <- as.integer(index[, 1] == 0 & index[, 2] == 1)
      pr[, 3] <- as.integer(index[, 1] == 1 & index[, 2] == 0)
      pr[, 4] <- as.integer(index[, 1] == 1 & index[, 2] == 1)
      colnames(pr) <- c("(Y1=0, Y2=0)",
                        "(Y1=0, Y2=1)",
                        "(Y1=1, Y2=0)",
                        "(Y1=1, Y2=1)")
      return(pr)
    }
    .make.match.table <- function(index, cols=NULL) {
      pr <- matrix(0, nrow=nrow(index), ncol=4)
      # assigns values by the rule:
      #   pr[j,1] = 1 iff index[j,1] == 0 && index[j,2] == 0
      #   pr[j,2] = 1 iff index[j,1] == 0 && index[j,2] == 1
      #   pr[j,3] = 1 iff index[j,1] == 1 && index[j,2] == 0
      #   pr[j,4] = 1 iff index[j,1] == 1 && index[j,2] == 1
      # NOTE: only one column can be true at a time, so as a result
      #       we can do a much more elegant one liner, that I'll code
      #       later.  In this current form, I don't think this actually
      #       explains what is going on.
      pr[, 1] <- as.integer(index[, 1] == 0 & index[, 2] == 0)
      pr[, 2] <- as.integer(index[, 1] == 0 & index[, 2] == 1)
      pr[, 3] <- as.integer(index[, 1] == 1 & index[, 2] == 0)
      pr[, 4] <- as.integer(index[, 1] == 1 & index[, 2] == 1)
      # assign column names
      colnames(pr) <- if (is.character(cols) && length(cols)==4)
        cols
      else
        c("(Y1=0, Y2=0)",
          "(Y1=0, Y2=1)",
          "(Y1=1, Y2=0)",
          "(Y1=1, Y2=1)")
      return(pr)
    }
    all.coef <- NULL
    coefs <- simparam
    cm <- constraints(.self$zelig.out$z.out[[1]])
    v <- vector("list", 3)
    for (i in 1:length(cm)) {
      if (ncol(cm[[i]]) == 1){
        for (j in 1:3)
          if (sum(cm[[i]][j, ]) == 1)
            v[[j]] <- c(v[[j]], names(cm)[i])
      }
      else {
        for (j in 1:3)
          if (sum(cm[[i]][j,]) == 1)
            v[[j]] <- c(v[[j]], paste(names(cm)[i], ":", j, sep=""))
      }
    }
    for(i in 1:3)
      all.coef[[i]] <- coefs[ , unlist(v[i]) ]
    col.names <- c("Pr(Y1=0, Y2=0)",
                   "Pr(Y1=0, Y2=1)",
                   "Pr(Y1=1, Y2=0)",
                   "Pr(Y1=1, Y2=1)"
    )
    ev <- .pp(.self$zelig.out$z.out[[1]], cm, all.coef, as.matrix(mm))
    pv <- .pr(ev)
    levels(pv) <- c(0, 1)
#     return(list("Predicted Probabilities: Pr(Y1=k|X)" = ev,
#                 "Predicted Values: Y=k|X" = pv))
    return(list(ev = ev, pv = pv))
  }
)

# zbinchoice$methods(
#   show = function() {
#     lapply(.self$zelig.out, function(x) print(VGAM::summary(x)))
#   }
# )

