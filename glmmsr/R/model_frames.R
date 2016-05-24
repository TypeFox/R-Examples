# convert formula to a model frame, using glFormula if there are random
# effects, and manually, using model.matrix, if not
parse_formula <- function(formula, data, family, weights, off) {
  formula_vars <- all.vars(formula)
  data_vars <- data[formula_vars]
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if(has_re(formula)) {
    if(is.null(off) && is.null(weights)) {
      modfr <- lme4::glFormula(formula, data = data_vars, family = family,
                               na.action = na.fail)
    } else {
      if(is.null(weights)) {
        modfr <- lme4::glFormula(formula, data = data_vars, family = family,
                                 na.action = na.fail, offset = off)
      } else {
        if(is.null(off)) {
          modfr <- lme4::glFormula(formula, data = data_vars, family = family,
                                   na.action = na.fail, weights = weights)
        } else {
          modfr <- lme4::glFormula(formula, data = data_vars, family = family,
                                   na.action = na.fail, weights = weights,
                                   offset = off)
        }
      }
    }
  } else {
    if(is.null(off) && is.null(weights)) {
      fr <- glm(formula, data = data_vars, family = family,
                method = "model.frame")
    } else {
      if(is.null(weights)) {
        data_vars_off <- data_vars
        data_vars_off$off <- off
        fr <- glm(formula, data = data_vars_off, family = family,
                  method = "model.frame", offset = off)
      } else {
        if(is.null(off)) {
          data_vars_weights <- data_vars
          data_vars_weights$weights <- weights
          fr <- glm(formula, data = data_vars_weights,
                    family = family, method = "model.frame",
                    weights = weights)
        } else {
          data_vars_off_weights <- data_vars
          data_vars_off_weights$off <- off
          data_vars_off_weights$weights <- weights
          fr <- glm(formula, data = data_vars_off_weights,
                    family = family, method = "model.frame",
                    offset = off, weights = weights)
        }

      }

    }
    X <- model.matrix(formula, data = fr)
    modfr <- list(fr = fr, X = X, family = family, formula = formula)
  }
  modfr
}

# convert formula to a model frame, using lFormula if there are random
# effects, and manually, using model.matrix, if not
parse_subformula <- function(formula, data) {
  formula_vars <- all.vars(formula)
  data_vars <- data[formula_vars]
  if(has_re(formula)) {
    modfr_tot <- lme4::lFormula(formula, data = data_vars,
                                control = lme4_control(), na.action = na.fail)
    modfr <- list(X = modfr_tot$X, reTrms = modfr_tot$reTrms)
  } else {

    fr <- lm(formula, data = data_vars, family = family,
             method = "model.frame")
    X <- model.matrix(formula, data = fr)
    modfr <- list(X = X)
  }
  modfr
}


has_reTrms <- function(modfr) {
  "reTrms" %in% names(modfr)
}

is_modfr <- function(x) {
  is.list(x) && ("X" %in% names(x))
}


`[fr` <- function(modfr, i) {
  out <- modfr
  i <- as.numeric(i)
  out$X <- modfr$X[i, , drop = FALSE]
  if(has_reTrms(modfr)) {
    reTrms <- modfr$reTrms
    out$reTrms$Zt <- reTrms$Zt[, i, drop = FALSE]
    out$reTrms$Ztlist <- lapply(reTrms$Ztlist, "[", j = i, drop = FALSE)
    out$reTrms$flist <- reTrms$flist[i, , drop = FALSE]
  }
  return(out)
}


# multiply a model frame by a constant
`*fr` <- function(modfr, a) {
  if(!is_modfr(modfr) || !is.numeric(a)) {
    if(is_modfr(a) && is.numeric(modfr)) {
      return(`*fr`(a, modfr))
    } else {
      stop("Invalid use of `*fr`", call. = FALSE)
    }
  }
  out <- modfr
  out$X <- modfr$X * a
  if(has_reTrms(modfr)) {
    out$reTrms$Zt <- modfr$reTrms$Zt * a
    out$reTrms$Ztlist <- lapply(modfr$reTrms$Ztlist, `*`, e2 = a)
  }
  return(out)
}

# divide a model frame by a constant
`/fr` <- function(modfr, a) {
  if(!is_modfr(modfr) || !is.numeric(a)) {
    stop("Invalid use of `/fr`", call. = FALSE)
  }
  `*fr`(modfr, 1/a)
}

# add two model frames
`+fr` <- function(modfr1, modfr2) {
  if(!is_modfr(modfr1) || !is_modfr(modfr2)) {
    stop("Invalid use of `+fr`: both arguments should be model frames",
         call. = FALSE)
  }
  if(has_reTrms(modfr1) != has_reTrms(modfr2)) {
    stop("Invalid use of `+fr`: Cannot add model frame with reTrms
         to model frame without",
         call. = FALSE)
  }
  out <- modfr1
  if(any(dim(modfr1$X) != dim(modfr2$X))) {
    stop("Invalid use of `+fr`: Dimension mismatch in model matrices",
         call. = FALSE)
  }
  out$X <- modfr1$X + modfr2$X
  if(has_reTrms(modfr1)) {
    if(any(dim(modfr1$reTrms$Zt) != dim(modfr2$reTrms$Zt))) {
      stop("Invalid use of `+fr`: Dimension mismatch in model matrices",
           call. = FALSE)
    }
    out$reTrms$Zt <- modfr1$reTrms$Zt + modfr2$reTrms$Zt
    out$reTrms$Ztlist <- mapply("+", modfr1$reTrms$Ztlist, modfr2$reTrms$Ztlist)
    # what about flist?
  }
  return(out)
}

# subtract one model frame from another
`-fr` <- function(modfr1, modfr2) {
  if(!is_modfr(modfr1) || !is_modfr(modfr2)) {
    stop("Invalid use of `-fr`: both arguments should be model frames",
         call. = FALSE)
  }
  minus_modfr2 <- `*fr`(modfr2, -1L)
  `+fr`(modfr1, minus_modfr2)
}

concatenate_Matrix <- function(M1, M2) {
  M1 <- as(M1, "dgTMatrix")
  M2 <- as(M2, "dgTMatrix")
  M_x <- c(M1@x, M2@x)
  dim1 <- dim(M1)
  dim2 <- dim(M2)
  M_i <- 1 + c(M1@i, dim1[1] + M2@i)
  M_j <- 1 + c(M1@j, dim1[2] + M2@j)
  Matrix::sparseMatrix(i = M_i, j = M_j, x = M_x, dims = dim1 + dim2)
}

concatenate_frames <- function(modfr1, modfr2) {
  out <- modfr1
  out$X <- cbind(modfr1$X, modfr2$X)
  if(has_reTrms(modfr1)) {
    if(has_reTrms(modfr2)) {
      reTrms1 <- modfr1$reTrms
      reTrms2 <- modfr2$reTrms
      out$reTrms$Zt <- Matrix::rBind(reTrms1$Zt, reTrms2$Zt)
      out$reTrms$theta <- c(reTrms1$theta, reTrms2$theta)
      out$reTrms$Lind <- c(reTrms1$Lind, reTrms2$Lind + length(reTrms1$theta))
      out$reTrms$Gp <- c(0L, reTrms1$Gp[2] + reTrms2$Gp[2])
      out$reTrms$lower <- c(reTrms1$lower, reTrms2$upper)
      out$reTrms$Lambdat <- concatenate_Matrix(reTrms1$Lambdat, reTrms2$Lambdat)
      out$reTrms$flist <- c(reTrms1$flist, reTrms2$flist)
      # need "assign" attribute to print out groups correctly
      attr(out$reTrms$flist, "assign") <- TRUE
      out$reTrms$cnms <- c(reTrms1$cnms, reTrms2$cnms)
      out$reTrms$Ztlist <- c(reTrms1$Ztlist, reTrms2$Ztlist)
    }
  } else {
    if(has_reTrms(modfr2)) {
      out$reTrms <- modfr2$reTrms
    }
  }
  out
}

attach_subframes <- function(modfr, subframes) {
  subframes_combined <- Reduce(concatenate_frames, subframes)
  concatenate_frames(modfr, subframes_combined)
}

