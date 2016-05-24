
sup.scores.rc <- function(model, tab, ass, rowsup, colsup,
                          symmetry=c("asymmetric", "symmetric", "skew-symmetric"), str="Mult") {
  symmetry <- match.arg(symmetry)
  stopifnot(!is.null(rowsup) || !is.null(colsup))
  stopifnot(symmetry =="asymmetric" || nrow(rowsup) == ncol(colsup))

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  phi <- ass$phi
  row <- ass$row
  col <- ass$col
  weighting <- ass$weighting

  row.sup <- NULL
  col.sup <- NULL
  row.weights.sup <- NULL
  col.weights.sup <- NULL

  args <- list()
  args$formula <- model$formula
  args$family <- model$family
  args$tolerance <- model$tolerance
  args$iterMax <- model$iterMax
  args$verbose <- FALSE


  # Diagonal does not make sense in this model since table is not square,
  # or has NA diagonal if both rowsup and colsup have been supplied
  if(grepl("Diag", deparse(args$formula)))
      args$formula <- as.formula(paste(deparse(args$formula),
                                       sprintf("- Diag(%s, %s)", vars[1], vars[2])))

  if(is.null(colsup)) {
      args$data <- as.table(rowsup)
  }
  else if(is.null(rowsup)) {
      args$data <- as.table(colsup)
  }
  else {
      # Block matrix that is symmetric if rowsup and colsup are transposes of each other
      args$data <- as.table(cbind(rbind(matrix(NA, nrow(colsup), ncol(rowsup)), rowsup),
                                  rbind(colsup, matrix(NA, nrow(rowsup), ncol(colsup)))))
      dimnames(args$data) <- list(c(rownames(model$data), rownames(rowsup)),
                                  c(colnames(model$data), colnames(colsup)))
  }

  names(dimnames(args$data)) <- names(dimnames(model$data))

  if(symmetry == "skew-symmetric") {
      hmterm <- sprintf("+ HMSkew(%s, %s)", vars[1], vars[2])

      args2 <- args
      args2$formula <- as.formula(sub(paste("\\Q", hmterm, "\\E", sep=""),
                                      "", deparse(args$formula)))

      base <- do.call("gnm", args2)
      args$start <- c(rep(NA, length(parameters(base))), residEVD(base, 1, skew=TRUE))

      # -1 is here to ensure all symmetric association coefficients are in the list
      args2$formula <- as.formula(paste("Freq ~ -1", hmterm))
      args2$method <- "coefNames"
      hmnames <- do.call("gnm", args2)

      cnames <- c(names(parameters(base)), hmnames)
  }
  else {
      args2 <- args
      args2$method <- "coefNames"

      cnames <- do.call("gnm", args2)
  }

  args$constrain <- which(cnames %in% names(parameters(model)) & grepl(str, cnames))

  if(symmetry != "asymmetric")
      args$constrainTo <- sweep(cbind(row[,, 1]), 2, sqrt(phi[1,]), "*")
  else if(length(rowsup) > 0 && length(colsup) > 0)
      args$constrainTo <- sweep(rbind(cbind(row[,, 1]), cbind(col[,, 1])), 2, sqrt(phi[1,]), "*")
  else if(length(rowsup) > 0)
      args$constrainTo <- sweep(cbind(col[,, 1]), 2, sqrt(phi[1,]), "*")
  else
      args$constrainTo <- sweep(cbind(row[,, 1]), 2, sqrt(phi[1,]), "*")

  msup <- do.call("gnm", args)

  sup <- parameters(msup)[setdiff(pickCoef(msup, str), args$constrain)]

  if(symmetry != "asymmetric") {
      dim(sup) <- c(sum(nrow(rowsup)), ncol(phi))

      if(weighting == "none")
          rpsup <- rep(1, nrow(rowsup))
      else if(weighting == "uniform")
          rpsup <- rep(1/nrow(rowsup), nrow(rowsup))
      else
          rpsup <- prop.table(apply(rowsup, 1, sum, na.rm=TRUE) + apply(colsup, 2, sum, na.rm=TRUE))

      rsup <- sup[seq(nrow(rowsup)), , drop=FALSE]

      # Linear constrain of zero sum on scores
      rsup <- sweep(rsup, 2, colSums(sweep(rsup, 1, rpsup/sum(rpsup), "*")), "-")

      # Supplementary scores do not need a zero sum of squares constrain,
      # but they have to be scaled the same way as standardized scores
      rsup <- sweep(rsup, 2, sqrt(phi), "/")

      row2 <- rbind(cbind(row[,,1]), rsup)
      dim(row2)[3] <- 1
      names(dimnames(row2)) <- names(dimnames(row))
      rownames(row2) <- c(rownames(row), rownames(rowsup))
      colnames(row2) <- colnames(row)
      row.sup <- col.sup <- seq(nrow(row) + 1, nrow(row2))
      row <- col <- row2

      if(length(dim(tab)) == 3) {
          row.weights.sup <- apply(rowsup, c(1, 3), sum, na.rm=TRUE)
          col.weights.sup <- apply(colsup, c(2, 3), sum, na.rm=TRUE)
      }
      else {
          row.weights.sup <- as.matrix(apply(rowsup, 1, sum, na.rm=TRUE))
          col.weights.sup <- as.matrix(apply(colsup, 2, sum, na.rm=TRUE))
      }
  }
  else {
      dim(sup) <- c(sum(nrow(rowsup), ncol(colsup)), ncol(phi))

      if(length(rowsup) > 0) {
          if(weighting == "none")
              rpsup <- rep(1, nrow(rowsup))
          else if(weighting == "uniform")
              rpsup <- rep(1/nrow(rowsup), nrow(rowsup))
          else
              rpsup <- prop.table(apply(rowsup, 1, sum, na.rm=TRUE))

          rsup <- sup[seq(nrow(rowsup)), , drop=FALSE]

          # Linear constrain of zero sum on scores
          rsup <- sweep(rsup, 2, colSums(sweep(rsup, 1, rpsup/sum(rpsup), "*")), "-")

          # Supplementary scores do not need a zero sum of squares constrain,
          # but they have to be scaled the same way as standardized scores
          rsup <- sweep(rsup, 2, sqrt(phi), "/")

          row2 <- rbind(cbind(row[,,1]), rsup)
          dim(row2)[3] <- 1
          names(dimnames(row2)) <- names(dimnames(row))
          rownames(row2) <- c(rownames(row), rownames(rowsup))
          colnames(row2) <- colnames(row)
          row.sup <- seq(nrow(row) + 1, nrow(row2))
          row <- row2

          if(length(dim(tab)) == 3)
              row.weights.sup <- apply(rowsup, c(1, 3), sum, na.rm=TRUE)
          else
              row.weights.sup <- as.matrix(apply(rowsup, 1, sum, na.rm=TRUE))
      }

      if(length(colsup) > 0) {
          if(weighting == "none")
              cpsup <- rep(1, ncol(colsup))
          else if(weighting == "uniform")
              cpsup <- rep(1/ncol(colsup), ncol(colsup))
          else
              cpsup <- prop.table(apply(colsup, 2, sum, na.rm=TRUE))

          csup <- sup[seq(NROW(rowsup) + 1, nrow(sup)), , drop=FALSE]

          # Linear constrain of zero sum on scores
          csup <- sweep(csup, 2, colSums(sweep(csup, 1, cpsup/sum(cpsup), "*")), "-")

          # Supplementary scores do not need a zero sum of squares constrain,
          # but they have to be scaled the same way as standardized scores
          csup <- sweep(csup, 2, sqrt(phi), "/")

          col2 <- rbind(cbind(col[,,1]), csup)
          dim(col2)[3] <- 1
          names(dimnames(col2)) <- names(dimnames(col))
          rownames(col2) <- c(rownames(col), colnames(colsup))
          colnames(col2) <- colnames(col)
          col.sup <- seq(nrow(col) + 1, nrow(col2))
          col <- col2

          if(length(dim(tab)) == 3)
              col.weights.sup <- apply(colsup, c(2, 3), sum, na.rm=TRUE)
          else
              col.weights.sup <- as.matrix(apply(colsup, 1, sum, na.rm=TRUE))
      }
  }

  list(row=row, col=col,
       row.weights=rbind(ass$row.weights, row.weights.sup),
       col.weights=rbind(ass$col.weights, col.weights.sup),
       row.sup=row.sup,
       col.sup=col.sup)
}
