print.rc <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)
  cat("\nNormalized row scores:\n")
  print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
  cat("\nNormalized column scores:\n")
  print(format(ass$col[,,1], digits=digits, ...), quote=FALSE)

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.rc.symm <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)
  cat("\nNormalized scores:\n")
  print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.hmskew <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x[["assoc"]]

  if(length(ass) > 0) {
      cat("Intrinsic symmetric association coefficients:\n")
      print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)
      cat("\nNormalized symmetric association scores:\n")
      print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
      cat("\n")
  }

  ass <- x$assoc.hmskew

  cat("Intrinsic skew association coefficients:\n")
  print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)
  cat("\nNormalized skew association scores:\n")
  print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.yrcskew <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x[["assoc"]]

  if(length(ass) > 0) {
      cat("Intrinsic symmetric association coefficients:\n")
      print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)
      cat("\nNormalized symmetric association scores:\n")
      print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
      cat("\n")
  }

  ass <- x$assoc.yrcskew

  cat("\nIntrinsic skew association coefficients:\n")
  print(format(ass$phi[1,], digits=digits, ...), quote=FALSE)

  cat("\nNormalized skew association scores:\n")
  print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.rcL <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(format(ass$phi, digits=3), quote=FALSE)

  if(dim(ass$row)[3] == 1) {
      cat("\nNormalized row scores for all layers:\n")
      print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
  }
  else {
      cat("\nNormalized row scores:\n")
      print(format(ass$row, digits=digits, ...), quote=FALSE)
  }

  if(dim(ass$col)[3] == 1) {
      cat("\nNormalized column scores for all layers:\n")
      print(format(ass$col[,,1], digits=digits, ...), quote=FALSE)
  }
  else {
      cat("\nNormalized column scores:\n")
      print(format(ass$col, digits=digits, ...), quote=FALSE)
  }

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.rcL.symm <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x$assoc

  cat("Intrinsic association coefficients:\n")
  print(format(ass$phi, digits=digits, ...), quote=FALSE)

  if(dim(ass$row)[3] == 1) {
      cat("\nNormalized scores for all layers:\n")
      print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
  }
  else {
      cat("\nNormalized scores:\n")
      print(format(ass$row, digits=digits, ...), quote=FALSE)
  }

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.hmskewL <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x[["assoc"]]

  if(length(ass) > 0) {
      cat("Intrinsic symmetric association coefficients:\n")
      print(format(ass$phi, digits=digits, ...), quote=FALSE)

      if(dim(ass$row)[3] == 1) {
          cat("\nNormalized symmetric association scores for all layers:\n")
          print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
          cat("\n")
      }
      else {
          cat("\nNormalized symmetric association scores:\n")
          print(format(ass$row, digits=digits, ...), quote=FALSE)
          cat("\n")
      }
  }

  ass <- x$assoc.hmskew

  cat("Intrinsic skew association coefficients:\n")
  print(format(ass$phi, digits=3), quote=FALSE)

  if(dim(ass$row)[3] == 1) {
      cat("\nNormalized skew association scores for all layers:\n")
      print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)
  }
  else {
      cat("\nNormalized skew association scores:\n")
      print(format(ass$row, digits=digits, ...), quote=FALSE)
  }

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}

print.rcL.trans <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  ass <- x$assoc

  cat("Transition coefficients:\n")
  print(format(ass$transition, digits=3), quote=FALSE)

  cat("\nIntrinsic association coefficients:\n")
  print(format(ass$phi, digits=3), quote=FALSE)


  cat("\nNormalized row scores for first layer:\n")
  print(format(ass$row[,,1], digits=digits, ...), quote=FALSE)

  cat("\nNormalized row scores for last layer:\n")
  print(format(ass$row[,,dim(ass$row)[3]], digits=digits, ...), quote=FALSE)

  cat("\nVariation of normalized row scores\nbetween first and layer layer:\n")
  print(format(ass$row[,,dim(ass$row)[3]] - ass$row[,,1], digits=digits, ...), quote=FALSE)


  cat("\nNormalized column scores for first layers:\n")
  print(format(ass$col[,,1], digits=digits, ...), quote=FALSE)

  cat("\nNormalized column scores for last layer:\n")
  print(format(ass$col[,,dim(ass$col)[3]], digits=digits, ...), quote=FALSE)

  cat("\nVariation of normalized column scores\nbetween first and layer layer:\n")
  print(format(ass$row[,,dim(ass$col)[3]] - ass$col[,,1], digits=digits, ...), quote=FALSE)

  if(length(ass$diag) > 0) {
    cat("\nDiagonal coefficients:\n")
    print(format(ass$diag[1:nrow(ass$diag),], digits=digits, ...), quote=FALSE)
  }

  cat("\nNormalization weights:", ass$weighting)
  printModelStats(x, digits=digits)
}