################################################################################
# print.ca(): Printing ca objects (ca package 0.70)
################################################################################
print.ca <- function(x, ...){
  obj <- x
  nd0 <- length(obj$sv)
  nd  <- obj$nd
  if (is.na(nd)){
    nd <- 2
    } else {
    if (nd > length(obj$sv)) nd <- length(obj$sv)
    }
 # Eigenvalues:
  Dimension <- 1:nd0
  Value <- round(obj$sv^2, 6)
  Percentage <- paste(as.character(round(100 * Value / sum(Value), 2)), "%", sep = "")
  tmp <- rbind(Value = as.character(Value), Percentage = as.character(Percentage))
  dimnames(tmp)[[2]] <- Dimension  
  Eigenvalues <- tmp
 # Row Profiles:
  tmp      <- rbind(obj$rowmass, obj$rowdist, obj$rowinertia, t(obj$rowcoord[,1:nd]))
  tmpnames <- obj$rownames
  if (!is.na(obj$rowsup[1])){
    tmpnames[obj$rowsup] <- paste(tmpnames[obj$rowsup],"(*)")
    }
  dimnames(tmp)[[2]] <- tmpnames
  dn <- paste("Dim.", 1:nd)
  dimnames(tmp)[[1]] <- c("Mass", "ChiDist", "Inertia", dn)
  Row.profiles <- tmp
 # Column Profiles:
  tmp <- rbind(obj$colmass, obj$coldist, obj$colinertia, t(obj$colcoord[,1:nd]))
  tmpnames <- obj$colnames
  if (!is.na(obj$colsup[1])){
    tmpnames[obj$colsup] <- paste(tmpnames[obj$colsup],"(*)")
    }
  dimnames(tmp)[[2]] <- tmpnames
  dn <- paste("Dim.", 1:nd)
  dimnames(tmp)[[1]] <- c("Mass", "ChiDist", "Inertia", dn)
  Column.profiles <- tmp
  cat("\n Principal inertias (eigenvalues):\n")
  print.table(Eigenvalues, width = 4)
  cat("\n\n Rows:\n")
  print(round(Row.profiles, 6))
  cat("\n\n Columns:\n")
  print(round(Column.profiles, 6))
  }
################################################################################
