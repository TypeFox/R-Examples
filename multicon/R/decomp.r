decomp=function(x, y=NULL, na.rm=TRUE, use="pair") {
	x=as.matrix(x)##
  GMx <- mean(x, na.rm=na.rm)
  rowEffx <- rowMeans(x) - GMx
  colEffx <- colMeans(x) - GMx
  partx <- x - GMx - matrix(rowEffx, ncol=ncol(x), nrow=nrow(x)) - matrix(colEffx, ncol=ncol(x), nrow=nrow(x), byrow=T)
  rowMS <- (ncol(x) * sum(rowEffx^2)) / (nrow(x) - 1)
  colMS <- (nrow(x) * sum(colEffx^2)) / (ncol(x) - 1)
  partMS <- sum(partx^2) / ((nrow(x)-1)*(ncol(x)-1))
  rowVar <- (rowMS - partMS) / ncol(x)
  colVar <- (colMS - partMS) / nrow(x)
  partVar <- partMS
  TotVar <- ifelse(rowVar < 0, 0, rowVar) + ifelse(colVar < 0, 0, colVar) + ifelse(partVar < 0, 0, partVar)
  rowVarPerc <- ifelse(rowVar < 0, 0, rowVar / TotVar)
  colVarPerc <- ifelse(colVar < 0, 0, colVar / TotVar)
  partVarPerc <- ifelse(partVar < 0, 0, partVar / TotVar)
  rowDFnum <- nrow(x) - 1
  rowDFden <- (ncol(x)-1) * (nrow(x)-1)
  colDFnum <- ncol(x) - 1
  colDFden <- (ncol(x)-1) * (nrow(x)-1)
  rowF <- rowMS / partMS
  colF <- colMS / partMS
  rowp <- 1 - pf(rowF, rowDFnum, rowDFden)
  colp <- 1 - pf(colF, colDFnum, colDFden)
  VarComp <- data.frame("Var"=c(rowVar, colVar, partVar), "VarPerc"=c(rowVarPerc, colVarPerc, partVarPerc), "NumDF"=c(rowDFnum, colDFnum, NA), "DenDF"=c(rowDFden, colDFden, NA), "F"=c(rowF, colF, NA), "p"=c(rowp, colp, NA), 
                row.names=c("Row Variance", "Column Variance", "Interaction Variance"))
                  
  if(!is.null(y)) {
  	y=as.matrix(y)##
    GMy <- mean(y, na.rm=na.rm)
    rowEffy <- rowMeans(y, na.rm=na.rm) - GMy
    colEffy <- colMeans(y, na.rm=na.rm) - GMy
    party <- y - GMy - matrix(rowEffy, ncol=ncol(y), nrow=nrow(y)) - matrix(colEffy, ncol=ncol(y), nrow=nrow(y), byrow=T)
    rowCorUniq <- diag(cor(partx, party))
    colCorUniq <- diag(cor(t(partx), t(party)))
    elevation <- GMx - GMy
    stereoAcc <- cor(rowEffx, rowEffy, use=use)
    diffElevation <- cor(colEffx, colEffy, use=use)
    diffAcc <- cor(as.vector(partx), as.vector(party), use=use)
    stats <- rbind(elevation, stereoAcc, diffElevation, diffAcc)
    colnames(stats) <- "Results"
    rownames(stats) <- c("Elevation", "Stereotype Accuracy", "Differential Elevation", "Differential Accuracy")
    out <- list(GMx, GMy, rowEffx, rowEffy, colEffx, colEffy, partx, party, rowCorUniq, colCorUniq, VarComp, stats)
    names(out) <- c("GrandMeanX", "GrandMeanY", "RowEffectX", "RowEffectY", "ColEffectX", "ColEffectY", "DecompositionX", "DecompositionY", "RowUniqueCor", "ColUniqueCor", "VarComp", "Stats")
  }
  if(is.null(y)) {
    out <- list(GMx, rowEffx, colEffx, partx, VarComp)
    names(out) <- c("GrandMeanX", "RowEffectX", "ColEffectX", "DecompositionX", "VarComp")    
  }
  out
}