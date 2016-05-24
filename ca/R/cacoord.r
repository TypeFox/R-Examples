################################################################################
# cacoord(): Extracting CA/MCA coordinates (ca 0.64)
# Arguments
#   - obj : A 'ca' or 'mjca' object.
#   - type: The type of coordinates ("standard" or "principal"); the remaining 
#           options ("symmetric", ..., "colgreen") return the corresponding 
#           row/column configuration for the map scaling options in plot.ca().
#   - dim : The dimensions to return; if NA, all available dimensions are 
#           returned.
#   - rows: If TRUE, the row coordinates are returned (see below for details).
#   - cols: If TRUE, the column coordinates are returned (see below for details).
#   -  ...: Further arguments (ignored).
# Value
#   - A list with the entries 'rows' and 'columns' containing the corresponding
#     row and column coordinates (for (rows=NA&cols=NA)|(rows=TRUE&cols=TRUE)).
#     For 'rows=TRUE' (and cols=NA or cols=FALSE) a matrix with the row 
#     coordinates is returned; and for 'cols=TRUE' (and cols=NA or cols=FALSE) 
#     a matrix with the column coordinates is returned.
################################################################################
cacoord <- function(obj, 
             type = c("standard", "principal", "symmetric", "rowprincipal", "colprincipal", 
                      "symbiplot", "rowgab", "colgab", "rowgreen", "colgreen"), 
             dim  = NA,
             rows = NA,
             cols = NA,
             ...){

  if (!inherits(obj, c("ca", "mjca"))){
    stop("'obj' must be a 'ca' or 'mjca' object")
    }
  map <- match.arg(type)
  if (is.na(rows) & is.na(cols)){
    rows <- TRUE
    cols <- TRUE
    } else{
    if (is.na(rows) | !rows){
      rows         <- FALSE
      cols         <- TRUE
      obj$rowcoord <- matrix(rep(0, ncol(obj$colcoord)), nrow = 1)
      obj$rowmass  <- 1
      }
    if (is.na(cols) | !cols){
      cols         <- FALSE
      rows         <- TRUE
      obj$colcoord <- matrix(rep(0, ncol(obj$rowcoord)), nrow = 1)
      obj$colmass  <- 1
      }
    }
 # Check row-/columnnames:
   if (is.null(rownames(obj$rowcoord))){
     x.rnames <- 1:nrow(obj$rowcoord)
     rownames(obj$rowcoord) <- x.rnames
     } else {
     x.rnames <- rownames(obj$rowcoord)
     }
   if (is.null(colnames(obj$rowcoord))){
     x.cnames <- paste("Dim", 1:ncol(obj$rowcoord), sep = "")
     colnames(obj$rowcoord) <- x.cnames
     } else {
     x.cnames <- colnames(obj$rowcoord)
     }

   if (is.null(rownames(obj$colcoord))){
     y.rnames <- 1:nrow(obj$colcoord)
     rownames(obj$colcoord) <- y.rnames
     } else {
     y.rnames <- rownames(obj$colcoord)
     }
   if (is.null(colnames(obj$colcoord))){
     y.cnames <- paste("Dim", 1:ncol(obj$colcoord), sep = "")
     colnames(obj$colcoord) <- y.cnames
     } else {
     y.cnames <- colnames(obj$colcoord)
     }
 # Extract dimensions
  if (is.na(dim)[1]){
    sv  <- obj$sv
    rsc <- obj$rowcoord
    csc <- obj$colcoord
    } else {
    sv  <- obj$sv[dim]
    rsc <- matrix(obj$rowcoord[,dim], ncol = length(dim))
    csc <- matrix(obj$colcoord[,dim], ncol = length(dim))
    rownames(rsc) <- x.rnames
    colnames(rsc) <- x.cnames[dim]
    rownames(csc) <- y.rnames
    colnames(csc) <- y.cnames[dim]
    }
 # Coordinates ([r,c]sc: standard; [r,c]pc: principal; sym[r,c]pc: biplot; [r,c]gab: gabriel; [r,c]green: greenacre):
  if (map == "standard"){
    x <- rsc
    y <- csc
    } else {
    I   <- nrow(rsc)
    J   <- nrow(csc)
    K   <- ncol(rsc)
    rpc <- rsc %*% diag(sv)
    cpc <- csc %*% diag(sv)
    if (map == "principal"){
      x <- rpc
      y <- cpc
      } else {
      symrpc <- rsc %*% diag(sqrt(sv))
      symcpc <- csc %*% diag(sqrt(sv))
      rgab   <- rsc * matrix(obj$rowmass, ncol = ncol(rsc), nrow = nrow(rsc))
      cgab   <- csc * matrix(obj$colmass, ncol = ncol(csc), nrow = nrow(csc))
      rgreen <- rsc * matrix(sqrt(obj$rowmass), ncol = ncol(rsc), nrow = nrow(rsc))
      cgreen <- csc * matrix(sqrt(obj$colmass), ncol = ncol(csc), nrow = nrow(csc))
     # Maptype LUT
      mt    <- c("symmetric", "rowprincipal", "colprincipal", "symbiplot", 
                 "rowgab", "colgab", "rowgreen", "colgreen")
      mti   <- 1:length(mt)
      mtlut <- list(symmetric    = list(x = rpc,    y = cpc),
                    rowprincipal = list(x = rpc,    y = csc),
                    colprincipal = list(x = rsc,    y = cpc),
                    symbiplot    = list(x = symrpc, y = symcpc), 
                    rowgab       = list(x = rpc,    y = cgab),
                    colgab       = list(x = rgab,   y = cpc), 
                    rowgreen     = list(x = rpc,    y = cgreen), 
                    rowgreen     = list(x = rgreen, y = cpc)
                    )
      x <- mtlut[[mti[mt == map]]][[1]]
      y <- mtlut[[mti[mt == map]]][[2]]
      } # End !"principal"
    } # End !"standard"
 # Fix row-/columnnames
  rownames(x) <- rownames(rsc)
  colnames(x) <- colnames(rsc)
  rownames(y) <- rownames(csc)
  colnames(y) <- colnames(csc)
 # Return rows and/or columns:
  if (rows & cols){ 
    out <- list(rows = x, columns = y)
    } else {
    if (rows){
      out <- x
      } else {
      out <- y
      }
    }
  return(out)
  }
################################################################################

