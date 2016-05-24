## 'value' is a unique named value (compatible with unlist())

addSimplexEq <- function(simplex, ## a single simplex 
                         value) { ## finds where a scalar equality intersects simplex edges 
  if (is.null(simplex)) {
    resu <- NULL
  } else {
    value <- unlist(value)
    var <- names(value)
    if (length(var)!=1) stop.redef("(!) From addVeqs: values should be a unique, named, scalar value.")
    left <- which(simplex[, var]<value)
    right <- which(simplex[, var]>value)
    if (length(left)==0) {
      left <- which(simplex[, var]==value)
      if (length(left)==0) {
        resu <- NULL ## unfeasible constraint
      } else {
        resu <- simplex[left, , drop=FALSE]
      }
    } else if (length(right)==0) {
      right <- which(simplex[, var]==value)
      if (length(right)==0) {
        resu <- NULL ## unfeasible constraint
      } else {
        resu <- simplex[right, , drop=FALSE]
      }
    } else { ## non-empty left and right
      ## greedy: find where left-right segments meet the constraint
      ## FR->FR in a more evolved version I should select facets that involve both left and right, then the points for these facets...
      left <- simplex[left, , drop=FALSE]
      right <- simplex[right, , drop=FALSE]
      blob  <-  proxy::dist(right[, var], left[, var], function(r, l, value) {(value-l)/(r-l)} , value=value)  ## nrow(right) X nrow(left)
      db  <-  dim(blob)
      blobij  <-  matrix(nrow=prod(db), ncol=2)
      blobij[, 1]  <-  rep.int(seq_len(db[1]), db[2])
      blobij[, 2]  <-  rep.int(seq_len(db[2]), rep.int(db[1], db[2])) ## cf expand.grid...
      resu  <-  (1-blob)[blobij]*left[blobij[, 2], ]+blob[blobij]*right[blobij[, 1], ]
      if( ! is.null(resu) && ! is.matrix(resu)) {
        ## resu is a data.frame and matrix(<df>) gives an object of class matrix but str() shows a list
        ## as.matrix() has no nrow and ncol args, (!) and  as.matrix(c(1,2), nrow = 1,ncol=2) gives the opposite of the expected result...
        # ? as.matrix suggests:
        dim(resu) <- c(1,length(resu))
        colnames(resu) <- colnames(simplex)
      }
    }
  }
  return(resu)
} ## end addVeq
