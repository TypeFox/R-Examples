"write.crd" <-
function(pdb = NULL, xyz = pdb$xyz, resno = NULL, resid = NULL,
                      eleno = NULL, elety = NULL, segid = NULL, resno2 = NULL,
                      b = NULL, verbose = FALSE, file = "R.crd") {
  
    if (is.null(xyz) || !is.numeric(xyz))
        stop("write.crd: please provide a 'pdb' object or numeric 'xyz' coordinates")
    if (any(is.na(xyz)))
        stop("write.crd: 'xyz' coordinates must have no NA's.")
    if (is.matrix(xyz) && nrow(xyz) == 1) xyz = as.vector(xyz)
    if (is.vector(xyz)) {
        natom <- length(xyz)/3
    } else {
      stop("write.crd: 'xyz' or 'pdb$xyz' must contain only one structure")
    }

    if (!is.null(pdb)) {
        if (is.null(resno))
          resno = pdb$atom[, "resno"]
        if (is.null(resid))
          resid = pdb$atom[, "resid"]
        if (is.null(eleno))
          eleno = pdb$atom[, "eleno"]
        if (is.null(elety))
          elety = pdb$atom[, "elety"]
        if (is.null(segid)) 
          segid = pdb$atom[, "chain"]
        if (is.null(resno2)) 
          resno2 = pdb$atom[, "resno"]
        if (is.null(b)) 
          b = pdb$atom[, "b"]
    } else {
      if (is.null(resno))
        resno = c(1:natom)
      if (is.null(resno2))
        resno2 = c(1:natom)
      if (is.null(resid))
        resid = rep("ALA", natom)
      if (is.null(eleno))
        eleno = c(1:natom)
      if (is.null(elety))
        elety = rep("CA", natom)
      if (is.null(segid))
        segid = rep("seg", natom)
      if (is.null(b))
        b = rep("0.00", natom)
    }
    if (length(as.vector(xyz))%%3 != 0) {
        stop("write.crd: 'length(xyz)' must be divisable by 3.")
    }
    check.lengths <- sum(length(resno), length(resid), length(eleno),
                         length(elety), length(resno2))
    if (check.lengths%%natom != 0) {
      stop("write.crd: the lengths of all input vectors != 'length(xyz)/3'.")
    }
    b <- as.numeric(b)
    eleno <- as.character(eleno)
    resno <- as.character(resno)
    resno2 <- as.character(resno2)
    
    crd.print <- function(eleno, resno, resid, elety,
                          x, y, z, segid="seg", resno2, b="0.00") {
      format <- "%5s%5s%4s  %-4s%10.5f%10.5f%10.5f %-4s %-4s%10.5f"
      sprintf(format, eleno, resno, resid, elety, x, y, z, segid, resno2, b)
    }
    coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, byrow = TRUE)
    if (verbose)
      cat(paste("Writing CRD file with", natom, "atoms "))
    
    lines <- NULL
    for (i in 1:natom) {
      lines <- rbind(lines, crd.print(
                    eleno = eleno[i], resno = resno[i], resid = resid[i],
                    elety = elety[i], 
                    x = coords[i,1], y = coords[i, 2], z = coords[i, 3],
                    segid = segid[i], resno2 = resno2[i], b = b[i]))
    }
    cat("* CRD from bio3d", file = file, "\n")
    cat("*", file = file, "\n", append = TRUE)
    cat(sprintf("%5g", natom), file = file, "\n", append = TRUE)
    cat(lines, file = file, sep = "\n", append = TRUE)
    if (verbose)
      cat(" DONE", "\n")
}

