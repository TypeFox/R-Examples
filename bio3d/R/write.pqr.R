`write.pqr` <-
function (pdb = NULL,
          xyz = pdb$xyz,
          resno = NULL,
          resid = NULL,
          eleno = NULL,
          elety = NULL,
          chain = NULL,
          o = NULL,
          b = NULL,
          het = FALSE,
          append = FALSE,
          verbose =FALSE,
          chainter = FALSE,
          file = "R.pdb") {

## For Testing:
## resno = NULL; resid = NULL;eleno = NULL;elety = NULL;chain = NULL;o = NULL;b = NULL; het = FALSE;append = FALSE;verbose =FALSE;chainter = FALSE
## pdb=mt; eleno=eleno.new; resno=resno.new; file="t3.pqr"
  
  if(is.null(xyz) || !is.numeric(xyz))
    stop("write.pqr: please provide a 'pdb' object or numeric 'xyz' vector")

  if(any(is.na(xyz)))
    stop("write.pqr: 'xyz' coordinates must have no NA's.")
  
  if(is.vector(xyz)) {
    natom <- length(xyz)/3
    nfile <- 1
  } else if (is.matrix(xyz)) {
    stop("write.pqr: no multimodel PQR support")
    ##natom <- ncol(xyz)/3 
    ##nfile <- nrow(xyz)
    ##if (verbose) {
    ##  cat("Multiple 'xyz' rows will be interperted as multimodels/frames\n")
    ##}
  } else {
    stop("write.pdb: 'xyz' or 'pdb$xyz' must be either a vector or matrix")
  }

  card <- rep("ATOM", natom)
  
  if(!is.null(pdb)) {
    if(natom == 1)
      ## make sure we are a matrix
      pdb$atom <- t(as.matrix(pdb$atom))
    
    if (het) 
      card <- c( rep("ATOM", nrow(pdb$atom)), rep("HETATM", nrow(pdb$het)) )
    if (is.null(resno)) {
      resno = pdb$atom[, "resno"]
      if (het) { resno = c(resno, pdb$het[, "resno"]) }}
    
    if (is.null(resid)) {
      resid = pdb$atom[, "resid"]
      if (het) { resid = c(resid, pdb$het[, "resid"]) }}
    
    if (is.null(eleno)) {
      eleno = pdb$atom[, "eleno"]
      if (het) { eleno = c(eleno, pdb$het[, "eleno"]) }}
    
    if (is.null(elety)) {
      elety = pdb$atom[, "elety"]
      if (het) { elety = c(elety, pdb$het[, "elety"]) }}
    
    if (is.null(chain)) {
      chain = pdb$atom[, "chain"]
      if (het) { chain = c(chain, pdb$het[, "chain"]) }}
    
    if (is.null(o)) {
      o = pdb$atom[, "o"]
      if (het) { o = c(o, pdb$het[, "o"]) }}
    
    if (is.null(b)) {
      b = pdb$atom[, "b"]
      if (het) { b = c(b, pdb$het[, "b"]) }}
    
    if (any(is.na(o))) {      o = rep("1.00", natom) }
    if (any(is.na(b))) {      b = rep("0.00", natom) }
    #if (any(is.na(chain))) { chain = rep(" ", natom) }
    chain[is.na(chain)]= " "
    
  } else {
    if (is.null(resno)) resno = c(1:natom)
    if (is.null(resid)) resid = rep("ALA", natom)
    if (is.null(eleno)) eleno = c(1:natom)
    if (is.null(elety)) elety = rep("CA", natom)
    if (is.null(chain)) chain = rep(" ", natom)
    if (is.null(o))         o = rep("1.00",natom)
    if (is.null(b))         b = rep("0.00", natom)
  }

  
  if (!is.logical(append)) 
    stop("write.pqr: 'append' must be logical TRUE/FALSE")
  
  if (length(as.vector(xyz))%%3 != 0) {
    stop("write.pqr: 'length(xyz)' must be divisable by 3.")
  }
  check.lengths <- sum(length(resno), length(resid), length(eleno),
                       length(elety), length(o), length(b))
  if (check.lengths%%natom != 0) {
    stop("write.pqr: the lengths of all input vectors != 'length(xyz)/3'.")
  }

  
  o <- as.numeric(o)
  b <- as.numeric(b)
  eleno <- as.character(eleno)
  resno <- as.character(resno)
  ## Inserted Jul 8th 2008 for adding TER between chains
  ter.lines <- (which(!duplicated(chain))[-1] - 1)

####  
####  ## Edit: Sat Aug  1 14:48:48 PDT 2009
####  ## for speed imporvment and for 
####  ## implementing 6 character atom numbers 
####
  
  if(nfile==1) {
    coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, byrow = TRUE)
    if (verbose) {
      cat(paste("Writing 1 frame with",natom,"atoms "))
    }
    
    coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, byrow = TRUE)
    lines <- matrix(, ncol=1, nrow=natom)

    ## Four format otions: regular; elety > 3; eleno > 5; eleno > 5 & elety > 3
    ## cases  nchar(elety) > 3; nchar(eleno) > 5
    cases <- matrix(1,ncol=2,nrow=natom)
    cases[(nchar(eleno) > 5) ,1] = 3
    cases[(nchar(elety) < 4) ,2] = 0
    cases <- rowSums(cases)
    
    ind.1 <- which(cases==1)
    ind.2 <- which(cases==2)
    ind.3 <- which(cases==3)
    ind.4 <- which(cases==4)

    atom.print.1 <- function(card = "ATOM", eleno, elety, alt = "",
        resid, chain = "", resno, insert = "", x, y, z, o = "1.00",
        b = "0.00", segid = "") {

      format <- "%-6s%5s  %-3s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%8.4f%7.4f%6s%4s"
      sprintf(format, card, eleno, elety, alt, resid, chain,
              resno, insert, "", x, y, z, o, b, "", segid)
    } 

    atom.print.2 <- function(card = "ATOM", eleno, elety, alt = "",
        resid, chain = "", resno, insert = "", x, y, z, o = "1.00",
        b = "0.00", segid = "") {

      format <- "%-6s%5s %-4s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%8.4f%7.4f%6s%4s"
      sprintf(format, card, eleno, elety, alt, resid, chain,
              resno, insert, "", x, y, z, o, b, "", segid)
    } 

    atom.print.3 <- function(card = "ATOM", eleno, elety, alt = "",
        resid, chain = "", resno, insert = "", x, y, z, o = "1.00",
        b = "0.00", segid = "") {

      format <- "%-4s%7s  %-3s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%8.4f%7.4f%6s%4s"
      sprintf(format, card, eleno, elety, alt, resid, chain,
              resno, insert, "", x, y, z, o, b, "", segid)
    } 
  
    atom.print.4 <- function(card = "ATOM", eleno, elety, alt = "",
        resid, chain = "", resno, insert = "", x, y, z, o = "1.00",
        b = "0.00", segid = "") {

      format <- "%-4s%7s %-4s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%8.4f%7.4f%6s%4s"
      sprintf(format, card, eleno, elety, alt, resid, chain,
              resno, insert, "", x, y, z, o, b, "", segid)
    } 


  
    if(length(ind.1)>0) {
      lines[ind.1,] <- atom.print.1( card = card[ind.1],
                                    eleno = eleno[ind.1],
                                    elety = elety[ind.1],
                                    resid = resid[ind.1],
                                    chain = chain[ind.1],
                                    resno = resno[ind.1],
                                    x = coords[ind.1, 1],
                                    y = coords[ind.1, 2],
                                    z = coords[ind.1, 3],
                                    o = o[ind.1], b = b[ind.1] )
    }
    if(length(ind.2)>0) {
      lines[ind.2,] <- atom.print.2( card = card[ind.2],
                                    eleno = eleno[ind.2],
                                    elety = elety[ind.2],
                                    resid = resid[ind.2],
                                    chain = chain[ind.2],
                                    resno = resno[ind.2],
                                    x = coords[ind.2, 1],
                                    y = coords[ind.2, 2],
                                    z = coords[ind.2, 3],
                                    o = o[ind.2], b = b[ind.2] )
    }
    if(length(ind.3)>0) {
      lines[ind.3,] <- atom.print.3( card = card[ind.3],
                                    eleno = eleno[ind.3],
                                    elety = elety[ind.3],
                                    resid = resid[ind.3],
                                    chain = chain[ind.3],
                                    resno = resno[ind.3],
                                    x = coords[ind.3, 1],
                                    y = coords[ind.3, 2],
                                    z = coords[ind.3, 3],
                                    o = o[ind.3], b = b[ind.3] )
    }
    if(length(ind.4)>0) {
      lines[ind.4,] <- atom.print.4( card = card[ind.4],
                                    eleno = eleno[ind.4],
                                    elety = elety[ind.4],
                                    resid = resid[ind.4],
                                    chain = chain[ind.4],
                                    resno = resno[ind.4],
                                    x = coords[ind.4, 1],
                                    y = coords[ind.4, 2],
                                    z = coords[ind.4, 3],
                                    o = o[ind.4], b = b[ind.4] )
    }
    
    write.table(lines, file=file, quote=FALSE,
                row.names=FALSE, col.names=FALSE, append=append)
    
####
#### End of Edit: removed big chunks of old code
####    
    
  } else {
    if (verbose) {
      cat(paste("Writing",nfile,"frames with",natom,"atoms"),"\n")
      cat("Frame Progress (x50) ")
    }
    stop("REMOVED code for multimodel PQR as these files dont have much support")
  }
  if (verbose) cat(" DONE","\n")
}

