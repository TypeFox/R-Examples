"dccm.enma" <- function(x, ncore=NULL, na.rm=FALSE, ...) {
 
  enma <- x
  if(!inherits(enma, "enma"))
    stop("input should be an 'enma' object as obtained from 'nma.pdbs'")

  ## Parallelized by parallel package
  ncore <- setup.ncore(ncore, bigmem = FALSE)

  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  mass <- TRUE
  if(!is.null(enma$call$mass))
    mass <- enma$call$mass
  
  pi <- 3.14159265359
  dims <- dim(enma$U.subspace)
  if(is.null(enma$full.nma)) {
    if((dims[1]-6)>dims[2])
      warning(paste(dims[2], "modes used in the calculation of the DCCMs"))
  }
  
  myCalcDCCM <- function(i, enma, na.rm=FALSE) {
    if(is.null(enma$full.nma)) {
      if(mass) {
        freqs <- sqrt(abs(enma$L[i,])) / (2 * pi)
        fcs <- NULL
      }
      else {
        freqs <- NULL
        fcs <- enma$L[i,]
      }

      if(na.rm) {
        inds <- which( !is.na(enma$U.subspace[,1,i]) )
        U <- enma$U.subspace[inds,,i]
      }
      else
        U <- enma$U.subspace[,,i]
      dummy.nma <- list(U=U,
                        L=enma$L[i,],
                        modes=NULL,
                        frequencies=freqs,
                        force.constants=fcs,
                        triv.modes=0,
                        natoms=nrow(U)/3)
                        ##natoms=nrow(enma$U.subspace[,,i])/3)
      class(dummy.nma) <- "nma"

      invisible(capture.output(
        cm.tmp <- dccm.nma(dummy.nma, ncore=1) ))
    }
    else {
      invisible(capture.output(
        cm.tmp <- dccm.nma(enma$full.nma[[i]], ncore=1) ))
    }

    setTxtProgressBar(pb, i)
    return(cm.tmp)
  }

  ## do the calc
  pb <- txtProgressBar(min=1, max=dims[3L], style=3)
  all.dccm <- mylapply(1:dims[3L], myCalcDCCM, enma, na.rm=na.rm)
  close(pb)

  if(any(is.na(enma$U.subspace)))
    arr <- FALSE
  else
    arr <- TRUE
  
  if(arr) {
    ## convert to a 3d-array
    dccm.arr <- array(0, dim=c(dims[1L]/3, dims[1L]/3, dims[3L]))
  
    ## collect data
    for(i in 1:length(all.dccm)) {
      tmp.cm <- all.dccm[[i]]
      dccm.arr[,,i] <- tmp.cm
    }
  }
  
  if(arr) {
    avg <- apply(dccm.arr, 1:2, mean)
    class(avg) <- c("matrix", "dccm")
    out <- list(all.dccm=dccm.arr, avg.dccm=avg)
  }
  else {
    out <- list(all.dccm=all.dccm, avg.dccm=NULL)
  }
  
  
  return(out)
}
