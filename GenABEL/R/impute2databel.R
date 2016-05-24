#'
#' converts IMPUTE-imputed files to DatABEL (filevector) format
#'
#' this function converts IMPUTE-imputed files to DatABEL (filevector) format
#' containing estimated dosages.
#' After conversion, two files (outfile.fvi and outfile.fvd), corresponding
#' to a single filevector object, will appear on the disk; a 'databel-class'
#' object connected to these files will be returned to R.
#'
#' @param genofile IMPUTE genotype file name
#' @param samplefile IMPUTE sample file name
#' @param outfile output file name
#' @param makeprob whether probability-files are also to be arranged
#' @param old for developers' use only
#' @param dataOutType the output data type, either "FLOAT" or "DOUBLE" (or
#'        another DatABEL/filevector type)
#'
#' @return 'databel-class' object
#'
#'


impute2databel <- function(genofile, samplefile, outfile,
                           makeprob = TRUE, old = FALSE, dataOutType = "FLOAT")
{
  if (!require(DatABEL))
    stop("this function requires the DatABEL package to be installed")

  if (missing(genofile))
    stop("genofile file must be specified")

  if (missing(outfile))
    outfile <- genofile

  if (dataOutType != "FLOAT")
      warning("The non-float dataOutType is not fully supported; your outputs may be in 'FLOAT'...",
              immediate. = TRUE);

  ## extract snp names (varnames)
  ## if (tmpname != "")
  ##     text2databel(infile=genofile,outfile=outfile,
  ##            colnames=tmpname,
  ##            rownames=2,skipcols=5,
  ##            #skiprows,
  ##            transpose=TRUE,R_matrix=FALSE,type = dataOutType)
  ## else
  tmpname <- get_temporary_file_name()
  tmp_fv  <- text2databel(infile=genofile,
                          outfile=tmpname,
                          rownames=2,
                          skipcols=5,
                          ## skiprows,
                          transpose=TRUE,
                          R_matrix=FALSE,
                          type = dataOutType)
  if ((dim(tmp_fv)[1] %% 3) != 0)
    stop("strange data in IMPUTE geno file: ",
         "number of columns - 5 not dividable by 3")

  makedose <- function(prob) {
    dose <- 2 * prob[c(FALSE, FALSE, TRUE)] + prob[c(FALSE, TRUE, FALSE)]
    bp   <- prob[c(TRUE, FALSE, FALSE)]
    miss <- which(abs(bp) < 1e-16 & abs(dose) < 1e-16)
    if (length(miss) > 0 ) dose[miss] <- NA
    return(dose)
  }

  pfun <- function(a) return(a)

  ## dosefile <- make_empty_fvf(paste(outfile, ".dose", sep=""),
  ##                            nvariables = dim(tmp_fv)[2],
  ##                            nobservations = round(dim(tmp_fv)[1]/3),
  ##                            type = dataOutType)

  ## print(dimnames(tmp_fv)[[2]])
  ## print(dim(tmp_fv))
  ## print(class(tmp_fv))
  saved_names <- get_dimnames(tmp_fv)
  varnames    <- saved_names[[2]]
  ## print("original names of tmp_fv")
  ## print(saved_names)
  ## print(class(tmp_fv))

  ## tmplst <- list( as.character(1:dim(tmp_fv)[1]) , as.character(1:dim(tmp_fv)[2]) )
  ## print(class(tmp_fv))
  ## print(dimnames(tmp_fv)[[2]])
  ## set_dimnames(tmp_fv) <- tmplst
  ## print(class(tmp_fv))
  ## print(dimnames(tmp_fv)[[2]])

  ## print("before apply2dfo")
  ## print("calling apply2dfo")
  if (old) {                             #
    ##          stop("'old' ways are no longer there...")
    dosefile <- apply2dfo(dfodata=tmp_fv,
                          anFUN = "makedose",
                          MAR = 2,
                          procFUN = "pfun",
                          prob=get("SNP"),
                          outclass="databel",
                          outfile=paste(outfile, ".dose", sep=""),
                          type=dataOutType,
                          transpose=FALSE)
    if (makeprob) {
      stop("makeprob is not possible with 'old' style")
    }
  } else {
    res <- .Call("iterator", tmp_fv@data, as.integer(0), as.integer(0),
                 as.character("databel_impute_prob_2_databel_mach_dose"),
                 paste(outfile, ".dose", sep=""), as.integer(2),
                 as.integer(0))
    dosefile <- databel(paste(outfile, ".dose", sep=""), 64, readonly=FALSE)
    if (makeprob) {
      res <- .Call("iterator", tmp_fv@data, as.integer(0), as.integer(0),
                   as.character("databel_impute_prob_2_databel_mach_prob"),
                   paste(outfile, ".prob", sep=""), as.integer(2),
                   as.integer(0))
      probfile <- databel(paste(outfile, ".prob", sep=""), 64, readonly=FALSE)
    }

    ## res <- .Call("databel_impute_prob_2_databel_mach_dose",
    ## tmp_fv@data, paste(outfile,".dose",sep=""), as.integer(64))
    ## dosefile <- databel_filtered_R(paste(outfile,".dose",sep=""),64)
    ## if (makeprob) {
    ##   res <- .Call("databel_impute_prob_2_databel_mach_prob",
    ##   tmp_fv@data, paste(outfile,".prob",sep=""), as.integer(64))
    ##   #print("AAAA")
    ##   probfile <- databel_filtered_R(paste(outfile, ".prob", sep=""),64)
    ##   #print(class(probfile))
    ##   #print(dim(probfile))
    ##   #print(dimnames(probfile))
    ##   #print(get_dimnames(probfile))
    ##}
  }                                #
  ##!print("after iterator")

  set_dimnames(dosefile) <- list(get_dimnames(dosefile)[[1]], varnames)
  ## print("dimnames dose:")
  ## print(get_dimnames(dosefile))


  if (!missing(samplefile)) {
    l1 <- scan(samplefile, what="character", nlines=1)
    l2 <- scan(samplefile, what="character", nlines=1, skip=1)
    temp <- scan(samplefile, what="character", nlines=1, skip=2)

    if( (length(l1)!=length(temp)) || (length(temp)!=length(l2)) ) {
      warning("Wrong number of words in sample file header.");
      }

    samnames <- scan(samplefile, what="character", skip=2)
    samnames <- samnames[c(FALSE, TRUE, rep(FALSE, (length(temp)-2)))]
    if (length(samnames) == dim(dosefile)[1]) {
      ##!print("before dose-file")
      ##print(samnames)
      ##print(varnames)
      set_dimnames(dosefile) <- list(samnames, varnames)
      ##!print("after dose-file")
      ##print(get_dimnames(dosefile))
    } else {
      warning("number of IDs specified in sample file does not match to ",
              "geno-dimension; dropping ID names")
    }
    if (makeprob) {
      if (length(samnames) == dim(probfile)[1] &&
          2*length(varnames)==dim(probfile)[2]) {
        dpnames <- rep(0, length(varnames)*2)
        dpnames[c(TRUE, FALSE)] <- paste(varnames, "_11", sep="")
        dpnames[c(FALSE, TRUE)] <- paste(varnames, "_01", sep="")
        ##!print("before prob-file")
        ##print(samnames)
        ##print(dpnames)
        set_dimnames(probfile) <- list(samnames, dpnames) #get_dimnames(probfile)[[2]])
        ##!print("after prob-file")
        ##print(get_dimnames(probfile))
      } else {
        warning("number of IDs/SNPs specified in sample file does not match to ",
                "geno-dimension; dropping IDs/SNPs names")
      }
    }
  } else
    warning("sample file not specified, you will not be able to use ID names ",
            "(only index)")

  rm(tmp_fv); gc();
  unlink(paste(tmpname, "*", sep=""))

  ##disconnect(dosefile)
  ##connect(dosefile)
  ##!print("at exit dosefile:")
  ##print(get_dimnames(dosefile))
  if (makeprob) {
    ##print("at exit probfile:")
    ##print(get_dimnames(probfile))
    disconnect(probfile)
  }
  ##!print("at return")
  return(dosefile)
}
