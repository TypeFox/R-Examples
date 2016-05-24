naImpute <- function(FMLA,DATA,impfn=median,na.rm=TRUE) {
  if (!all("na.rm" %in% names(formals(impfn)))){stop("The imputation function requires a na.rm argument like that of mean.default() or median()")}
  fmla.rhs <- terms.formula(if (length(FMLA)>2) FMLA[-2] else FMLA,
                            data = DATA, keep.order=TRUE)
  dat <- get_all_vars(fmla.rhs,DATA)
  badfactor <- sapply(dat,function(x) nlevels(x)==1)
  dat[badfactor] <- lapply(dat[badfactor], as.integer) ##Is this right? shouldn't it be dat[,badfactor] ?
  factor.dat <- sapply(dat,is.factor)
  ordered.dat <- sapply(dat,is.ordered)
  dat.NA <- as.data.frame(lapply(dat, is.na))
  impute.dat <- sapply(dat.NA,any)
  dat.NA <- dat.NA[impute.dat& (!factor.dat | ordered.dat)]
  if (any(impute.dat & (!factor.dat | ordered.dat))) {
    names(dat.NA) <- paste(names(dat.NA), 'NA',sep='.')
  }
  if (any(impute.dat))  {
    dat <- lapply(dat, function(x){
      if (is.factor(x) & !is.ordered(x)) {
        if (any(is.na(x))) {
          levels(x) <- c(levels(x),'.NA') ##There is a weird problem here with factors with no missing data
          x[is.na(x)] <- '.NA'
        }
      } else {
        if (is.ordered(x)) {
          x[is.na(x)] <- levels(x)[1]
        } else {
          if (is.logical(x)) {
            x[is.na(x)] <- impfn(x,na.rm=na.rm)>.5
          } else {
            x[is.na(x)] <- impfn(x,na.rm=na.rm)
          }
        }
      }
      x
    }
    )
    dat <- as.data.frame(dat)
    dat <- data.frame(DATA[setdiff(names(DATA),
                                   c(names(dat),names(dat.NA)))],
                      dat, dat.NA)
    TFMLA <- if (length(dat.NA)) update.formula(FMLA, as.formula(paste(".~.+",paste(names(dat.NA), collapse=" + ")))) else FMLA
    TFMLA <- terms.formula(TFMLA,data = dat, keep.order=TRUE)
    return(structure(dat, terms=TFMLA))
  } else return(structure(DATA,terms=terms.formula(FMLA,data = DATA, keep.order=TRUE)))
}
