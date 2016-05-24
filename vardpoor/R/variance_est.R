
variance_est <- function(Y, H, PSU, w_final, N_h=NULL, fh_zero=FALSE, PSU_level=TRUE, period=NULL, dataset=NULL) {

  ### Checking

    if (!is.logical(fh_zero)) stop("'fh_zero' must be the logical value")
    if (!is.logical(PSU_level)) stop("'PSU_level' must be the logical value")

    if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE] 
      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }
      if(!is.null(PSU)) {
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) PSU <- dataset[, PSU, with=FALSE] }
      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }
       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE] }
      }

  # Y
  Y <- data.table(Y, check.names=TRUE)
  n <- nrow(Y)
  m <- ncol(Y)
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numeric values")
  if (any(is.na(Y))) print("'Y' has unknown values")
  if (is.null(names(Y))) stop("'Y' must be the column names")
  
  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' length must be equal with 'Y' row count")
  if (ncol(H) != 1) stop("'H' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colnames")
  H[, (names(H)):=lapply(.SD, as.character)]
  
  # PSU
  PSU <- data.table(PSU)
  if (nrow(PSU) != n) stop("'PSU' length must be equal with 'Y' row count")
  if (ncol(PSU) != 1) stop("'PSU' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(PSU))) stop("'PSU' has unknown values")
  PSU[, (names(PSU)):=lapply(.SD, as.character)]
  
  # w_final
  w_final <- data.frame(w_final)
  if (nrow(w_final) != n) stop("'w_final' must be equal with 'Y' row count")
  if (ncol(w_final) != 1) stop("'w_final' must be vector or 1 column data.frame, matrix, data.table")
  w_final <- w_final[, 1]
  if (!is.numeric(w_final)) stop("'w_final' must be numerical")
  if (any(is.na(w_final))) stop("'w_final' has unknown values")

  # period     
  if (!is.null(period)) {
       period <- data.table(period)
       if (any(duplicated(names(period)))) 
                 stop("'period' are duplicate column names: ", 
                      paste(names(period)[duplicated(names(period))], collapse = ","))
       if (nrow(period) != n) stop("'period' must be the same length as 'inc'")
       if(any(is.na(period))) stop("'period' has unknown values")  
  }   
  np <- sum(ncol(period))
  
  # N_h
  if (!is.null(N_h)) {
      N_h <- data.table(N_h)
      if (ncol(N_h) != np+2) stop(paste0("'N_h' should be ",toString(np+2)," columns"))
      if (!is.numeric(N_h[[ncol(N_h)]])) stop("The last column of 'N_h' should be numerical")
      if (any(is.na(N_h))) stop("'N_h' has unknown values") 
      if (is.null(names(N_h))) stop("'N_h' must be colnames")
      if (all(names(H) %in% names(N_h))) {N_h[, (names(H)):=lapply(.SD, as.character), .SDcols=names(H)]
             } else stop("All strata titles of 'H' have not in 'N_h'")
      if (is.null(period)) {
             if (names(H) != names(N_h)[1]) stop("Strata titles for 'H' and 'N_h' is not equal")
             if (any(is.na(merge(unique(H), N_h, by=names(H), all.x=TRUE)))) stop("'N_h' is not defined for all stratas")
             if (any(duplicated(N_h[, head(names(N_h),-1), with=FALSE]))) stop("Strata values for 'N_h' must be unique")
       } else { pH <- data.table(period, H)
                if (any(names(pH) != names(N_h)[c(1:(1+np))])) stop("Strata titles for 'period' with 'H' and 'N_h' is not equal")
                nperH <- names(period)
                if (pH[, class(get(nperH))]!=N_h[, class(get(nperH))]) 
                                                       stop("Period class for 'period' and 'N_h' is not equal ")
                if (any(is.na(merge(unique(pH), N_h, by=names(pH), all.x=TRUE)))) stop("'N_h' is not defined for all stratas and periods")
                if (any(duplicated(N_h[, head(names(N_h),-1), with=FALSE]))) stop("Strata values for 'N_h' must be unique in all periods")
                }
    setnames(N_h, names(N_h)[ncol(N_h)], "N_h")
    setkeyv(N_h, names(N_h)[c(1:(1+np))])
  } else {
    Nh <- data.table(H, w_final)
    if (!is.null(period)) Nh <- data.table(period, Nh)
    N_h <- Nh[, .(N_h = sum(w_final, na.rm=TRUE)), keyby=c(names(Nh)[1:(1+np)])]
  }
  pH <- NULL  

  ### Calculation
  
  # z_hi
  .SD <- .N <- NULL
  hpY <- data.table(H, PSU, Y*w_final)
  if (!is.null(period)) hpY <- data.table(period, hpY)
  z_hi <- hpY[, lapply(.SD, sum, na.rm=TRUE), 
                       keyby = c(names(hpY)[1:(2+np)]),
                       .SDcols = names(hpY)[-(1:(2+np))]]
  setkeyv(z_hi, names(z_hi)[c(1:(1+np))]) 

  # n_h
  n_h <- data.table(z_hi[, c(1:(1+np)), with=FALSE])
  n_h <- n_h[,.N, keyby=c(names(n_h)[1:(1+np)])]
  setnames(n_h, "N", "n_h")

  # var_z_hi
  var_z_hi <- z_hi[, lapply(.SD, var, na.rm=FALSE), keyby = c(names(z_hi)[1:(1+np)]),
    .SDcols = names(z_hi)[-(1:(2+np))]]
  
  # f_h
  F_h <- merge(N_h, n_h, by = names(hpY)[c(1:(1+np))], sort=TRUE)
  F_h[, f_h:=n_h/N_h]

  if (nrow(F_h[n_h==1 & f_h != 1])>0) {
    print("There is stratas, where n_h == 1 and f_h <> 1")
    print("Not possible to estimate the variance in these stratas!")
    print("At these stratas estimation of variance was not calculated")
    nh <- F_h[n_h==1 & f_h != 1]
    print(nh)
  }
  
  if (nrow(F_h[f_h > 1])>0) {    
     print("There is stratas, where f_h > 1")
     print("At these stratas estimation of variance will be 0")
     print(F_h[f_h > 1])
     F_h[f_h > 1, f_h:=1]
   }

  # fh1
  if (!(PSU_level)) {
       n_h1 <- NULL
       fh1 <- data.table(hpY[, c(1:(1+np)), with=FALSE], w_final)
       fh1[, n_h1:=1]
       fh1 <- fh1[, lapply(.SD, sum, na.rm=TRUE), keyby = c(names(fh1)[c(1:(1+np))]),
	           .SDcols = c("n_h1", "w_final")]
       fh1[, fh1:=n_h1/w_final]
       f_h <- fh1[, "fh1", with=FALSE]
     } else  f_h <- F_h[,"f_h", with=FALSE]

  # var_h
  var_h <- data.table(matrix((1 - f_h*(1-fh_zero)) * n_h$n_h) * var_z_hi[, -c(1:(1+np)), with=FALSE])
  var_h2 <- copy(var_h)
  if (np>0) var_h2 <- data.table(var_z_hi[, c(1:np), with=FALSE], var_h)

  # Variance_est 

  if (np==0) {var_est <- data.table(t(colSums(var_h, na.rm=TRUE)))
           } else   var_est <- var_h2[, lapply(.SD, sum, na.rm=TRUE), 
                                        keyby = c(names(var_h2)[c(1:np)]),
                                       .SDcols = names(var_h)]


  return(var_est)
}



