
variance_othstr <- function(Y, H, H2, w_final, N_h=NULL, N_h2, period=NULL, dataset=NULL) {

  ### Checking
    if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE] 


      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }
      if(!is.null(H2)) {
          if (min(H2 %in% names(dataset))!=1) stop("'H2' does not exist in 'dataset'!")
          if (min(H2 %in% names(dataset))==1) H2 <- dataset[, H2, with=FALSE] }


      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }


       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE]}
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


  # H2
  H2 <- data.table(H2)
  if (nrow(H2) != n) stop("'H2' length must be equal with 'Y' row count")
  if (ncol(H2) != 1) stop("'H2' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H2))) stop("'H2' has unknown values")
  if (is.null(names(H2))) stop("'H2' must be colnames")
  H2[, (names(H2)):=lapply(.SD, as.character)]


  # w_final
  w_final <- data.frame(w_final)
  if (nrow(w_final) != n) stop("'w_final' must be equal with 'Y' row count")
  if (ncol(w_final) != 1) stop("'w_final' must be vector or 1 column data.frame, matrix, data.table")
  w_final <- w_final[,1]
  if (!is.numeric(w_final)) stop("'w_final' must be numerical")
  if (any(is.na(w_final))) stop("'w_final' has unknown values")


  # period     
  if (!is.null(period)) {
       period <- data.table(period)
       if (any(duplicated(names(period)))) 
                 stop("'period' are duplicate column names: ", 
                      paste(names(period)[duplicated(names(period))], collapse = ","))
       if (nrow(period) != n) stop("'period' must be the same length as 'Y'")
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
             if (any(is.na(merge(unique(H), N_h, by=names(H), all.x = TRUE)))) stop("'N_h' is not defined for all stratas")
       } else { pH <- data.table(period, H)
                if (any(names(pH) != names(N_h)[c(1:(1+np))])) stop("Strata titles for 'period' with 'H' and 'N_h' is not equal")
                nperH <- names(period)
                if (pH[, class(get(nperH))]!=N_h[, class(get(nperH))])  stop("Period class for 'period' and 'N_h' is not equal")
                if (any(is.na(merge(unique(pH), N_h, by=names(pH), all.x = TRUE)))) stop("'N_h' is not defined for all stratas and periods")
               } 
     setkeyv(N_h, names(N_h)[c(1:(1+np))])
  } else {
    Nh <- data.table(H, w_final)
    if (!is.null(period)) Nh <- data.table(period, Nh)
    N_h <- Nh[, .(N_h = sum(w_final, na.rm=TRUE)), keyby=c(names(Nh)[1:(1+np)])]
  }
  Nh1 <- names(N_h)[ncol(N_h)]


  # N_h2
  if (!is.null(N_h2)) {
      N_h2 <- data.table(N_h2)
      if (ncol(N_h2) != np+2) stop(paste0("'N_h2' should be ",toString(np+2)," columns"))
      if (!is.numeric(N_h2[[ncol(N_h2)]])) stop("The last column of 'N_h2' should be numerical")
      if (any(is.na(N_h2))) stop("'N_h2' has unknown values") 
      if (is.null(names(N_h2))) stop("'N_h2' must be colnames")
      if (all(names(H2) %in% names(N_h2))) {N_h2[, (names(H2)):=lapply(.SD, as.character), .SDcols=names(H2)]
             } else stop("All strata titles of 'H2' have not in 'N_h2'")
      if (is.null(period)) {
             if (names(H2) != names(N_h2)[1]) stop("Strata titles for 'H2' and 'N_h2' is not equal")
             if (any(is.na(merge(unique(H2), N_h2, by=names(H2), all.x = TRUE)))) stop("'N_h2' is not defined for all stratas")
       } else { pH2 <- data.table(period, H2)
                if (any(names(pH2) != names(N_h2)[c(1:(1+np))])) stop("Strata titles for 'period' with 'H2' and 'N_h2' is not equal")
                nperH <- names(period)
                if (pH2[, class(get(nperH))]!=N_h2[, class(get(nperH))])  stop("Period class for 'period' and 'N_h2' is not equal")
                if (any(is.na(merge(unique(pH2), N_h2, by=names(pH2), all.x = TRUE)))) stop("'N_h2' is not defined for all stratas and periods")
                } 
    setkeyv(N_h2, names(N_h2)[c(1:(1+np))])
  } else stop ("N_h2 is not defined!")
  Nh2 <- names(N_h2)[ncol(N_h2)]

  if (all(names(H)==names(H2))) {
      if (!is.null(N_h2))  setnames(N_h2, names(H), paste0(names(H), "2"))  
      setnames(H2, names(H2), paste0(names(H), "2"))  }



  ### Calculation
  
  # z_hi
  f_h1 <- .SD <- .N <- NULL
  Ys <- copy(Y)
  Ys[, paste0(names(Y),"_sa") := lapply(Y, function(x) w_final * x^2)]
  Ys[, paste0(names(Y),"_sb") := lapply(Y, function(x) x * w_final)]
  Ys[, paste0(names(Y),"_sc") := lapply(Y, function(x) x ^ 2)]
  Ys[, paste0(names(Y),"_sd") := Y]


  Ys <- data.table(H, H2, Ys)
  if (!is.null(period)) Ys <- data.table(period, Ys)


  # n_h1
  n_h1 <- data.table(H)
  if (!is.null(period))   n_h1 <- data.table(period, n_h1)
  n_h1 <- n_h1[, .(n_h1=.N), keyby=c(names(n_h1))]

  F_h1 <- merge(N_h, n_h1, keyby = c(names(N_h)[1:(1+np)]))
  F_h1[, f_h1:=n_h1/get(Nh1)]

  if (nrow(F_h1[n_h1==1 & f_h1 != 1])>0) {
    print("There is stratas, where n_h1 == 1 and f_h1 <> 1")
    print("Not possible to estimate the variance in these stratas!")
    print("At these stratas estimation of variance was not calculated")
    nh1 <- F_h1[n_h1==1 & f_h1 != 1]
    print(nh1)
  }

  # n_h2
  n_h2 <- data.table(H2)
  if (!is.null(period)) n_h2 <- data.table(period, n_h2)
  nn_h2 <- names(n_h2)
  n_h2 <- n_h2[, .(n_h2=.N), keyby=nn_h2]

  F_h2 <- merge(N_h2, n_h2, keyby=nn_h2)
  F_h2[, f_h2:=n_h2/get(Nh2)]

  if (nrow(F_h2[n_h2==1 & f_h2 != 1])>0) {
    print("There is stratas, where n_h2 == 1 and f_h2 <> 1")
    print("Not possible to estimate the variance in these stratas!")
    print("At these stratas estimation of variance was not calculated")
    nh2 <- F_h2[n_h2==1 & f_h2 != 1]
    print(nh2)
  }
  
  if (nrow(F_h2[f_h2 > 1])>0) {    
      print("There is stratas, where f_h2 > 1")
      print("At these stratas estimation of variance will be 0")
      print(F_h2[f_h2 > 1])
      F_h2[f_h2 > 1, f_h2:=1]
   }

  z_h_h2 <- Ys[, lapply(.SD, sum, na.rm=TRUE), keyby = c(names(Ys)[1:(2+np)]),
                      .SDcols = names(Ys)[-(0:(ncol(Y)+2+np))]]

  z_h_h2 <- merge(z_h_h2, F_h1, keyby = names(z_h_h2)[c(1:(1+np))])

  pop <- z_h_h2[[Nh1]]

  z_h_h2[, paste0(names(Y),"_sc") := lapply(.SD[, 
           paste0(names(Y),"_sc"), with=FALSE], function(x)
           x * pop ^ 2 * ( 1/n_h1 - 1/pop)/(n_h1-1))]


  z_h_h2[, paste0(names(Y),"_sd") := lapply(.SD[,
           paste0(names(Y),"_sd"), with=FALSE], function(x) 
                 (1/n_h1) * x ^ 2  * pop ^ 2 * (1/n_h1 - 1/pop)/(n_h1-1))]


  z_h_h2[n_h1==1, paste0(names(Y),"_sc"):=NA]
  z_h_h2[n_h1==1, paste0(names(Y),"_sd"):=NA]


  nameszh2 <- names(H2)
  if (!is.null(period)) nameszh2 <- c(names(period), nameszh2)
  
  zh2 <- z_h_h2[, lapply(.SD, sum, na.rm=TRUE), keyby = nameszh2,
                      .SDcols = names(z_h_h2)[-(1:(2+np))]] 


  zh2 <- merge(zh2, F_h2, by=nn_h2)
  pop2 <- zh2[[names(N_h2)[ncol(N_h2)]]]
  nh2 <- zh2[["n_h2"]]
  f_h2 <- zh2[["f_h2"]]

  # s2
  s2_g <- zh2[,mapply(function(sa, sb, sc, sd) sa/(pop2-1)-pop2/(pop2-1)*((sb/pop2)^2-(sc-sd)/pop2^2),
              zh2[, paste0(names(Y),"_sa"), with=FALSE], 
              zh2[, paste0(names(Y),"_sb"), with=FALSE],
              zh2[, paste0(names(Y),"_sc"), with=FALSE],
              zh2[, paste0(names(Y),"_sd"), with=FALSE])]


  # var_g 
  if (is.null(nrow(s2_g))) s2_g <- t(s2_g)
  s2_g <- data.table(s2_g)
  setnames(s2_g, names(s2_g), names(Y))

  s2g <- data.table(zh2[, nn_h2, with=FALSE], s2_g)

  s2_g <- matrix(pop2^2 * 1/nh2 * (1 - f_h2)) * s2_g

  if (np>0) s2_g <- data.table(zh2[, names(period), with=FALSE], s2_g)

  # Variance_est
  if (np==0) {var_est <- data.table(t(colSums(s2_g, na.rm=TRUE)))
           } else var_est <- s2_g[, lapply(.SD, sum, na.rm=TRUE), 
                                       keyby = c(names(s2_g)[c(1:np)]),
                                      .SDcols = names(Y)]
  list(s2g=s2g,
       var_est=var_est)
}