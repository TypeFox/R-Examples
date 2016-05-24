 vardom_othstr <- function(Y, H, H2, PSU, w_final,
                   id = NULL,  
                   Dom = NULL,
                   period = NULL,
                   N_h = NULL,
                   N_h2 = NULL,
                   Z = NULL,
                   X = NULL,
                   g = NULL,
                   q = NULL,
                   dataset = NULL, 
                   confidence = .95, 
                   percentratio = 1,
                   outp_lin = FALSE,
                   outp_res = FALSE) {
 
  ### Checking

  if (length(outp_lin) != 1 | !any(is.logical(outp_lin))) stop("'outp_lin' must be the logical value")
  if (length(outp_res) != 1 | !any(is.logical(outp_res))) stop("'outp_res' must be the logical value")
  if (length(percentratio) != 1 | !any(is.integer(percentratio) | percentratio > 0)) stop("'percentratio' must be the positive integer value")
  if(length(confidence) != 1 | any(!is.numeric(confidence) | confidence < 0 | confidence > 1)) {
         stop("'confidence' must be a numeric value in [0, 1]")  }

  if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE] 
      if(!is.null(id)) {
          if (min(id %in% names(dataset))!=1) stop("'id' does not exist in 'dataset'!")
          if (min(id %in% names(dataset))==1) id <- dataset[, id, with=FALSE]  }

      if (!is.null(period)) {
           if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
           if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE] }
      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }
      if(!is.null(H2)) {
          if (min(H2 %in% names(dataset))!=1) stop("'H2' does not exist in 'dataset'!")
          if (min(H2 %in% names(dataset))==1) H2 <- dataset[, H2, with=FALSE] }
      if(!is.null(PSU)) {
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) PSU <- dataset[, PSU, with=FALSE] }
      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }
      if(!is.null(Z)) {
          if (min(Z %in% names(dataset))!=1) stop("'Z' does not exist in 'dataset'!")
          if (min(Z %in% names(dataset))==1) Z <- dataset[, Z, with=FALSE] }
      if(!is.null(X)) {
          if (min(X %in% names(dataset))!=1) stop("'X' does not exist in 'dataset'!")
          if (min(X %in% names(dataset))==1) X <- dataset[, X, with=FALSE] }
      if(!is.null(g)) {
          if (min(g %in% names(dataset))!=1) stop("'g' does not exist in 'dataset'!")
          if (min(g %in% names(dataset))==1) g <- dataset[, g, with=FALSE] }
      if(!is.null(q)) {
          if (min(q %in% names(dataset))!=1) {
              if (length(q)!=nrow(dataset))  stop("'q' does not exist in 'dataset'!") }
          if (min(q %in% names(dataset))==1) q <- dataset[, q, with=FALSE] } 
      if (!is.null(Dom)) {
          if (min(Dom %in% names(dataset))!=1) stop("'Dom' does not exist in 'data'!")
          if (min(Dom %in% names(dataset))==1) Dom <- dataset[, Dom, with=FALSE]   }
    }

  # Y
  Y <- data.table(Y, check.names=TRUE)
  n <- nrow(Y)
  m <- ncol(Y)
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numeric values")
  if (any(is.na(Y))) stop("'Y' has unknown values")
  if (is.null(names(Y))) stop("'Y' must be colnames")
  
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

  # PSU
  PSU <- data.table(PSU)
  if (any(is.na(PSU))) stop("'PSU' has unknown values")
  if (nrow(PSU) != n) stop("'PSU' length must be equal with 'Y' row count")
  if (ncol(PSU) != 1) stop("'PSU' has more than 1 column")
  PSU[, (names(PSU)):=lapply(.SD, as.character)]
  
  # id
  if (is.null(id)) id <- PSU
  id <- data.table(id)
  if (any(is.na(id))) stop("'id' has unknown values")
  if (nrow(id) != n) stop("'id' length must be equal with 'Y' row count")
  if (ncol(id) != 1) stop("'id' must be 1 column data.frame, matrix, data.table")
  if (is.null(names(id))||(names(id)=="id")) setnames(id,names(id),"ID")

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

  # w_final 
  w_final <- data.frame(w_final)
  if (nrow(w_final) != n) stop("'w_final' must be equal with 'Y' row count")
  if (ncol(w_final) != 1) stop("'w_final' must be vector or 1 column data.frame, matrix, data.table")
  w_final <- w_final[,1]
  if (!is.numeric(w_final)) stop("'w_final' must be numerical")
  if (any(is.na(w_final))) stop("'w_final' has unknown values") 

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
             if (any(duplicated(N_h[, head(names(N_h),-1), with=FALSE]))) stop("Strata values for 'N_h' must be unique")
       } else { pH <- data.table(period, H)
                if (any(names(pH) != names(N_h)[c(1:(1+np))])) stop("Strata titles for 'period' with 'H' and 'N_h' is not equal")
                nperH <- names(period)
                if (pH[, class(get(nperH))]!=N_h[, class(get(nperH))]) 
                                                       stop("Period class for 'period' and 'N_h' is not equal ")
                if (any(is.na(merge(unique(pH), N_h, by=names(pH), all.x=TRUE)))) stop("'N_h' is not defined for all stratas and periods")
                if (any(duplicated(N_h[, head(names(N_h),-1), with=FALSE]))) stop("Strata values for 'N_h' must be unique in all periods")
               pH <- NULL
     }
    setkeyv(N_h, names(N_h)[c(1:(1+np))])
  } 

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
                if (pH2[, class(get(nperH))]!=N_h2[, class(get(nperH))]) 
                                                       stop("Period class for 'period' and 'N_h2' is not equal ")
                if (any(is.na(merge(unique(pH2), N_h2, by=names(pH2), all.x=TRUE)))) stop("'N_h2' is not defined for all stratas and periods")
                } 
    setkeyv(N_h2, names(N_h2)[c(1:(1+np))])
  } else stop ("N_h2 is not defined!")

  if (all(names(H)==names(H2))) {
      if (!is.null(N_h2))  setnames(N_h2, names(N_h2), paste0(names(N_h2),"2"))  
     setnames(H2, names(H2), paste0(names(H),"2"))  }

  # Dom
  namesDom <- NULL
  if (!is.null(Dom)) {
    Dom <- data.table(Dom)
    if (any(duplicated(names(Dom)))) 
           stop("'Dom' are duplicate column names: ", 
                 paste(names(Dom)[duplicated(names(Dom))], collapse = ","))
    if (nrow(Dom) != n) stop("'Dom' and 'Y' must be equal row count")
    if (any(is.na(Dom))) stop("'Dom' has unknown values")
    if (is.null(names(Dom))) stop("'Dom' must be colnames")
    Dom[, (names(Dom)):=lapply(.SD, as.character)]
    namesDom <- names(Dom)
  }
  
  # Z
  if (!is.null(Z)) {
    Z <- data.table(Z)
    if (nrow(Z) != n) stop("'Z' and 'Y' must be equal row count")
    if (ncol(Z) != m) stop("'Z' and 'Y' must be equal column count")
    if (!all(sapply(Z, is.numeric))) stop("'Z' must be numeric values")
    if (any(is.na(Z))) stop("'Z' has unknown values")
    if (is.null(names(Z))) stop("'Z' must be colnames")
  }
      
  # X
  if (!is.null(X)) {
    X <- data.table(X)
    if (any(is.na(X))) stop("'X' has unknown values")
    if (!all(sapply(X, is.numeric))) stop("'X' must be numeric values")
    if (nrow(X) != n) stop("'X' and 'Y' must be equal row count")
  }
      
  # g
  if (!is.null(X)) {
    if (is.null(class(g)) | all(class(g)=="function")) stop("'g' must be numerical")
    g <- data.frame(g)
    if (nrow(g) != nrow(X)) stop("'g' length must be equal with 'X' row count")
    if (ncol(g) != 1) stop("'g' must be 1 column data.frame, matrix, data.table")
    g <- g[,1]
    if (!is.numeric(g)) stop("'g' must be numerical")
    if (any(is.na(g))) stop("'g' has unknown values")
    if (any(g == 0)) stop("'g' value can not be 0")
   }
    
  # q
  if (!is.null(X)) {
    if (is.null(q))  q <- rep(1, nrow(X))
    if (is.null(class(q)) | all(class(q)=="function")) stop("'q' must be numerical")
    q <- data.frame(q)
    if (nrow(q) != nrow(X)) stop("'q' length must be equal with 'X' row count")
    if (ncol(q) != 1) stop("'q' must be 1 column data.frame, matrix, data.table")
    q <- q[,1]
    if (!is.numeric(q)) stop("'q' must be numerical")
    if (any(is.na(q))) stop("'q' has unknown values")
    if (any(is.infinite(q))) stop("'q' value can not be infinite")
  }


  ### Calculation
      
  # Domains
  if (!is.null(Dom)) Y1 <- domain(Y, Dom) else Y1 <- Y

  n_nonzero <- copy(Y1)
  if (!is.null(period)){ n_nonzero <- data.table(period, n_nonzero) 
                         n_nonzero <- n_nonzero[, lapply(.SD, function(x) 
                                                         sum(as.integer(abs(x)> .Machine$double.eps))),
                                                         keyby=names(period),
                                                         .SDcols = names(Y1)]
                  } else n_nonzero <- n_nonzero[, lapply(.SD, function(x) 
                                                         sum(as.integer(abs(x)> .Machine$double.eps))),
                                                         .SDcols = names(Y1)]

  respondent_count <- pop_size <- NULL
  nhs <- data.table(respondent_count=1, pop_size=w_final)
  if (!is.null(period)) nhs <- data.table(period, nhs)
  if (!is.null(Dom)) nhs <- data.table(Dom, nhs)
  if (!is.null(c(Dom, period))) {nhs <- nhs[, lapply(.SD, sum, na.rm=TRUE),
                                                       keyby=eval(names(nhs)[0:1-ncol(nhs)]),
                                                      .SDcols=c("respondent_count", "pop_size")]
                          } else nhs <- nhs[, lapply(.SD, sum, na.rm=TRUE),
                                                     .SDcols=c("respondent_count", "pop_size")]


  # Design weights
  if (!is.null(X)) w_design <- w_final / g else w_design <- w_final
      
  # Ratio of two totals
  linratio_outp <- per <- variableZ <- estim <- deff_sam <- NULL
  deff_est <- deff <- var_est2 <- se <- rse <- cv <- NULL
  absolute_margin_of_error <- relative_margin_of_error <- NULL
  sar_nr <- CI_lower <- CI_upper <- variable <- n_eff <- NULL

  idper <- id
  if (!is.null(period)) idper <- data.table(idper, period)

  Z1 <- NULL
  if (!is.null(Z)) {
    if (!is.null(Dom)) Z1 <- domain(Z, Dom) else Z1 <- Z
    if (is.null(period)) {
          Y2 <- lin.ratio(Y1, Z1, w_final, Dom=NULL)
          Y2a <- lin.ratio(Y1, Z1, w_design, Dom=NULL)
        } else {
            periodap <- do.call("paste", c(as.list(period), sep="_"))
            lin1 <- lapply(split(Y1[, .I], periodap), function(i)
                            data.table(sar_nr=i, 
                                   lin.ratio(Y1[i], Z1[i], w_final[i],
                                     Dom=NULL, percentratio=percentratio)))
            Y2 <- rbindlist(lin1)
            setkeyv(Y2, "sar_nr")
            lin2 <- lapply(split(Y1[, .I], periodap), function(i)
                            data.table(sar_nr=i, 
                                       lin.ratio(Y1[i], Z1[i], w_design[i],
                                       Dom=NULL, percentratio=percentratio)))
            Y2a <- rbindlist(lin2)
            setkeyv(Y2a, "sar_nr")
            Y2[, sar_nr:=NULL]
            Y2a[, sar_nr:=NULL]
        }
    if (any(is.na(Y2))) print("Results are calculated, but there are cases where Z = 0")
    if (outp_lin) linratio_outp <- data.table(idper, PSU, Y2) 
  } else {
          Y2 <- Y1
          Y2a <- Y1
         }
  Y <- Z <- NULL

  # Calibration
  res_outp <- NULL
  if (!is.null(X)) {
        ind_gr <- data.table(nsk=rep(1, nrow(X)))
        if (!is.null(period)) ind_gr <- data.table(ind_gr, period)
        ind_gr <- do.call("paste", c(as.list(ind_gr), sep="_"))

        lin1 <- lapply(split(Y2[,.I], ind_gr), function(i) 
                        data.table(sar_nr=i, residual_est(Y=Y2[i],
                                    X=X[i], weight=w_design[i], q=q[i])))
        Y3 <- rbindlist(lin1)
        setkeyv(Y3, "sar_nr")
        Y3[, sar_nr:=NULL] 
      if (outp_res) res_outp <- data.table(idper, PSU, Y3)
  } else Y3 <- Y2
  Y2 <- NULL

  var_est <- variance_othstr(Y=Y3, H=H, H2=H2,  
                             w_final=w_final, N_h=N_h,
                             N_h2=N_h2, period=period, dataset=NULL)
  s2g <- var_est$s2g
  var_est <- var_est$var_est
  var_est <- transpos(var_est, is.null(period), "var_est", names(period))
  all_result <- var_est

  n_nonzero <- transpos(n_nonzero, is.null(period), "n_nonzero", names(period))
  all_result <- merge(all_result, n_nonzero, all=TRUE)

  # Variance of HT estimator under current design
  var_cur_HT <- variance_othstr(Y=Y2a, H=H, H2=H2, 
                                w_final=w_design, N_h=N_h,
                                N_h2=N_h2, period=period, dataset=NULL)
  var_cur_HT <- var_cur_HT$var_est
  var_cur_HT <- transpos(var_cur_HT, is.null(period), "var_cur_HT", names(period))
  all_result <- merge(all_result, var_cur_HT)
  n_nonzero <- var_est <- var_cur_HT <- NULL

  # Variance of HT estimator under SRS
  if (is.null(period)) {
           var_srs_HT <- var_srs(Y2a, w = w_design)
       } else {
           period_agg <- unique(period)
           lin1 <- lapply(1:nrow(period_agg), function(i) {
                          per <- period_agg[i,][rep(1, nrow(Y2a)),]
                          ind <- (rowSums(per == period) == ncol(period))
                          data.table(period_agg[i,], 
                                     var_srs(Y2a[ind], w = w_design[ind]))
                        })
           var_srs_HT <- rbindlist(lin1)
      }
  var_srs_HT <- transpos(var_srs_HT, is.null(period), "var_srs_HT", names(period))
  all_result <- merge(all_result, var_srs_HT)


  # Variance of calibrated estimator under SRS
  if (is.null(period)) {
           var_srs_ca <- var_srs(Y3, w = w_final)
      } else {
           period_agg <- unique(period)
           lin1 <- lapply(1:nrow(period_agg), function(i) {
                          per <- period_agg[i,][rep(1, nrow(Y2a)),]
                          ind <- (rowSums(per == period) == ncol(period))
                          data.table(period_agg[i,], 
                                     var_srs(Y3[ind], w = w_final[ind]))
                        })
           var_srs_ca <- rbindlist(lin1)
        }
  Y3 <- Y2a <- NULL
  var_srs_ca <- transpos(var_srs_ca, is.null(period), "var_srs_ca", names(period))
  all_result <- merge(all_result, var_srs_ca)


  # Total estimation
  Y_nov <- Z_nov <- .SD <- NULL

  hY <- data.table(Y1*w_final)
  if (is.null(period)) { Y_nov <- hY[, lapply(.SD, sum, na.rm=TRUE), .SDcols = names(Y1)]
                } else { hY <- data.table(period, hY)
                         Y_nov <- hY[, lapply(.SD, sum, na.rm=TRUE), keyby=names(period), .SDcols = names(Y1)]
                       }
  Y_nov <- transpos(Y_nov, is.null(period), "Y_nov", names(period))

  all_result <- merge(all_result, Y_nov)
  
  if (!is.null(Z1)) {
         YZnames <- data.table(variable=names(Y1), variableDZ=names(Z1))
         setkeyv(YZnames, "variable")
         setkeyv(all_result, "variable")
         all_result <- merge(all_result, YZnames)
         
         hZ <- data.table(Z1*w_final)
         if (is.null(period)) { Z_nov <- hZ[, lapply(.SD, sum, na.rm=TRUE), .SDcols = names(Z1)]
                       } else { hZ <- data.table(period, hZ)
                                Z_nov <- hZ[, lapply(.SD, sum, na.rm=TRUE), keyby=names(period), .SDcols = names(Z1)]
                              }
         Z_nov <- transpos(Z_nov, is.null(period), "Z_nov", names(period), "variableDZ")
         setkeyv(all_result, "variableDZ")
         all_result <- merge(all_result, Z_nov)                                            
      }

  vars <- data.table(variable=names(Y1), nr_names=1:ncol(Y1))
  setkey(vars, "variable")
  setkey(all_result, "variable")
  all_result <- merge(vars, all_result)
                        
  vars <- idper <- Y1 <- Z1 <- Y_nov <- NULL
  Z_nov <- hY <- hZ <- YZnames <- dati <- NULL                            

  
  all_result[, estim:=Y_nov]   
  if (!is.null(all_result$Z_nov)) all_result[, estim:=Y_nov/Z_nov]

  if (nrow(all_result[var_est < 0])>0) print("Estimation of variance are negative!")
 
  # Design effect of sample design
  all_result[, deff_sam:=var_cur_HT / var_srs_HT]
  
  # Design effect of estimator
  all_result[, deff_est:= var_est / var_cur_HT]
  
  # Overall effect of sample design and estimator
  all_result[, deff:= deff_sam * deff_est]

  all_result[, var_est2:=var_est]
  all_result[xor(is.na(var_est2), var_est2 < 0), var_est2:=NA]
  all_result[, se:=sqrt(var_est2)]
  all_result[(estim!=0) & !is.nan(estim), rse:= se/estim]
  all_result[estim==0 | is.nan(estim), rse:=NA]
  all_result[, cv:= rse*100]

  tsad <- qnorm(0.5*(1+confidence))
  all_result[, absolute_margin_of_error:= tsad*se]
  all_result[, relative_margin_of_error:= tsad*cv]
  all_result[, CI_lower:= estim - tsad*se]
  all_result[, CI_upper:= estim + tsad*se]

  setnames(all_result, c("variable", "var_est"), c("variableD", "var"))
  if (!is.null(all_result$Z_nov)) {
                         nosrZ <- all_result$variableDZ
                         nosrZ <- nosrZ[!duplicated(nosrZ)]
                         nosrZ1 <- data.table(variableZ=t(data.frame(strsplit(nosrZ, "__")))[,c(1)])
                         nosrZ <- data.table(variableDZ=nosrZ, nosrZ1)
                         setkeyv(nosrZ, "variableDZ")
                         setkeyv(all_result, "variableDZ")
                         all_result <- merge(all_result, nosrZ)
                         nosrZ <- nosrZ1 <- NULL
                       }

  nosr <- data.table(variableD=all_result$variableD, t(data.frame(strsplit(all_result$variableD, "__"))))
  nosr <- nosr[!duplicated(nosr)]
  nosr <- nosr[, lapply(nosr, as.character)]
  setnames(nosr, names(nosr)[2], "variable")

  namesDom1 <- namesDom
  if (!is.null(Dom)) {
       setnames(nosr, names(nosr)[3:ncol(nosr)], paste0(namesDom, "_new"))
       nhs[, (paste0(namesDom, "_new")):=lapply(namesDom, function(x) make.names(paste0(x,".", get(x))))]
       namesDom1 <- paste0(namesDom, "_new")
    }

  setkeyv(nosr, "variableD")
  setkeyv(all_result, "variableD")
  all_result <- merge(nosr, all_result)
  namesDom <- nosr <- NULL
  
  if (!is.null(all_result$Z_nov)) {
       all_result[, variable:=paste("R", get("variable"), get("variableZ"), sep="__")] }
  setkeyv(all_result, c(namesDom1, names(period)))
  setkeyv(nhs, c(namesDom1, names(period)))

  if (!is.null(c(Dom, period))) { all_result <- merge(all_result, nhs, all=TRUE)
                         } else { all_result[, respondent_count:=nhs$respondent_count]
                                  all_result[, pop_size:=nhs$pop_size]} 

  variab <- c("respondent_count", "n_nonzero", "pop_size", "estim", "var", "se", 
              "rse", "cv", "absolute_margin_of_error", "relative_margin_of_error",
              "CI_lower", "CI_upper", "var_srs_HT",  "var_cur_HT", 
              "var_srs_ca", "deff_sam", "deff_est", "deff")

  setkeyv(all_result, c("nr_names", names(Dom), names(period)))
  all_result <- all_result[, c("variable", names(Dom), names(period), variab), with=FALSE]

  list(lin_out = linratio_outp,
       res_out = res_outp,
       s2g = s2g,
       all_result = all_result)
}
