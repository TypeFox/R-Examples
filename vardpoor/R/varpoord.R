varpoord <- function(Y, w_final,
                     age=NULL,
                     pl085=NULL,
                     month_at_work=NULL,
                     Y_den=NULL,
                     Y_thres = NULL,
                     wght_thres = NULL,                    
                     ID_household,
                     id = NULL, 
                     H, PSU, N_h,
                     fh_zero = FALSE,
                     PSU_level=TRUE,
                     sort = NULL,
                     Dom = NULL,
                     period = NULL,
                     gender = NULL,
                     dataset = NULL,
                     X = NULL,
                     periodX = NULL,
                     X_ID_household = NULL,
                     ind_gr = NULL,
                     g = NULL,
                     q = NULL,
                     datasetX = NULL,
                     percentage=60,
                     order_quant=50,
                     alpha = 20,
                     confidence = .95,
                     outp_lin = FALSE,
                     outp_res = FALSE,
                     type="linrmpg") {

  ### Checking

  if (length(fh_zero) != 1 | !any(is.logical(fh_zero))) stop("'fh_zero' must be the logical value")
  if (length(PSU_level) != 1 | !any(is.logical(PSU_level))) stop("'PSU_level' must be the logical value")
  if (length(outp_lin) != 1 | !any(is.logical(outp_lin))) stop("'outp_lin' must be the logical value")
  if (length(outp_res) != 1 | !any(is.logical(outp_res))) stop("'outp_res' must be the logical value")

  all_choices <- c("linarpr","linarpt","lingpg","linpoormed",
                   "linrmpg","lingini","lingini2", "linqsr", "linrmir", "linarr")
  choices <- c("all_choices", all_choices)
  type <- tolower(type)
  type <- match.arg(type, choices, length(type)>1) 
  if (any(type == "all_choices")) type <- all_choices 

  # check 'p'
  p <- percentage
   if(length(p) != 1 |  any(!is.numeric(p) | p < 0 | p > 100)) {
          stop("'percentage' must be a numeric value in [0, 100]")  }

  # check 'order_quant'
  oq <- order_quant
   if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]")  }

  if(length(alpha) != 1 | any(!is.numeric(alpha) | alpha < 0 | alpha > 100)) {
         stop("'alpha' must be a numeric value in [0,100]")  }

  if(length(confidence) != 1 | any(!is.numeric(confidence) | confidence < 0 | confidence > 1)) {
         stop("'confidence' must be a numeric value in [0, 1]")  }

  if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE]
      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }
      if(!is.null(age)) {
          if (min(age %in% names(dataset))!=1) stop("'age' does not exist in 'dataset'!")
          if (min(age %in% names(dataset))==1) age <- dataset[, age, with=FALSE] }
      if(!is.null(pl085)) {
          if (min(pl085 %in% names(dataset))!=1) stop("'pl085' does not exist in 'dataset'!")
          if (min(pl085 %in% names(dataset))==1) pl085 <- dataset[, pl085, with=FALSE] }
      if(!is.null(month_at_work)) {
          if (min(month_at_work %in% names(dataset))!=1) stop("'month_at_work' does not exist in 'dataset'!")
          if (min(month_at_work %in% names(dataset))==1) month_at_work <- dataset[, month_at_work, with=FALSE] }
      if(!is.null(Y_den)) {
          if (min(Y_den %in% names(dataset))!=1) stop("'Y_den' does not exist in 'dataset'!")
          if (min(Y_den %in% names(dataset))==1) Y_den <- dataset[, Y_den, with=FALSE] }
      if(!is.null(Y_thres)) {
          if (min(Y_thres %in% names(dataset))!=1) stop("'Y_thres' does not exist in 'dataset'!")
          if (min(Y_thres %in% names(dataset))==1) Y_thres <- dataset[, Y_thres, with=FALSE] }    
      if(!is.null(wght_thres)) {
          if (min(wght_thres %in% names(dataset))!=1) stop("'wght_thres' does not exist in 'dataset'!")
          if (min(wght_thres %in% names(dataset))==1) wght_thres <- dataset[, wght_thres, with=FALSE] }
      if(!is.null(id)) {
          if (min(id %in% names(dataset))!=1) stop("'id' does not exist in 'dataset'!")
          if (min(id %in% names(dataset))==1) id <- dataset[, id, with=FALSE]  }
      if(!is.null(ID_household)) {
          if (min(ID_household %in% names(dataset))!=1) stop("'ID_household' does not exist in 'dataset'!")
          if (min(ID_household %in% names(dataset))==1) ID_household <- dataset[, ID_household, with=FALSE] }
      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }
      if(!is.null(PSU)) {
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) PSU <- dataset[, PSU, with=FALSE]  }
      if(!is.null(gender)) {
          if (min(gender %in% names(dataset))!=1) stop("'gender' does not exist in 'dataset'!")
          if (min(gender %in% names(dataset))==1) gender <- dataset[, gender, with=FALSE] }
      if(!is.null(sort)) {
          if (min(sort %in% names(dataset))!=1) stop("'sort' does not exist in 'dataset'!")
          if (min(sort %in% names(dataset))==1) sort <- dataset[, sort, with=FALSE] }
      if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE] }
      if (!is.null(Dom)) {
          if (min(Dom %in% names(dataset))!=1) stop("'Dom' does not exist in 'dataset'!")
          if (min(Dom %in% names(dataset))==1) Dom <- dataset[, Dom, with=FALSE] }
    }

  if(!is.null(datasetX)) {
      dataset <- data.table(datasetX)
       if (!is.null(periodX)) {
            if (min(periodX %in% names(datasetX))!=1) stop("'periodX' does not exist in 'datasetX'!")
            if (min(periodX %in% names(datasetX))==1) periodX <- datasetX[, periodX, with=FALSE] }     
      if(!is.null(X_ID_household)) {
          if (min(X_ID_household %in% names(datasetX))!=1) stop("'X_ID_household' does not exist in 'datasetX'!")
          if (min(X_ID_household %in% names(datasetX))==1) X_ID_household <- datasetX[, X_ID_household, with=FALSE] }
      if(!is.null(X)) {
          if (min(X %in% names(datasetX))!=1) stop("'X' does not exist in 'datasetX'!")
          if (min(X %in% names(datasetX))==1) X <- datasetX[, X, with=FALSE] }
      if(!is.null(ind_gr)) {
          if (min(ind_gr %in% names(datasetX))!=1) stop("'ind_gr' does not exist in 'datasetX'!")
          if (min(ind_gr %in% names(datasetX))==1) ind_gr <- datasetX[, ind_gr, with=FALSE] }
      if(!is.null(g)) {
          if (min(g %in% names(datasetX))!=1) stop("'g' does not exist in 'datasetX'!")
          if (min(g %in% names(datasetX))==1) g <- datasetX[, g, with=FALSE] }
      if(!is.null(q)) {
          if (min(q %in% names(datasetX))!=1) {
               if (length(q) != nrow(datasetX))  stop("'q' does not exist in 'datasetX'!") }
          if (min(q %in% names(datasetX))==1) q <- datasetX[, q, with=FALSE]  }
    }
  N <- dataset <- datasetX <- NULL

  # Y
  Y <- data.frame(Y)
  n <- nrow(Y)
  if (ncol(Y) != 1) stop("'Y' must have vector or 1 column data.frame, matrix, data.table")
  Y <- Y[,1]
  if (!is.numeric(Y)) stop("'Y' must be numerical")
  if (any(is.na(Y))) stop("'Y' has unknown values")
  
   if (!is.null(Y_den)) {
          Y_den <- data.frame(Y_den)
          if (ncol(Y_den) != 1) stop("'Y_den' must be vector or 1 column data.frame, matrix, data.table")
          if (nrow(Y_den) != n) stop("'Y_den' must be the same length as 'Y'")
          Y_den <- Y_den[,1]
          if(!is.numeric(Y_den)) stop("'Y_den' must be numerical")
          if (any(is.na(Y_den))) stop("'Y_den' has unknown values")
  }
 
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

  # id
  if (is.null(id)) id <- 1:n
  id <- data.table(id)
  if (any(is.na(id))) stop("'id' has unknown values")
  if (ncol(id) != 1) stop("'id' must be 1 column data.frame, matrix, data.table")
  if (nrow(id) != n) stop("'id' must be the same length as 'Y'")
  if (is.null(names(id))||(names(id)=="id")) setnames(id,names(id),"ID")
  if (is.null(period)){ if (any(duplicated(id))) stop("'id' are duplicate values") 
                       } else {dd <- data.table(period, id)
                                  if (any(duplicated(dd, by=names(dd)))) stop("'id' by period are duplicate values")
                                  dd <- NULL}

  # age
  if (!is.null(age)) {
       age <- data.frame(age)
       if (nrow(age) != n) stop("'age' must be the same length as 'Y'")
       if (ncol(age) != 1) stop("'age' must be vector or 1 column data.frame, matrix, data.table")
      age <- age[, 1]
      if (!is.numeric(age)) stop("'age' must be numerical")
      if (any(is.na(age))) stop("'age' has unknown values")
   }

   # pl085
   if (!is.null(pl085)) {
       pl085 <- data.frame(pl085)
       if (nrow(pl085) != n) stop("'pl085' must be the same length as 'Y'")
       if (ncol(pl085) != 1) stop("'pl085' must be vector or 1 column data.frame, matrix, data.table")
       pl085 <- pl085[, 1]
       if (!is.numeric(pl085)) stop("'pl085' must be numerical")
       if (any(is.na(pl085))) stop("'pl085' has unknown values")
   }

   # month_at_work
   if (!is.null(month_at_work)) {
        month_at_work <- data.frame(month_at_work)
        if (nrow(month_at_work) != n) stop("'month_at_work' must be the same length as 'Y'")
        if (ncol(month_at_work) != 1) stop("'month_at_work' must be vector or 1 column data.frame, matrix, data.table")
        month_at_work <- month_at_work[, 1]
        if (!is.numeric(pl085)) stop("'month_at_work' must be numerical")
        if (any(is.na(pl085))) stop("'month_at_work' has unknown values")
  }

  # ID_household
  if (is.null(ID_household)) stop("'ID_household' must be defined")
  ID_household <- data.table(ID_household)
  if (ncol(ID_household) != 1) stop("'ID_household' must be 1 column data.frame, matrix, data.table")
  if (nrow(ID_household) != n) stop("'ID_household' must be the same length as 'Y'")
  if (is.null(names(ID_household))) setnames(ID_household,names(ID_household),"ID_household")
  if (names(id)==names(ID_household)) setnames(id,names(id),paste(names(id),"_id",sep=""))

  # w_final 
  w_final <- data.frame(w_final)
  if (nrow(w_final) != n) stop("'w_final' must have the same length as 'Y'")
  if (ncol(w_final) != 1) stop("'w_final' must have vector or 1 column data.frame, matrix, data.table")
  w_final <- w_final[,1]
  if (!is.numeric(w_final)) stop("'w_final' must be numerical")
  if (any(is.na(w_final))) stop("'w_final' has unknown values") 
  
  # Y_thres
  if (!is.null(Y_thres)) {
       Y_thres <- data.frame(Y_thres)
       if (nrow(Y_thres) != n) stop("'Y_thres' must have the same length as 'Y'")
       if (ncol(Y_thres) != 1) stop("'Y_thres' must have vector or 1 column data.frame, matrix, data.table")
       Y_thres <- Y_thres[,1]
       if (!is.numeric(Y_thres)) stop("'Y_thres' must be numerical")
       if (any(is.na(Y_thres))) stop("'Y_thres' has unknown values") 
     } else Y_thres <- Y

  # wght_thres
  if (is.null(wght_thres)) wght_thres <- w_final
  wght_thres <- data.frame(wght_thres)
  if (nrow(wght_thres) != n) stop("'wght_thres' must have the same length as 'Y'")
  if (ncol(wght_thres) != 1) stop("'wght_thres' must have vector or 1 column data.frame, matrix, data.table")
  wght_thres <- wght_thres[,1]
  if (!is.numeric(wght_thres)) stop("'wght_thres' must be a numeric vector")
 
  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' must have the same length as 'Y'")
  if (ncol(H) != 1) stop("'H' must have 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colnames")
  H[, (names(H)):=lapply(.SD, as.character)]

  # PSU
  PSU <- data.table(PSU)
  if (nrow(PSU) != n) stop("'PSU' must have the same length as 'Y'")
  if (ncol(PSU) != 1) stop("'PSU' must have vector or 1 column data.frame, matrix, data.table")
  if (any(is.na(PSU))) stop("'PSU' has unknown values")
  PSU[, (names(PSU)):=lapply(.SD, as.character)]

  # gender
  if (!is.null(gender)) {
      gender <- data.frame(gender)
      if (nrow(gender) != n) stop("'gender' must be the same length as 'Y'")
      if (ncol(gender) != 1) stop("'gender' must be vector or 1 column data.frame, matrix, data.table")
      gender <- gender[,1]
      if (!is.numeric(gender)) stop("'gender' must be numerical")
      if (length(unique(gender)) != 2) stop("'gender' must be exactly two values")
      if (!all.equal(unique(gender),c(1, 2))) stop("'gender' must be value 1 for male, 2 for females")
   }

  # N_h
  if (!is.null(N_h)) {
      N_h <- data.table(N_h)
      if (ncol(N_h) != np+2) stop(paste0("'N_h' should be ", np+2, " columns"))
      if (!is.numeric(N_h[[ncol(N_h)]])) stop("The last column of 'N_h' should be numerical")
      if (any(is.na(N_h))) stop("'N_h' has unknown values") 
      if (is.null(names(N_h))) stop("'N_h' must be colnames")
      if (all(names(H) %in% names(N_h))) {N_h[, (names(H)):=lapply(.SD, as.character), .SDcols=names(H)]
             } else stop("All strata titles of 'H' have not in 'N_h'")
      if (is.null(period)) {
             if (names(H) != names(N_h)[1]) stop("Strata titles for 'H' and 'N_h' is not equal")
             if (any(is.na(merge(unique(H), N_h, by=names(H), all.x = T)))) stop("'N_h' is not defined for all stratas")
             if (any(duplicated(N_h[, head(names(N_h),-1), with=F]))) stop("Strata values for 'N_h' must be unique")
       } else { pH <- data.table(period, H)
                if (any(names(pH) != names(N_h)[c(1:(1+np))])) stop("Strata titles for 'period' with 'H' and 'N_h' is not equal")
                nperH <- names(period)
                if (pH[, class(get(nperH))]!=N_h[, class(get(nperH))]) 
                                                       stop("Period class for 'period' and 'N_h' is not equal ")
                if (any(is.na(merge(unique(pH), N_h, by=names(pH), all.x = TRUE)))) stop("'N_h' is not defined for all stratas and periods")
                if (any(duplicated(N_h[, head(names(N_h),-1), with=FALSE]))) stop("Strata values for 'N_h' must be unique in all periods")
                pH <- NULL
     }
    setkeyv(N_h, names(N_h)[c(1:(1+np))])
  }

  # sort
  if (!is.null(sort) && !is.vector(sort) && !is.ordered(sort)) {
        stop("'sort' must be a vector or ordered factor") }
  if (!is.null(sort) && length(sort) != n) stop("'sort' must have the same length as 'Y'")     

  # Dom
  if (!is.null(Dom)) {
    Dom <- data.table(Dom)
    if (any(duplicated(names(Dom)))) 
           stop("'Dom' are duplicate column names: ", 
                 paste(names(Dom)[duplicated(names(Dom))], collapse = ","))
    if (nrow(Dom) != n) stop("'Dom' and 'Y' have different row count")
    if (any(is.na(Dom))) stop("'Dom' has unknown values")
    if (is.null(names(Dom))) stop("'Dom' must be colnames")
    Dom <- Dom[, lapply(.SD, as.character), .SDcols = names(Dom)]
  }

  # X
  if (!is.null(X)) {
    X <- data.table(X, check.names=TRUE)
    if (!all(sapply(X, is.numeric))) stop("'X' must be numeric values")
  }

  # periodX
  if (!is.null(X)) {
     if(!is.null(periodX)) {
        periodX <- data.table(periodX)
        periodX[, (names(periodX)):=lapply(.SD, as.character)]
        periX <- data.table(unique(periodX))
        setkeyv(periX, names(periX))
        peri <- data.table(unique(period))
        setkeyv(peri, names(peri))
        if (any(duplicated(names(periodX)))) 
                    stop("'periodX' are duplicate column names: ", 
                         paste(names(periodX)[duplicated(names(periodX))], collapse = ","))
        if (nrow(periodX) != nrow(X)) stop("'periodX' length must be equal with 'X' row count")
        if (ncol(periodX) != ncol(period)) stop("'periodX' length must be equal with 'period' column count")
        if (names(periodX) != names(period)) stop("'periodX' must be equal with 'period' names")
        if (any(is.na(periodX))) stop("'periodX' has unknown values")
        if (any(peri != periX)) stop("'unique(period)' and 'unique(periodX)' records have different")
        if (peri[, class(get(names(peri)))]!=periX[, class(get(names(periX)))])  stop("Class for 'periodX' and class for 'period' must be equal")
      } else if (!is.null(period)) stop("'periodX' must be defined")
   } 

 # X_ID_household
  if (!is.null(X)) {
    X_ID_household <- data.table(X_ID_household)
    X_ID_household[, (names(X_ID_household)):=lapply(.SD, as.character)]
    if (nrow(X) != nrow(X_ID_household)) stop("'X' and 'X_ID_household' have different row count")
    if (ncol(X_ID_household) != 1) stop("'X_ID_household' must be 1 column data.frame, matrix, data.table")
    if (any(is.na(X_ID_household))) stop("'X_ID_household' has unknown values")

    IDh <- data.table(unique(ID_household))
    if (!is.null(period)) { X_ID_household <- data.table(periodX, X_ID_household)
                            IDh <- data.table(period, ID_household)
                            IDh <- IDh[, .N, by=names(IDh)][, N:=NULL] }
    if (nrow(X_ID_household[,.N,by=names(X_ID_household)][N>1])>0) stop("'X_ID_household' have duplicates")
    setkeyv(X_ID_household, names(X_ID_household))
    setkeyv(IDh, names(IDh))

    nperIDh <- names(IDh)
    if (nperIDh != names(X_ID_household)) stop("'X_ID_household' and 'ID_household' must be equal names")
    if (IDh[, class(get(nperIDh))]!=X_ID_household[, class(get(nperIDh))])  stop("Class for 'X_ID_household' and class for 'ID_household' must be equal ")

    if (!is.null(period)) {
        if (nrow(IDh) != nrow(X_ID_household)) stop("'periodX' with 'X_ID_household' and 'unique(period, ID_household)' have different row count")
        if (any(IDh != X_ID_household)) stop("'periodX' with 'X_ID_household' and 'unique(period, ID_household)' records have different")
      } else {
        if (nrow(IDh) != nrow(X_ID_household)) stop("'X_ID_household' and 'unique(ID_household)' have different row count")
        if (any(IDh != X_ID_household)) stop("'X_ID_household' and 'unique(ID_household)' records have different")
    }}

  # ind_gr
  if (!is.null(X)) {
     if(is.null(ind_gr)) ind_gr <- rep.int(1, nrow(X)) 
     ind_gr <- data.table(ind_gr, check.names=TRUE)
     if (nrow(ind_gr) != nrow(X)) stop("'ind_gr' length must be equal with 'X' row count")
     if (ncol(ind_gr) != 1) stop("'ind_gr' must be 1 column data.frame, matrix, data.table")
     if (any(is.na(ind_gr))) stop("'ind_gr' has unknown values")
   }

  # X
  if (!is.null(X)) {
       X1 <- data.table(X, check.names=TRUE)
       nX1 <- names(X1)
       ind_gr1 <- copy(ind_gr) 
       if (!is.null(periodX)) ind_gr1 <- data.table(periodX, ind_gr1, check.names=TRUE)
       X2 <- data.table(ind_gr1, X1)
       X1 <- X2[, .N, keyby=names(ind_gr1)][[ncol(ind_gr1)+1]]
       X2 <- X2[,lapply(.SD, function(x) sum(!is.na(x))), keyby=names(ind_gr1), .SDcols=nX1]
       X2 <- X2[, !(names(X2) %in% names(ind_gr)), with=FALSE]
       if (!all(X2==0 | X1==X2)) stop("X has unknown values")
       ind_gr1 <- nX1 <- X1 <- X2 <- NULL
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

  # Design weights
  if (!is.null(X)) {
             idh <- data.frame(ID_household)
             if (!is.null(period)) idh <- data.table(period, idh)
             idhx <- data.table(X_ID_household, g)
             setnames(idhx, names(idhx)[c(1:(ncol(idhx)-1))], names(idh))
             idg <- data.table(merge(idh, idhx, by=names(idh), sort=FALSE))
             w_design <- w_final / idg[[ncol(idg)]]
             idg <- data.table(idg, w_design=w_design)
             idh <- idg[, .N, keyby=c(names(idh), "w_design")]
             if (nrow(X) != nrow(idh))  stop("Aggregated 'w_design' length must the same as matrix 'X'")
             idg <- idhx <- idh <- NULL
      } else w_design <- w_final

  ### Calculation
  sar_nr <- respondent_count <- pop_size <- n_nonzero <- NULL
  nhs <- data.table(respondent_count=1, pop_size=w_final, 
                               n_nonzero=as.integer(abs(Y)> .Machine$double.eps))
  if (!is.null(period)) nhs <- data.table(period, nhs)
  if (!is.null(Dom)) nhs <- data.table(Dom, nhs)
  if (!is.null(c(Dom, period))) {nhs <- nhs[, lapply(.SD, sum, na.rm=TRUE),
                                                       keyby=eval(names(nhs)[0:2-ncol(nhs)]),
                                                      .SDcols=c("respondent_count", "pop_size", "n_nonzero")]
                          } else nhs <- nhs[, lapply(.SD, sum, na.rm=TRUE),
                                                     .SDcols=c("respondent_count", "pop_size", "n_nonzero")]

  estim <- c()
  aH <- names(H)
  idper <- copy(id)
  Y1sort <- Y1asort <- NULL
  aPSU <- names(PSU)
  if (!is.null(period)) idper <- data.table(idper, period)
  Y1 <- data.table(idper, ID_household, H, PSU, w_final, check.names=TRUE)
  Y1a <- data.table(idper, ID_household, H, PSU, w_design, check.names=TRUE)
  Y1[, Y1sort:=.I]
  Y1a[, Y1asort:=.I]
  setkeyv(Y1, names(idper))
  setkeyv(Y1a, names(idper))
  value <- NULL

  if ("linarpt" %in% type) {
       varpt <- linarpt(Y=Y, id=id, weight=w_final,
                        sort=sort, Dom=Dom, period=period,
                        dataset=NULL, percentage=percentage,
                        order_quant=order_quant, var_name="lin_arpt")

       varpta <- linarpt(Y=Y, id=id, weight=w_design,
                         sort=sort, Dom=Dom, period=period,
                         dataset=NULL, percentage=percentage,
                         order_quant=order_quant, var_name="lin_arpt")

       Y1 <- merge(Y1, varpt$lin, all.x=TRUE)
       Y1a <- merge(Y1a, varpta$lin, all.x=TRUE)

       esti <- data.table("ARPT", varpt$value, NA)
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       varpt <- varpta <- esti <- NULL
     }
  if ("linarpr" %in% type) {
       varpr <- linarpr(Y=Y, id=id, weight=w_final,
                        Y_thres=Y_thres,
                        wght_thres=wght_thres, sort=sort, 
                        Dom=Dom, period=period, dataset=NULL, 
                        percentage=percentage,
                        order_quant=order_quant,
                        var_name="lin_arpr")
       varpra <- linarpr(Y=Y, id=id, weight=w_design,
                         Y_thres=Y_thres,
                         wght_thres=wght_thres, sort=sort,
                         Dom=Dom, period=period, dataset=NULL, 
                         percentage=percentage,
                         order_quant=order_quant,
                         var_name="lin_arpr")

       Y1 <- merge(Y1, varpr$lin, all.x=TRUE)
       Y1a <- merge(Y1a, varpra$lin, all.x=TRUE)

       esti <- data.table("ARPR", varpr$value, NA)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       varpr <- varpra <- esti <- NULL
     }
  if (("lingpg" %in% type)&&(!is.null(gender))) {
        vgpg <- lingpg(Y=Y, gender=gender, id=id,
                       weight=w_final, sort=sort,
                       Dom=Dom, period=period, dataset=NULL, 
                       var_name="lin_gpg")
        vgpga <- lingpg(Y=Y, gender=gender, id=id,
                        weight=w_design, sort=sort,
                        Dom=Dom, period=period, dataset=NULL, 
                        var_name="lin_gpg")

        Y1 <- merge(Y1, vgpg$lin, all.x=TRUE)
        Y1a <- merge(Y1a, vgpga$lin, all.x=TRUE)
     
        esti <- data.table("GPG", vgpg$value, NA)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vgpg <- vgpga <- esti <- NULL
     }
  if ("linpoormed" %in% type) {
        vporm <- linpoormed(Y=Y, id=id, weight=w_final,
                            sort=sort, Dom=Dom, period=period, 
                            dataset=NULL, percentage=percentage,
                            order_quant=order_quant, var_name="lin_poormed")
        vporma <- linpoormed(Y=Y, id=id, weight=w_design,
                             sort=sort, Dom=Dom, period=period, 
                             dataset=NULL, percentage=percentage,
                             order_quant=order_quant, var_name="lin_poormed")
        Y1 <- merge(Y1, vporm$lin, all.x=TRUE)
        Y1a <- merge(Y1a, vporma$lin, all.x=TRUE)

        esti <- data.table("POORMED", vporm$value, NA)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vporm <- vporma <- esti <- NULL
     }
  if ("linrmpg" %in% type) {
        vrmpg <- linrmpg(Y=Y, id=id, weight=w_final,
                         sort=sort, Dom=Dom, period=period,
                         dataset=NULL, percentage=percentage,
                         order_quant=order_quant, var_name="lin_rmpg")

        vrmpga <- linrmpg(Y=Y, id=id, weight=w_design,
                          sort=sort, Dom=Dom, period=period,
                          dataset=NULL, percentage=percentage,
                          order_quant=order_quant, var_name="lin_rmpg")

        Y1 <- merge(Y1, vrmpg$lin, all.x=TRUE)
        Y1a <- merge(Y1a, vrmpga$lin, all.x=TRUE)

        esti <- data.table("RMPG", vrmpg$value, NA)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu")) 
       estim <- rbind(estim, esti)
       vrmpg <- vrmpga <- esti <- NULL
      }
  if ("linqsr" %in% type) {
       vqsr <- linqsr(Y=Y, id=id, weight=w_final, 
                      sort=sort, Dom=Dom, period=period,
                      dataset=NULL, alpha=alpha, var_name="lin_qsr") 
       vqsra <- linqsr(Y=Y, id=id, weight=w_design,
                      sort=sort, Dom=Dom, period=period,
                      dataset=NULL, alpha=alpha, var_name="lin_qsr") 

       Y1 <- merge(Y1, vqsr$lin, all.x=TRUE)
       Y1a <- merge(Y1a, vqsra$lin, all.x=TRUE)

       esti <- data.table("QSR", vqsr$value)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       vqsr <- vqsra <- esti <- NULL
    }
  if ("lingini" %in% type) {
       vgini <- lingini(Y=Y, id=id, weight=w_final,
                        sort=sort, Dom=Dom, period=period,
                        dataset=NULL, var_name="lin_gini")
       vginia <- lingini(Y=Y, id=id, weight=w_design,
                         sort=sort, Dom=Dom, period=period,
                         dataset=NULL, var_name="lin_gini")

       Y1 <- merge(Y1, vgini$lin, all.x=TRUE)
       Y1a <- merge(Y1a, vginia$lin, all.x=TRUE)

       esti <- data.table("GINI", vgini$value)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       vgini <- vginia <- esti <- NULL
     }
  if ("lingini2" %in% type) {
       vgini2 <- lingini2(Y=Y, id=id, weight=w_final,
                          sort=sort, Dom=Dom, period=period,
                          dataset=NULL, var_name="lin_gini2")
       vgini2a <- lingini2(Y=Y, id=id, weight=w_design,
                          sort=sort, Dom=Dom, period=period,
                          dataset=NULL, var_name="lin_gini2")

       Y1 <- merge(Y1, vgini2$lin, all.x=TRUE)
       Y1a <- merge(Y1a, vgini2a$lin, all.x=TRUE)

       esti <- data.table("GINI2", vgini2$value)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       vgini2 <- vgini2a <- esti <- NULL
     }
  if (("linrmir" %in% type)&&(!is.null(age))) {
       vrmir <- linrmir(Y=Y, id=id, age=age, weight=w_final, 
                      sort=sort, Dom=Dom, period=period,
                      dataset=NULL,  order_quant=order_quant,
                      var_name="lin_rmir") 
       vrmira <- linrmir(Y=Y, id=id, age=age, weight=w_design,
                      sort=sort, Dom=Dom, period=period,
                      dataset=NULL, order_quant=order_quant,
                       var_name="lin_rmir") 

       Y1 <- merge(Y1, vrmir$lin, all.x=TRUE)
       Y1a <- merge(Y1a, vrmira$lin, all.x=TRUE)

       esti <- data.table("RMIR", vrmir$value, NA)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       vrmir <- vrmira <- esti <- NULL
    }
  if (("linarr" %in% type)&&(!is.null(age))
                &&(!is.null(pl085))&&(!is.null(month_at_work))) {

       varr <- linarr(Y=Y, Y_den=Y_den, id=id, age=age, pl085=pl085, 
                             month_at_work=month_at_work, weight=w_final, 
                             sort=sort, Dom=Dom, period=period, dataset=NULL,
                             order_quant=order_quant,  var_name="lin_arr") 
       varra <- linarr(Y=Y, Y_den=Y_den, id=id, age=age, pl085=pl085, 
                             month_at_work=month_at_work, weight=w_design, 
                             sort=sort, Dom=Dom, period=period, dataset=NULL,
                             order_quant=order_quant,  var_name="lin_arr") 

       Y1 <- merge(Y1, varr$lin, all.x=TRUE)
       Y1a <- merge(Y1a, varra$lin, all.x=TRUE)

       esti <- data.table("ARR", varr$value, NA)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       varr <- varra <- esti <- NULL
    }


  setkey(Y1, Y1sort)
  setkey(Y1a, Y1asort)
  Y1[, Y1sort:=NULL]
  Y1a[, Y1asort:=NULL]

  .SD <- lin_outp <- NULL
  if (outp_lin) lin_outp <- Y1[, c(-(3:5)-np), with=FALSE]

  Y2 <- Y1[, lapply(.SD, sum, na.rm = TRUE), by = c(names(Y1)[c(2:(5+np))]), .SDcols = names(Y1)[-(1:(5+np))]]
  Y2a <- Y1a[, lapply(.SD, sum, na.rm = TRUE), by = c(names(Y1a)[c(2:(5+np))]), .SDcols = names(Y1a)[-(1:(5+np))]]
  
  Y3 <- Y2[, c(-(1:(4+np))), with=FALSE]
  Y3a <- Y2a[, c(-(1:(4+np))), with=FALSE]
  
  idper <- period <- NULL
  if (np>0) period <- Y2[, c(1:np), with=FALSE]

  IDh <- Y2[, np+1, with=FALSE]
  H <- Y2[, np+2, with=FALSE]
  setnames(H, names(H), aH)

  PSU <- Y2[, np+3, with=FALSE]
  setnames(PSU, names(PSU), aPSU)

  w_final2 <- Y2[[np+4]]
  w_design2 <- Y2a[[np+4]]
  
  Y1 <- Y1a <- NULL
  Y2 <- Y2a <- NULL

  # Calibration

  res_outp <- variable <- NULL
  if (!is.null(X)) {
       if (np>0) IDh <- data.table(period, IDh)
       setnames(IDh, names(IDh), names(X_ID_household))
       X0 <- data.table(X_ID_household, ind_gr, q, g, X)
       D1 <- merge(IDh, X0, by=names(IDh))
       ind_gr <- D1[, np+2, with=FALSE]
       if (!is.null(period)) ind_gr <- data.table(D1[, names(periodX), with=FALSE], ind_gr)
       ind_period <- do.call("paste", c(as.list(ind_gr), sep="_"))
    
       lin1 <- lapply(split(Y3[, .I], ind_period), function(i) 
                      data.table(sar_nr=i, 
                             residual_est(Y=Y3[i],
                                          X=D1[i, (np+5):ncol(D1), with=FALSE],
                                          weight=w_design2[i],
                                          q=D1[i, np+3, with=FALSE])))
       Y4 <- rbindlist(lin1)
       setkeyv(Y4, "sar_nr")
       Y4[, sar_nr:=NULL]
       if (outp_res) res_outp <- data.table(IDh, PSU, w_final2, Y4)
   } else Y4 <- Y3

  var_est <- variance_est(Y=Y4, H=H, PSU=PSU, w_final=w_final2,
                          N_h=N_h, fh_zero=fh_zero, PSU_level=PSU_level,
                          period=period, dataset=NULL)   
  var_est <- transpos(var_est, is.null(period), "var_est", names(period))
  all_result <- var_est

    
  # Variance of HT estimator under current design
  var_cur_HT <- variance_est(Y=Y3a, H=H, PSU=PSU, w_final=w_design2, 
                             N_h=N_h, fh_zero=fh_zero, PSU_level=PSU_level,
                             period=period, dataset=NULL)                          
  var_cur_HT <- transpos(var_cur_HT, is.null(period), "var_cur_HT", names(period))
  all_result <- merge(all_result, var_cur_HT)
  var_est <- var_cur_HT <- NULL
  H <- PSU <- N_h <- NULL

  # Variance of HT estimator under SRS
  if (is.null(period)) {
           var_srs_HT <- var_srs(Y3a, w = w_design2)
       } else {
           period_agg <- unique(period)
           lin1 <- lapply(1:nrow(period_agg), function(i) {
                          per <- period_agg[i,][rep(1, nrow(Y3a)),]
                          ind <- (rowSums(per == period) == ncol(period))
                          data.table(period_agg[i,], 
                                     var_srs(Y3a[ind], w = w_design2[ind]))
                        })
           var_srs_HT <- rbindlist(lin1)
      }
  var_srs_HT <- transpos(var_srs_HT, is.null(period), "var_srs_HT", names(period))
  all_result <- merge(all_result, var_srs_HT)


  # Variance of calibrated estimator under SRS
   if (is.null(period)) {
           var_srs_ca <- var_srs(Y4, w = w_final2)
      } else {
           period_agg <- unique(period)
           lin1 <- lapply(1:nrow(period_agg), function(i) {
                          per <- period_agg[i,][rep(1, nrow(Y3a)),]
                          ind <- (rowSums(per == period) == ncol(period))
                          data.table(period_agg[i,], 
                                     var_srs(Y4[ind], w = w_final2[ind]))
                        })
           var_srs_ca <- rbindlist(lin1)
        }
  var_srs_ca <- transpos(var_srs_ca, is.null(period), "var_srs_ca", names(period))
  all_result <- merge(all_result, var_srs_ca)
  var_srs_HT <-  var_srs_ca <- NULL
  Y3a <- Y4 <- NULL

  estim <- data.table(estim)
  estim[, variable:=paste0("lin_", tolower(type))]
  nDom <- names(copy(Dom))
  if (!is.null(nDom)) estim[, (paste0(nDom,"at1at")):=lapply(nDom, function(x) paste(x, get(x), sep="."))]

  Dom <- estim[, "variable", with=F]
  if (!is.null(nDom)) Dom <- estim[, c("variable", paste0(nDom,"at1at")), with=F]

  estim$variable <- do.call("paste", c(as.list(Dom), sep="__"))
  estim[, variable:=str_replace_all(variable, "[ ]", ".")]
  if (!is.null(nDom)) estim[, (paste0(nDom,"at1at")):=NULL]
  
  if (nrow(all_result[var_est < 0])>0) stop("Estimation of variance are negative!")
  
  variables <- "variable"
  if (!is.null(period)) variables <- c(variables, names(period))
  setkeyv(estim, variables)
  setkeyv(all_result, variables)
  all_result <- merge(estim, all_result, all=TRUE)
  
  all_result[, variable:=NULL]
  deff_sam <- deff_est <- deff <- n_eff <- var_est2 <- NULL
  se <- rse <- cv <- absolute_margin_of_error <- NULL
  relative_margin_of_error <- CI_lower <- CI_upper <- NULL

  # Design effect of sample design
  all_result[, deff_sam:=var_cur_HT / var_srs_HT]
  
  # Design effect of estimator
  all_result[, deff_est:= var_est / var_cur_HT]
  
  # Overall effect of sample design and estimator
  all_result[, deff:= deff_sam * deff_est]
 
  all_result[, var_est2:=var_est]
  all_result[xor(is.na(var_est2), var_est2 < 0), var_est2:=0]
  all_result[, se:=sqrt(var_est2)]
  all_result[xor(is.na(var_est2), var_est2 < 0), se:=NA]
  all_result[(value!=0) & (!is.nan(value)), rse:= se/value]
  all_result[value==0 | is.nan(value), rse:=NA]
  all_result[, cv:= rse*100]

  tsad <- qnorm(0.5*(1+confidence))
  all_result[, absolute_margin_of_error:= tsad*se]
  all_result[, relative_margin_of_error:= tsad*cv]
  all_result[, CI_lower:= value - tsad*se]
  all_result[, CI_upper:= value + tsad*se]
  
  setnames(all_result, "var_est", "var")
  
  setkeyv(all_result, c(nDom, names(period)))

  if (!is.null(c(nDom, period))) { all_result <- merge(all_result, nhs, all=TRUE)
                         } else { all_result[, respondent_count:=nhs$respondent_count]
                                  all_result[, pop_size:=nhs$pop_size]
                                  all_result[, n_nonzero:=nhs$n_nonzero]} 

  variabl <- c("respondent_count", "n_nonzero", "pop_size", 
                      "value", "value_eu", "var", "se", "rse",
                      "cv", "absolute_margin_of_error",
                      "relative_margin_of_error", "CI_lower",  
                      "CI_upper", "var_srs_HT", "var_cur_HT", 
                      "var_srs_ca", "deff_sam", "deff_est", "deff")

  type <- "type"
  if (!is.null(period)) type <- c(type, names(period))
  setkeyv(all_result, c(type, nDom))
  list(lin_out = lin_outp,
       res_out = res_outp,
       all_result = all_result[, c(type, nDom, variabl), with=FALSE])
}