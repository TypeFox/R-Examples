#******************************************************************************************
#******************************************************************************************
#******************************************************************************************
#                                                                                
#                                                                                
#                       LINEARIZATION OF THE AGGREGATE REPLACEMENT RATIO     
#
#
#******************************************************************************************
#******************************************************************************************
#******************************************************************************************

linarr <- function(Y, Y_den, id=NULL, age, pl085, month_at_work, weight=NULL,  sort=NULL, 
                             Dom=NULL, period=NULL, dataset = NULL, 
                             order_quant=50, var_name="lin_arr") {
 
   ## initializations
   if (min(dim(data.table(var_name))==1)!=1) {
       stop("'var_name' must have defined one name of the linearized variable")}

   # check 'order_quant'
   oq <- order_quant
   if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]")  }

   if(!is.null(dataset)) {
       dataset <- data.table(dataset)
       if (checker(Y, dataset, "Y"))  Y <- dataset[, Y, with=FALSE] 
       if (checker(Y_den, dataset, "Y_den"))  Y_den <- dataset[, Y_den, with=FALSE] 

       if(!is.null(id)) {
          if (checker(id, dataset, "id")) id <- dataset[, id, with=FALSE] }

       if(!is.null(age)) {
           if (checker(age, dataset, "age")) age <- dataset[, age, with=FALSE] }

       if(!is.null(pl085)) {
           if (checker(pl085, dataset, "pl085")) pl085 <- dataset[, pl085, with=FALSE] }

       if(!is.null(month_at_work)) {
           if (checker(month_at_work, dataset, "month_at_work")) {
                  month_at_work <- dataset[, month_at_work, with=FALSE] }}

       if(!is.null(weight)) {
           if (checker(weight, dataset, "weight")) weight <- dataset[, weight, with=FALSE] }

       if(!is.null(sort)) {
           if (checker(sort, dataset, "sort")) sort <- dataset[, sort, with=FALSE] }

       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE] }

       if (!is.null(Dom)) {
            if (checker(Dom, dataset, "Dom")) Dom <- dataset[, Dom, with=FALSE] }
    }

   # check vectors
   # Y
   Y <- data.frame(Y)
   n <- nrow(Y)
   if (ncol(Y) != 1) stop("'Y' must be vector or 1 column data.frame, matrix, data.table")
   Y <- Y[,1]
   if(!is.numeric(Y)) stop("'Y' must be numerical")
   if (any(is.na(Y))) stop("'Y' has unknown values")

   Y_den <- data.frame(Y_den)
   if (ncol(Y_den) != 1) stop("'Y_den' must be vector or 1 column data.frame, matrix, data.table")
   if (nrow(Y_den) != n) stop("'Y_den' must be the same length as 'Y'")
   Y_den <- Y_den[,1]
   if(!is.numeric(Y_den)) stop("'Y_den' must be numerical")
   if (any(is.na(Y_den))) stop("'Y_den' has unknown values")

   # weight
   weight <- data.frame(weight)
   if (nrow(weight) != n) stop("'weight' must be the same length as 'Y'")
   if (ncol(weight) != 1) stop("'weight' must be vector or 1 column data.frame, matrix, data.table")
   weight <- weight[, 1]
   if (!is.numeric(weight)) stop("'weight' must be numerical")
   if (any(is.na(weight))) stop("'weight' has unknown values")
   
   # age
   age <- data.frame(age)
   if (nrow(age) != n) stop("'age' must be the same length as 'Y'")
   if (ncol(age) != 1) stop("'age' must be vector or 1 column data.frame, matrix, data.table")
   age <- age[, 1]
   if (!is.numeric(age)) stop("'age' must be numerical")
   if (any(is.na(age))) stop("'age' has unknown values")
   
   # pl085
   pl085 <- data.frame(pl085)
   if (nrow(pl085) != n) stop("'pl085' must be the same length as 'Y'")
   if (ncol(pl085) != 1) stop("'pl085' must be vector or 1 column data.frame, matrix, data.table")
   pl085 <- pl085[, 1]
   if (!is.numeric(pl085)) stop("'pl085' must be numerical")
   if (any(is.na(pl085))) stop("'pl085' has unknown values")

   # month_at_work
   month_at_work <- data.frame(month_at_work)
   if (nrow(month_at_work) != n) stop("'month_at_work' must be the same length as 'Y'")
   if (ncol(month_at_work) != 1) stop("'month_at_work' must be vector or 1 column data.frame, matrix, data.table")
   month_at_work <- month_at_work[, 1]
   if (!is.numeric(pl085)) stop("'month_at_work' must be numerical")
   if (any(is.na(pl085))) stop("'month_at_work' has unknown values")

   # sort
   if (!is.null(sort)) {    
        sort <- data.frame(sort)
        if (length(sort) != n) stop("'sort' must have the same length as 'Y'")
        if (ncol(sort) != 1) stop("'sort' must be vector or 1 column data.frame, matrix, data.table")
        sort <- sort[, 1]
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
   
   # id
   if (is.null(id)) id <- 1:n
   id <- data.table(id)
   if (any(is.na(id))) stop("'id' has unknown values")
   if (ncol(id) != 1) stop("'id' must be 1 column data.frame, matrix, data.table")
   if (nrow(id) != n) stop("'id' must be the same length as 'Y'")
   if (is.null(names(id))||(names(id)=="id")) setnames(id,names(id),"ID")
   if (is.null(period)){ if (any(duplicated(id))) stop("'id' are duplicate values") 
                       } else {
                          id1 <- data.table(period, id)
                          if (any(duplicated(id1))) stop("'id' by period are duplicate values")
                         }

   # Dom     
   if (!is.null(Dom)) {
             Dom <- data.table(Dom)
             if (any(duplicated(names(Dom)))) 
                 stop("'Dom' are duplicate column names: ", 
                      paste(names(Dom)[duplicated(names(Dom))], collapse = ","))
             if (is.null(names(Dom))) stop("'Dom' must be colnames")
             if (nrow(Dom) != n) stop("'Dom' must be the same length as 'Y'")
             Dom[, (names(Dom)):=lapply(.SD, as.character)]
       }

    ## computations
    ind0 <- rep.int(1, n)
    period_agg <- period1 <- NULL
    if (!is.null(period)) { period1 <- copy(period)
                            period_agg <- data.table(unique(period))
                        } else period1 <- data.table(ind=ind0)
    period1_agg <- data.table(unique(period1))

    # RMI by domain (if requested)
    age_65_74pl <- data.table(age_65_74pl=as.integer(65 <= age & age <= 74 & pl085==12))
    age_50_59mo <- data.table(age_50_59mo=as.integer(50 <= age & age <= 59 & month_at_work==12))


    if (!is.null(Dom)) { age_65_74pl <- data.table(age_65_74pl, Dom)
                         age_50_59mo <- data.table(age_50_59mo, Dom) }

    quantile1 <- incPercentile(Y = Y, weights = weight,
                                sort = sort, Dom = age_65_74pl,
                                period = period, k = order_quant,
                                dataset = NULL)
    quantile2 <- incPercentile(Y = Y_den, weights = weight,
                                sort = sort, Dom = age_50_59mo,
                                period = period, k = order_quant,
                                dataset = NULL)
    quantile1 <- data.table(quantile1)
    quantile2 <- data.table(quantile2)

    quantile1 <- quantile1[age_65_74pl==1][, age_65_74pl:=NULL]
    quantile2 <- quantile2[age_50_59mo==1][, age_50_59mo:=NULL]
    setnames(quantile1, names(quantile1)[ncol(quantile1)], "quantile_65_74pl")
    setnames(quantile2, names(quantile2)[ncol(quantile2)], "quantile_50_59mo")
    sk <- length(names(quantile2))-1
    if (sk > 0) {
               quantile <- merge(quantile1, quantile2, all=TRUE,
                                 by=names(quantile1)[1:sk])
        } else quantile <- data.table(quantile1, quantile2)
    
    arr_id <- id
    quantile1 <- quantile2 <- NULL
    age_65_74pl <- age_65_74pl[["age_65_74pl"]]
    age_50_59mo <- age_50_59mo[["age_50_59mo"]]
    if (!is.null(period)) arr_id <- data.table(arr_id, period)

  if (!is.null(Dom)) {
        Dom_agg <- data.table(unique(Dom))
        setkeyv(Dom_agg, names(Dom_agg))
          
        arr_v <- c()
        arr_m <- copy(arr_id)
        for(i in 1:nrow(Dom_agg)) {
              g <- c(var_name, paste(names(Dom), as.matrix(Dom_agg[i,]), sep = "."))
              var_nams <- do.call(paste, as.list(c(g, sep="__")))
              ind <- as.integer(rowSums(Dom == Dom_agg[i,][ind0,]) == ncol(Dom))

              arrl <- lapply(1:nrow(period1_agg), function(j) {
                              if (!is.null(period)) { 
                                       rown <- cbind(period_agg[j], Dom_agg[i])
                                       setkeyv(rown, names(rown))
                                       rown2 <- copy(rown)
                                       rown <- merge(rown, quantile, all.x=TRUE)
                                     } else {rown <- quantile[i]
                                             rown2 <- Dom_agg[i] }

                               indj <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))

                               arr_l <- arrlinCalc(Y_num=Y[indj],
                                                   Y_den=Y_den[indj],
                                                   ids=arr_id[indj],
                                                   wght=weight[indj],
                                                   indicator=ind[indj],
                                                   order_quants=order_quant,
                                                   age_65_74pl=age_65_74pl[indj],
                                                   age_50_59mo=age_50_59mo[indj],
                                                   quant_65_74pls=rown[["quantile_65_74pl"]],
                                                   quant_50_59mon=rown[["quantile_50_59mo"]])

                      list(arr=data.table(rown2, arr=arr_l$arr_val), lin=arr_l$lin)
                      })
                 arrs <- rbindlist(lapply(arrl, function(x) x[[1]]))
                 arrlin <- rbindlist(lapply(arrl, function(x) x[[2]]))

                 setnames(arrlin, names(arrlin), c(names(arr_id), var_nams))
                 arr_m <- merge(arr_m, arrlin, all.x=TRUE, by=names(arr_id))
                 arr_v <- rbind(arr_v, arrs) 
           }
     } else { arrl <- lapply(1:nrow(period1_agg), function(j) {  
         
             if (!is.null(period)) { rown <- period_agg[j]
                                     rown <- merge(rown, quantile, all.x=TRUE, 
                                                   by=names(rown))
                                 } else rown <- quantile
                           ind2 <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
      
                           arr_l <- arrlinCalc(Y_num=Y[ind2],
                                                       Y_den=Y_den[ind2],
                                                       ids=arr_id[ind2],
                                                       wght=weight[ind2],
                                                       indicator=ind0[ind2],
                                                       order_quants=order_quant,
                                                       age_65_74pl=age_65_74pl[ind2],
                                                       age_50_59mo=age_50_59mo[ind2],
                                                       quant_65_74pls=rown[["quantile_65_74pl"]],
                                                       quant_50_59mon=rown[["quantile_50_59mo"]])

                          if (!is.null(period)) { 
                                   arrs <- data.table(period_agg[j], arr=arr_l$arr_val)
                             } else arrs <- data.table(arr=arr_l$arr_val)
                          list(arr=arrs, lin=arr_l$lin)
                       })
               arr_v <- rbindlist(lapply(arrl, function(x) x[[1]]))
               arr_m <- rbindlist(lapply(arrl, function(x) x[[2]]))
               setnames(arr_m, names(arr_m), c(names(arr_id), var_name))
            } 
     arr_m[is.na(arr_m)] <- 0
     setkeyv(arr_m, names(arr_id))
     return(list(value=arr_v, lin=arr_m))
}



## workhorse
arrlinCalc <- function(Y_num, Y_den, ids, wght, indicator, order_quants,
                       age_65_74pl, age_50_59mo, quant_65_74pls, quant_50_59mon) {

    dom1 <- (age_65_74pl==1) * indicator
    dom2 <- (age_50_59mo==1) * indicator
 
   # Size of the domains
    N1 <- sum(wght * dom1)   
    N2 <- sum(wght * dom2) 
				
    arr_val <- quant_65_74pls/quant_50_59mon  # Estimated aggregate replacement ratio

    # Bandwith parameter - h=S/N^(1/5) (calculated over the whole population) 

    h1 <- sqrt((sum(wght*Y_num*Y_num)-sum(wght*Y_num)*sum(wght*Y_num)/sum(wght))/sum(wght))/exp(0.2*log(sum(wght))) 
    h2 <- sqrt((sum(wght*Y_den*Y_den)-sum(wght*Y_den)*sum(wght*Y_den)/sum(wght))/sum(wght))/exp(0.2*log(sum(wght))) 

    #---- 1. Linearization of the median income of people aged below 65 ----

    u1 <- (quant_65_74pls - Y_num) * dom1/h1
    vect_f1 <- exp(-(u1^2)/2)/sqrt(2 * pi)
    f_quant1 <- sum(vect_f1 * wght)/(N1 * h1)   # Estimate of F'(quantile)

    lin_quant_65_74pl <- -(1/N1) * dom1 * ((Y_num<=quant_65_74pls) - order_quants/100)/f_quant1  # Linearized variable

    #---- 2. Linearization of the median income of people aged above 65 -----

    u2 <- (quant_50_59mon - Y_den) * dom2/h2
    vect_f2 <- exp(-(u2^2)/2)/sqrt(2 * pi)
    f_quant2 <- sum(vect_f2 * wght)/(N2 * h2)   # Estimate of F'(quantile)

    lin_quant_50_59mon <- -(1/N2) * dom2 * ((Y_den<=quant_50_59mon) - order_quants/100)/f_quant2  # Linearized variable

   #********************************************************************************
   #         3. Linearization of the relative median income ratio                  *
   #********************************************************************************
    lin <- (quant_50_59mon * lin_quant_65_74pl - quant_65_74pls * lin_quant_50_59mon)/(quant_50_59mon*quant_50_59mon)

    lin_id <- data.table(ids, lin)
    return(list(arr_val=arr_val, lin=lin_id))
}

