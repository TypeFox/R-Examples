#******************************************************************************************
#******************************************************************************************
#******************************************************************************************
#***                                                                                    ***
#***                                                                                    ***
#***                    LINEARIZATION OF THE AT-RISK-OF-POVERTY RATE                    ***
#***                                                                                    ***
#***                                                                                    ***
#******************************************************************************************
#******************************************************************************************
#******************************************************************************************

linarpr <- function(Y, id=NULL, weight=NULL, Y_thres=NULL, wght_thres=NULL,
                 sort=NULL, Dom=NULL, period=NULL, dataset = NULL, percentage=60, 
                 order_quant=50, var_name="lin_arpr") {
 
   ## initializations
   if (min(dim(data.table(var_name))==1)!=1) {
       stop("'var_name' must have defined one name of the linearized variable")}

   # check 'p'
   p <- percentage
   if(length(p) != 1 |  any(!is.numeric(p) | p < 0 | p > 100)) {
          stop("'percentage' must be a numeric value in [0, 100]") }

   # check 'order_quant'
   oq <- order_quant
   if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]") }

   if(!is.null(dataset)) {
       dataset <- data.table(dataset)
       if (checker(Y, dataset, "Y"))  Y <- dataset[, Y, with=FALSE] 

       if(!is.null(id)) {
          if (checker(id, dataset, "id")) id <- dataset[, id, with=FALSE]}

       if(!is.null(weight)) {
           if (checker(weight, dataset, "weight")) weight <- dataset[, weight, with=FALSE] }

       if(!is.null(Y_thres)) {
           if (checker(Y_thres, dataset, "Y_thres")) Y_thres <- dataset[, Y_thres, with=FALSE] }
     
       if(!is.null(wght_thres)) {
           if (checker(wght_thres, dataset, "wght_thres")) wght_thres <- dataset[, wght_thres, with=FALSE] }

       if(!is.null(sort)) {
           if (checker(sort, dataset, "sort")) sort <- dataset[, sort, with=FALSE] }

       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE]}

       if (!is.null(Dom)) {
            if (checker(Dom, dataset, "Dom")) Dom <- data.table(dataset[, Dom, with=FALSE])}
    }

   # check vectors
   # Y
   Y <- data.frame(Y)
   n <- nrow(Y)
   if (ncol(Y) != 1) stop("'Y' must be vector or 1 column data.frame, matrix, data.table")
   Y <- Y[,1]
   if(!is.numeric(Y)) stop("'Y' must be numerical")
   if (any(is.na(Y))) stop("'Y' has unknown values")

   # weight
   weight <- data.frame(weight)
   if (nrow(weight) != n) stop("'weight' must be the same length as 'inc'")
   if (ncol(weight) != 1) stop("'weight' must be vector or 1 column data.frame, matrix, data.table")
   weight <- weight[,1]
   if (!is.numeric(weight)) stop("'weight' must be numerical")
   if (any(is.na(weight))) stop("'weight' has unknown values")

   # Y_thres
   if (is.null(Y_thres)) Y_thres <- Y
   Y_thres <- data.frame(Y_thres)
   if (!is.null(Y_thres)) {
        if (nrow(Y_thres) != n) stop("'Y_thres' must be the same length as 'Y'")
        if (ncol(Y_thres) != 1) stop("'Y_thres' must be vector or 1 column data.frame, matrix, data.table")
        Y_thres <- Y_thres[,1]
        if(!is.numeric(Y_thres)) stop("'Y_thres' must be numerical")
        if (any(is.na(Y_thres))) stop("'Y_thres' has unknown values") 
     } 

   # wght_thres
   if (is.null(wght_thres)) wght_thres <- weight
   wght_thres <- data.frame(wght_thres)
   if (nrow(wght_thres) != n) stop("'wght_thres' must be the same length as 'Y'")
   if (ncol(wght_thres) != 1) stop("'wght_thres' must be vector or 1 column data.frame, matrix, data.table")
   wght_thres <- wght_thres[,1]
   if (!is.numeric(wght_thres)) stop("'wght_thres' must be numerical")
   if (any(is.na(wght_thres))) stop("'wght_thres' has unknown values") 
    
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

    # ARPR by domain (if requested)
    quantile <- incPercentile(Y = Y_thres,
                               weights = wght_thres,
                               sort = sort, Dom = NULL,
                               period = period,
                               k = order_quant,
                               dataset = NULL)
    quantile <- data.table(quantile)
    setnames(quantile, names(quantile)[ncol(quantile)], "quantile")
    if (ncol(quantile)>1) setkeyv(quantile, head(names(quantile), -1))
    threshold <- copy(quantile)
    threshold[, threshold:=p/100 * quantile]
    threshold[, quantile:=NULL]

    arpr_id <- id
    if (!is.null(period)) arpr_id <- data.table(arpr_id, period)

    if (!is.null(Dom)) {
        Dom_agg <- data.table(unique(Dom))
        setkeyv(Dom_agg, names(Dom_agg))
          
        arpr_v <- c()
        arpr_m <- copy(arpr_id)
        for(i in 1:nrow(Dom_agg)) {
              g <- c(var_name, paste(names(Dom), as.matrix(Dom_agg[i,]), sep = "."))
              var_nams <- do.call(paste, as.list(c(g, sep="__")))
              ind <- as.integer(rowSums(Dom == Dom_agg[i,][ind0,]) == ncol(Dom))

              arprl <- lapply(1:nrow(period1_agg), function(j) {
                               if (!is.null(period)) { 
                                       rown <- cbind(period_agg[j], Dom_agg[i])
                                       setkeyv(rown, names(rown))
                                       rown2 <- copy(rown)
                                       rown <- merge(rown, quantile, all.x=TRUE)
                                     } else {rown <- quantile
                                             rown2 <- Dom_agg[i] }

                               indj <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
      
                               arpr_l <- arprlinCalc(Y1=Y[indj],
                                                     ids=arpr_id[indj],
                                                     wght1=weight[indj],
                                                     indicator=ind[indj],
                                                     Y_thresh=Y_thres[indj],
                                                     wght_thresh=wght_thres[indj],
                                                     percent=p,
                                                     order_quants=order_quant,
                                                     quant_val=rown[["quantile"]])
                      list(arpr=data.table(rown2, arpr=arpr_l$rate_val_pr), lin=arpr_l$lin)
                      })
                 arprs <- rbindlist(lapply(arprl, function(x) x[[1]]))
                 arprlin <- rbindlist(lapply(arprl, function(x) x[[2]]))

                 setnames(arprlin, names(arprlin), c(names(arpr_id), var_nams))
                 arpr_m <- merge(arpr_m, arprlin, all.x=TRUE, by=names(arpr_id))
                 arpr_v <- rbind(arpr_v, arprs) 
           }
     } else { arprl <- lapply(1:nrow(period1_agg), function(j) {
                           if (!is.null(period)) { 
                                         rown <- period_agg[j]
                                         rown <- merge(rown, quantile, all.x=TRUE,
                                                       by=names(rown))
                                       } else rown <- quantile
                           ind2 <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
      
                           arpr_l <- arprlinCalc(Y1=Y[ind2],
                                                 ids=arpr_id[ind2],
                                                 wght1=weight[ind2],
                                                 indicator=ind0[ind2],
                                                 Y_thresh=Y_thres[ind2],
                                                 wght_thresh=wght_thres[ind2],
                                                 percent=p,
                                                 order_quants=order_quant,
                                                 quant_val=rown[["quantile"]]) 
                          if (!is.null(period)) { 
                                   arprs <- data.table(period_agg[j], arpr=arpr_l$rate_val_pr)
                             } else arprs <- data.table(arpr=arpr_l$rate_val_pr)
                          list(arpr=arprs, lin=arpr_l$lin)
                       })
               arpr_v <- rbindlist(lapply(arprl, function(x) x[[1]]))
               arpr_m <- rbindlist(lapply(arprl, function(x) x[[2]]))
               setnames(arpr_m, names(arpr_m), c(names(arpr_id), var_name))
            } 
     arpr_m[is.na(arpr_m)] <- 0
     setkeyv(arpr_m, names(arpr_id))
     return(list(quantile=quantile, threshold=threshold, value=arpr_v, lin=arpr_m))
}




## workhorse
arprlinCalc <- function(Y1, ids, wght1, indicator, Y_thresh, wght_thresh, percent, order_quants=NULL, quant_val) {

    inc2 <- Y_thresh
    wght2 <- wght_thresh
  
    wt <- indicator * wght1
    thres_val <- percent/100 * quant_val
    N0 <- sum(wght2)   # Estimated whole population size 
    N <- sum(wt)       # Estimated (sub)population size

    poor <- (Y1<=thres_val)
    rate_val <- sum(wt*poor)/N  # Estimated poverty rate */
    rate_val_pr <- 100*rate_val

    h <- sqrt((sum(wght2*inc2*inc2)-sum(wght2*inc2)*sum(wght2*inc2)/sum(wght2))/sum(wght2))/exp(0.2*log(sum(wght2))) 
    # h=S/N^(1/5)

    #---- 1. Linearization of the poverty threshold ----

    u1 <- (quant_val-inc2)/h
    vect_f1 <- exp(-(u1^2)/2)/sqrt(2*pi)
    f_quant1 <- sum(vect_f1*wght2)/(N0*h)   # Estimate of F'(quantile)

    lin_thres <- -(percent/100)*(1/N0)*((inc2<=quant_val)-order_quants/100)/f_quant1  # Linearized variable

    #---- 2. Linearization of the poverty rate -----

    u2 <-   (thres_val-Y1)/h
    vect_f2 <- exp(-(u2^2)/2)/sqrt(2*pi)
    f_quant2 <- sum(vect_f2*wt)/(N*h)      # Estimate of F'(beta*quantile)

 #****************************************************************************************
 #                       LINEARIZED VARIABLE OF THE POVERTY RATE (IN %)                  *
 #****************************************************************************************
    lin <- 100*((1/N)*indicator*((Y1<=thres_val)-rate_val)+f_quant2*lin_thres)

    lin_id <- data.table(ids,lin)
    return(list(rate_val=rate_val, rate_val_pr=rate_val_pr, lin=lin_id))
}

