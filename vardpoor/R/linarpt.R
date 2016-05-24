#***********************************************************************************************
#***********************************************************************************************
#***********************************************************************************************
#***                                                                                         ***
#***                                                                                         ***
#***                    LINEARIZATION OF THE AT-RISK-OF-POVERTY THRESHOLD                    ***
#***                                                                                         ***
#***                                                                                         ***
#***********************************************************************************************
#***********************************************************************************************
#***********************************************************************************************

linarpt <- function(Y, id = NULL, weight = NULL, sort = NULL, 
        Dom = NULL, period=NULL, dataset = NULL, percentage = 60,
        order_quant=50, var_name="lin_arpt") {

   ## initializations
   if (min(dim(as.data.frame(var_name))==1)!=1) {
       stop("'var_name' must have defined name of the linearized variable")}

   # check 'p'
   p <- percentage
   if(length(p) != 1 |  any(!is.numeric(p) | p < 0 | p > 100)) {
          stop("'percentage' must be a numeric value in [0, 100]")  }

   # check 'order_quant'
   oq <- order_quant
   if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]")  }

   if(!is.null(dataset)) {
       dataset <- data.table(dataset)
       if (checker(Y, dataset, "Y")) Y <- dataset[, Y, with=FALSE] 

       if(!is.null(id)) {
          if (checker(id, dataset, "id")) id <- dataset[, id, with=FALSE]}

       if(!is.null(weight)) {
           if (checker(weight, dataset, "weight")) weight <- dataset[, weight, with=FALSE] }

       if(!is.null(sort)) {
           if (checker(sort, dataset, "sort")) sort <- dataset[, sort, with=FALSE] }

       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE] }

       if(!is.null(Dom)) {
            if (checker(Dom,dataset,"Dom")) Dom <- dataset[, Dom, with=FALSE] }
      }

   # check vectors
   # Y
   Y <- data.frame(Y)
   n <- nrow(Y)
   if (ncol(Y) != 1) stop("'Y' must be a vector or 1 column data.frame, matrix, data.table")
   Y <- Y[,1]
   if(!is.numeric(Y)) stop("'Y' must be numerical")
   if (any(is.na(Y))) stop("'Y' has unknown values")

   # weight
   weight <- data.frame(weight)
   if (nrow(weight) != n) stop("'weight' must be the same length as 'Y'")
   if (ncol(weight) != 1) stop("'weight' must be vector or 1 column data.frame, matrix, data.table")
   weight <- weight[,1]
   if (!is.numeric(weight)) stop("'weight' must be numerical")
   if (any(is.na(weight))) stop("'weight' has unknown values")

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
       period[, (names(period)):=lapply(.SD, as.character)]
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
    ind0 <- rep.int(1,n)
    period_agg <- period1 <- NULL
    if (!is.null(period)) { period1 <- copy(period)
                            period_agg <- data.table(unique(period))
                        } else period1 <- data.table(ind=ind0)
    period1_agg <- data.table(unique(period1))

    # ARPT by domain (if requested)  

    quantile <- incPercentile(Y = Y, weights = weight, sort = sort,
                              Dom = Dom, period=period,
                              k = order_quant, dataset = NULL)
    quantile <- data.table(quantile)
    setnames(quantile, names(quantile)[ncol(quantile)], "quantile")
    if (ncol(quantile)>1) setkeyv(quantile, head(names(quantile), -1))
    threshold <- copy(quantile)
    threshold[, threshold:=p/100 * quantile]
    threshold[, quantile:=NULL]
    
    arpt_id <- id       
    if (!is.null(period)) arpt_id <- data.table(arpt_id, period)
    
    if(!is.null(Dom)) {
        Dom_agg <- data.table(unique(Dom))
        setkeyv(Dom_agg, names(Dom_agg))

        arpt_m <- copy(arpt_id)
        for(i in 1:nrow(Dom_agg)) {
              g <- c(var_name, paste(names(Dom), as.matrix(Dom_agg[i,]), sep = "."))
              var_nams <- do.call(paste, as.list(c(g, sep="__")))
              ind <- as.integer(rowSums(Dom == Dom_agg[i,][ind0,]) == ncol(Dom))
              arpt_l <- lapply(1:nrow(period1_agg), function(j) {
                               if (!is.null(period)) { 
                                       rown <- cbind(period_agg[j], Dom_agg[i])
                                       } else rown <- Dom_agg[i]

                               setkeyv(rown, names(rown))
                               rown2 <- copy(rown)
                               rown <- merge(rown, quantile, all.x=TRUE)
                               ind2 <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
                               
                               arptl <- arptlinCalc(inco=Y[ind2], 
                                                    ids=arpt_id[ind2],
                                                    wght=weight[ind2],
                                                    indicator=ind[ind2], 
                                                    order_quan=order_quant,
                                                    quant_val=rown[["quantile"]],
                                                    percentag=p)
                               })
             arptl <- rbindlist(arpt_l)
             setnames(arptl, names(arptl), c(names(arpt_id), var_nams))
             arpt_m <- merge(arpt_m, arptl, all.x=TRUE, by=names(arpt_id))
          }
      } else { arptl <- lapply(1:nrow(period1_agg), function(j) {
                           if (!is.null(period)) { 
                                         rown <- period_agg[j]
                                         setkeyv(rown, names(rown))
                                         rown <- merge(rown, quantile, all.x=TRUE)
                                       } else rown <- quantile
                           ind2 <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
 
                           arptl <- arptlinCalc(inco=Y[ind2], 
                                                ids=arpt_id[ind2],
                                                wght=weight[ind2],
                                                indicator=ind0[ind2], 
                                                order_quan=order_quant,
                                                quant_val=rown[["quantile"]],
                                                percentag=p)
                       })
               arpt_m <- rbindlist(arptl)
               setnames(arpt_m, names(arpt_m), c(names(arpt_id), var_name))
            }
    arpt_m[is.na(arpt_m)] <- 0
    setkeyv(arpt_m, names(arpt_id))
    return(list(quantile=quantile, value=threshold, lin=arpt_m)) 
 }

    ## workhorse
arptlinCalc <- function(inco, ids, wght, indicator, order_quan, quant_val, percentag) {
    wt <- wght * indicator
    N <- sum(wt); # Estimated (sub)population size

    # h=S/N^(1/5) 
    h <- sqrt((sum(wght*inco*inco)-sum(wght*inco)*sum(wght*inco)/sum(wght))/sum(wght))/exp(0.2*log(sum(wght))) 

    u <- (quant_val-inco)/h
    vect_f <- exp(-(u^2)/2)/sqrt(2*pi)
    f_quant <- sum(vect_f*wt)/(N*h) # Estimate of F'(quantile)

 #****************************************************************************************
 #*                    LINEARIZED VARIABLE OF THE POVERTY THRESHOLD                      *
 #****************************************************************************************
    lin <- -(percentag/100)*(1/N)*indicator*((inco<=quant_val)-order_quan/100)/f_quant
    lin_id <- data.table(ids, lin)
    return(lin_id)
}

