# ************************************************************************
# ************************************************************************
# ************************************************************************
# ***                                                                  ***
# ***                                                                  ***
# ***            LINEARIZATION OF THE GINI COEFFICIENT II              ***
# ***                                                                  ***
# ***                                                                  ***
# ************************************************************************
# ************************************************************************
# ************************************************************************
# ************************************************************************


lingini2 <- function(Y, id = NULL, weight = NULL, sort = NULL, 
                     Dom = NULL, period=NULL, dataset = NULL,
                     var_name="lin_gini2") {

   ## initializations
   if (min(dim(as.data.frame(var_name))==1)!=1) {
       stop("'var_name' must have defined name of the linearized variable")}

   if (!is.null(dataset)) {
        dataset <- data.table(dataset)
        if (checker(Y,dataset,"Y")) Y <- dataset[, Y, with=FALSE] 

        if (!is.null(id)) {
             if (checker(id,dataset,"id")) id <- dataset[, id, with=FALSE]}

        if (!is.null(weight)) {
             if (checker(weight,dataset,"weight")) weight <- dataset[, weight, with=FALSE] }

        if (!is.null(sort)) {
             if (checker(sort,dataset,"sort")) sort <- dataset[, sort, with=FALSE] }

        if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) period <- dataset[, period, with=FALSE]}

        if (!is.null(Dom)) {
             if (checker(Dom,dataset,"Dom")) Dom <- dataset[, Dom, with=FALSE] }
    } 

   # check vectors
   # Y
   Y <- data.frame(Y)
   n <- nrow(Y)
   if (ncol(Y) != 1) stop("'Y' must be a vector or 1 column data.frame, matrix, data.table")
   Y <- Y[,1]
   if(!is.numeric(Y)) stop("'Y' must be a numerical")
   if (any(is.na(Y))) stop("'Y' has unknown values")
 
   # weight
   weight <- data.frame(weight)
   if (is.null(weight)) weight <- data.frame(rep.int(1, n))
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
                          id1 <- NULL
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


   # Gini by domain (if requested)
   gini_id <- id
   if (!is.null(period)) gini_id <- data.table(period, gini_id)

   if (!is.null(Dom)) {
        Dom_agg <- data.table(unique(Dom))
        setkeyv(Dom_agg, names(Dom_agg))

        Gini <- c()
        gini_m <- copy(gini_id)
        for(i in 1:nrow(Dom_agg)) {
            g <- c(var_name, paste(names(Dom), as.matrix(Dom_agg[i,]), sep = "."))
            var_nams <- do.call(paste, as.list(c(g, sep="__")))
            indi <- (rowSums(Dom == Dom_agg[i,][ind0,]) == ncol(Dom))

            gini_l <- lapply(1:nrow(period1_agg), function(j) {
               indj <- ((rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))&(indi))
               if (!is.null(period)) { rown <- cbind(period_agg[j], Dom_agg[i])
                                     } else rown <- Dom_agg[i] 
               
               ginil <- lingini2Calc(x=Y[indj], ids=gini_id[indj],
                                     weights=weight[indj], sort=sort[indj])
               list(data.table(rown, ginil$Gini), ginil$lin)
              })

            giniv <- rbindlist(lapply(gini_l, function(x) x[[1]]))
            ginilin <- rbindlist(lapply(gini_l, function(x) x[[2]]))

            setnames(ginilin, names(ginilin), c(names(gini_id), var_nams))
            gini_m <- merge(gini_m, ginilin, all=TRUE, by=names(gini_id))
            Gini <- rbind(Gini, giniv)
         }
     } else { gini_l <- lapply(1:nrow(period1_agg), function(j) {
                           indj <- (rowSums(period1 == period1_agg[j,][ind0,]) == ncol(period1))
                           ginil <- lingini2Calc(x=Y[indj], ids=gini_id[indj],
                                                 weights=weight[indj], sort=sort[indj])                                                 
                           if (!is.null(period)) {
                                  list(data.table(period_agg[j], ginil$Gini), ginil$lin)
                                }  else ginil
                         })
            Gini <- rbindlist(lapply(gini_l, function(x) x[[1]]))
            gini_m <- rbindlist(lapply(gini_l, function(x) x[[2]]))
            setnames(gini_m, names(gini_m), c(names(gini_id), var_name))   
     }
    gini_m[is.na(gini_m)] <- 0
    setkeyv(gini_m, names(gini_id))
    return(list(value=Gini, lin=gini_m))
}


# workhorse 
lingini2Calc <- function(x, ids, weights = NULL, sort = NULL) { 

    # sort values and weights 
    order <- if(is.null(sort)) order(x) else order(x, sort) 
    x <- x[order]  # order values 
    ids <- ids[order]  # order values 
    if (is.null(weights)) { weights <- rep.int(1, length(x))  # equal weights 
     } else weights <- weights[order]  # order weights 

    ## calculations 
    wx <- weights * x       # weighted values 
    N <- sum(weights)     # Estimated population size 
    cw <- cumsum(weights)   # cumulative sum of weights 
    T <- sum(wx)             # Estimated total income 

    dt <- data.table(x = x, weights = weights, wx = wx, key = "x")
    weights0 <- wx0 <- NULL
    dt1 <- dt[, list(weights0 = sum(weights),
                wx0 = sum(wx)), keyby = x][,
                Nk := cumsum(weights0)][,
                wx1 := cumsum(wx0)]
    dt <- merge(dt, dt1)

    # Nk - estimation of the cumulative distribution function 
    Nk <- dt[["Nk"]]

    # wx - weighted partial sum
    wx1 <- dt[["wx1"]]
    dt <- dt1 <- NULL

    Num_eu <- 2 * sum(wx*cw) - sum(weights^2 * x)
    Num <- 2 * sum(wx*Nk) - sum(weights^2 * x)

    Den <- N*T;

    Gini_eu <- 100*(Num_eu/Den-1)
    Gini <- Num/Den-1

    Gini_pr <- 100*Gini
        

    # LINEARIZED VARIABLE OF THE GINI COEFFICIENT (IN %)

    lin <- 100*(2*Nk*(x-wx1/Nk)+T-N*x-Gini*(T+N*x))/(N*T)
    
    if (is.nan(Gini))  Gini_pr <- lin <- 0

    Gini_pr <- data.table(Gini=Gini_pr, Gini_eu=Gini_eu)
    lin_id <- data.table(ids, lin=lin)
    
    return(list(Gini=Gini_pr, lin=lin_id))
}

