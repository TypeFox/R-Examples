
checker <- function(variables, datasets, varname) {
      vars <- 0
      if (sum(variables %in% names(datasets)) > 0) { vars <- vars + 100 
         } else vars <- vars + 1
      if (is.numeric(variables)) {
         if ((variables <= ncol(datasets))*(0 < variables)) { vars <- vars + 100
            } else vars <- vars + 10 }
      if (is.logical(variables)) {
         if (length(variables) == ncol(datasets)) { vars <- vars + 100
            } else vars <- vars + 20 }
       
      if (vars == 1)  stop(paste(variables,"does not exist in 'dataset'!", sep = " "))
      if (vars == 11) stop(paste("Column",variables,"does not exist in 'dataset'!", sep = " "))
      if (vars == 21) stop(paste0("'",varname,"' logical vector must be the same length as 'dataset' column count!"))

      return(ifelse(vars >= 100, TRUE, FALSE))
 }

incPercentile <- function(Y, weights = NULL, sort = NULL, 
        Dom = NULL, period=NULL, k = c(20, 80), dataset = NULL) {
   
   ## initializations
   if(length(k) == 0 | any(!is.numeric(k) | k < 0 | k > 100)) {
        stop("'k' must be a vector of integers between 0 and 100")
    } else k <- round(k)
   
   if(!is.null(dataset)) {
       dataset <- data.table(dataset)
       if (checker(Y, dataset, "Y")) Y <- dataset[, Y, with=FALSE] 

       if(!is.null(weights)) {
           if (checker(weights, dataset, "weights")) weights <- dataset[, weights, with=FALSE] }
     
       if(!is.null(sort)) {
           if (checker(sort, dataset, "sort")) sort <- dataset[, sort, with=FALSE] }

       if (!is.null(period)) {
            if (min(period %in% names(dataset))!=1) stop("'period' does not exist in 'dataset'!")
            if (min(period %in% names(dataset))==1) {
                                period <- dataset[, period, with=FALSE]
                                period[, (names(period)):=lapply(.SD, as.character)] }}

       if (!is.null(Dom)) {
            if (checker(Dom, dataset, "Dom")) {
                                Dom <- dataset[, Dom, with=FALSE]
                                Dom[, (names(Dom)):=lapply(.SD, as.character)] }}
      }

   # check vectors
   # Y
   Y <- data.frame(Y)
   n <- nrow(Y)
   if (ncol(Y) != 1) stop("'Y' must be a vector or 1 column data.frame, matrix, data.table")
   Y <- Y[, 1]
   if(!is.numeric(Y)) stop("'Y' must be numerical")
   if (any(is.na(Y))) stop("'Y' has unknown values")

   # weights
   weights <- data.table(weights)
   if (nrow(weights) != n) stop("'weights' must be the same length as 'Y'")
   if (ncol(weights) != 1) stop("'weights' must be vector or 1 column data.frame, matrix, data.table")
   weights <- weights[[1]]
   if(!is.numeric(weights)) stop("'weights' must be numerical")
   if (any(is.na(weights))) stop("'weights' has unknown values")

   # sort  
   if(!is.null(sort)) {
         sort <- data.frame(sort)
         if (nrow(sort) != n) stop("'sort' must be the same length as 'Y'")
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

   # Dom
   namesDom <- NULL
   if(!is.null(Dom)) {
             Dom <- data.table(Dom)
             namesDom <- names(Dom)
             if (is.null(names(Dom))) stop("'Dom' must be colnames")
             if (nrow(Dom) != n) stop("'Dom' must be the same length as 'Y'")
             Dom[, (namesDom):=lapply(.SD, as.character)]
       }
    
    if (!is.null(period)) {
       if (!is.null(Dom)) { Dom <- data.table(period, Dom)
        } else Dom <- period } 
 
    # Percentiles by domain (if requested)
    
    N <- NULL
    if(!is.null(Dom)) {
        Dom_app <- do.call("paste", c(as.list(Dom), sep="__"))
        q1 <- lapply(split(Dom[, .I], Dom_app), function(i) {
               Yind <- Y[i]
               weightsind <- weights[i]
               sortind <- sort[i]
               order <- if(is.null(sortind)) order(Yind) else order(Yind, sortind)
               Yind <- Yind[order]
               weightsind <- weightsind[order]  # also works if 'weights' is NULL                               
               percentile <- weightedQuantile(Yind, weightsind, probs=k/100, sorted=FALSE, na.rm=FALSE)               
               q <- data.table(Dom[i][1], t(percentile))})
        q <- rbindlist(q1)
        setnames(q, names(q)[ncol(Dom)+1:length(k)], paste0("x",k))
        if (!is.null(period)&!is.null(namesDom)) {
              q1 <-  q[,.N, keyby=namesDom][,N:=NULL]
              q2 <- q[, .N, by=names(period)][,N:=NULL]
              qrs <- rbindlist(lapply(1:nrow(q2), function(i) {
                                   data.table(q2[i], q1) }))
              qrs[, (c(paste0("x", k))):=0]
              qrs <- rbind(q, qrs)
              q <- qrs[,lapply(.SD, sum), keyby=names(Dom), .SDcols=paste0("x", k)]              
            }
         setkeyv(q, names(Dom))
     } else {  order <- if(is.null(sort)) order(Y) else order(Y, sort)
               Y <- Y[order]
               weights <- weights[order]  # also works if 'weights' is NULL
               percentile <- weightedQuantile(Y, weights, probs=k/100, sorted=TRUE, na.rm=FALSE)
               q <- data.table(t(percentile))
               setnames(q, names(q)[1:length(k)], paste0("x",k))
     }
     ## return results
    return(q)
}

