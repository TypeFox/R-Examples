cmx = function (DATA, threshold = 0.5, which.model = 1, na.rm = FALSE) 
{

### check logicals

    if (is.logical(na.rm) == FALSE) {
        stop("'na.rm' must be of logical type")
    }

### check for and deal with NA values:

    if (sum(is.na(DATA)) > 0) {
        if (na.rm == TRUE) {
            NA.rows <- apply(is.na(DATA), 1, sum)
            warning(length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
            DATA <- DATA[NA.rows == 0, ]
        }
        else {
            return(NA)
        }
    }

### Check that which.model is a single integer, and not greater than number of models in DATA

    if (length(which.model) != 1) {
        stop("this function will only work for a single model, 'which.model' must be of length one")
    }
    if (which.model < 1 || round(which.model) != which.model) {
        stop("'which.model' must be a positive integer")
    }
    if (which.model + 2 > ncol(DATA)) {
        stop("'which.model' must not be greater than number of models in 'DATA'")
    }

### Check that 'threshold' is valid

    if (length(threshold) != 1) {
        stop("'threshold' must be a single number between zero and one")
    }
    if (max(threshold) > 1) {
        stop("'threshold' must be a single number between zero and one")
    }
    if (min(threshold) < 0) {
        stop("'threshold' must be a single number between zero and one")
    }

### translate observations from values to presence/absence

    OBS.ind <- DATA[,2]>0 # This will makes negatives = absences

### Pull out data from single model and apply threshold

    if (threshold == 0) {
        PRED.ind = DATA[,which.model+2] >= threshold
    }
    else {
        PRED.ind = DATA[,which.model+2] > threshold
    }

### Calculate confusion matrix

    C = c( sum( PRED.ind & OBS.ind ), sum( !PRED.ind & OBS.ind ),
      sum( PRED.ind & !OBS.ind ), sum( !PRED.ind & !OBS.ind ) )
    C = as.table( matrix(C,nrow=2) )
    dimnames(C) = list(predicted=c(1,0),observed=c(1,0))
    storage.mode(C)="double"
    
    return(C)
}
