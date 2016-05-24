## Apply a collection of statistics with only one call
tdStats <- function(m, o,
                   functions = c('mo', 'mm', 'sdo', 'sdm',
                       'mbe', 'mae', 'rmse',
                       'nmbe', 'cvmbe',
                       'nmae', 'cvmae',
                       'nrmse', 'cvrmse',
                       'r2','tStone')){
    mo <- function(m, o) mean(o, na.rm = TRUE)
    mm <- function(m, o) mean(m, na.rm = TRUE)

    sdo <- function(m, o) sd(o, na.rm=TRUE)
    sdm <- function(m, o) sd(m, na.rm=TRUE)
    
    mbe <- function(m, o) mean(m - o, na.rm=TRUE)
    mae <- function(m, o) mean(abs(m - o), na.rm=TRUE)
    rmse <- function(m, o) sqrt(mean((m - o)^2, na.rm=TRUE))
    
    nmbe <- function(m, o) mbe(m, o) / diff(range(o, na.rm = TRUE))
    cvmbe <- function(m, o) mbe(m, o) / mean(o, na.rm=TRUE)
    
    nmae <- function(m, o) mae(m, o) / diff(range(o, na.rm = TRUE))
    cvmae <- function(m, o) mae(m, o) / mean(o, na.rm=TRUE)
    
    nrmse <- function(m, o) rmse(m, o) / diff(range(o, na.rm = TRUE))
    cvrmse <- function(m, o) rmse(m, o) / mean(o, na.rm=TRUE)
    
    
    ## Stone1993
    tStone <- function(m, o) {
        N <- NROW(m)
        MBE <- mbe(m, o)
        RMSE <- rmse(m, o)
        sqrt((N-1) * MBE^2 /(RMSE^2 - MBE^2))
    }
    
    r2 <- function(m, o) cor(m, o)^2
    
    ss <- sapply(functions,
                 FUN=function(f) do.call(f, list(m, o)))
    as.data.frame(t(ss))
}


applyStats <- function(models, o,
                       functions = c('mo', 'mm', 'sdo', 'sdm',
                           'mbe', 'mae', 'rmse',
                           'nmbe', 'cvmbe',
                           'nmae', 'cvmae',
                           'nrmse', 'cvrmse',
                           'r2','tStone')){
    nModels <- ncol(models) 
    nms <- names(models)
    
    errModel <- lapply(seq_len(nModels),
                       FUN = function(i){
                           err <- tdStats(models[,i], o)
                           err$model <- nms[i]
                           err
                       })
    
    do.call(rbind, errModel)
}
