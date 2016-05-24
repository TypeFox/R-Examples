gv.model <- function(gv, model='spherical', sill=NA, 
                     range=NA, nugget=0, slope=NA, ctrl=nls.control()) {

    dt <- data.frame(x=gv$lag[!is.na(gv$lag)], 
                     y=gv$gamma[!is.na(gv$gamma)])

    models <- c("gaussian", "exponential", "spherical", "linear")
    model  <- agrep(model, models)

    # add additional info to object gv
    if(!mtest.gv(gv)) gv$model <- list()
    gv$model$type <- model

    fit <- NULL

    if (model %in% 1:3) {
        # Prepare start list for nls
        start <- list()
        if (is.na(sill)) 
            start[['S']] <- max(dt$y) else S <- sill
        if (is.na(range))
            start[['R']] <- min(dt$x[dt$y == max(dt$y)]) else R <- range
        if (is.na(nugget))
            start[['N']] <- 0 else N <- nugget

        if (length(start) > 0) {
            if (model == 1) {
                # Fit Gaussian model
                fit <- nls(y ~ N + (S-N) * (1-exp(-3*(x**2)/(R**2))), 
                           data=dt, control=ctrl, 
                           start=start)
            } else if (model == 2) {
                # Fit exponential model (decay increasing form)
                fit <- nls(y ~ N + (S-N)*(1-exp(-x/(R))), 
                           data=dt, control=ctrl, 
                           start=start)
            } else if (model == 3) {
                # Fit spherical
                fit <- nls(y ~ N+(S-N)*(((3*(x))/(2*R))-((x)**3/(2*R**3))), 
                           data=dt, control=ctrl, 
                           start=start)
            }
        }

        cf <- coef(fit)
        if ('S' %in% names(cf)) 
            gv$model$sill <- cf[['S']] else gv$model$sill <- S
        if ('N' %in% names(cf)) 
            gv$model$nugget <- cf[['N']] else gv$model$nugget <- N
        if ('R' %in% names(cf)) 
            gv$model$range <- cf[['R']] else gv$model$range <- R

    } else {
        start = list()
        if (is.na(nugget)) 
            start[['N']] <- 0 else N <- nugget
        if (is.na(slope))  
            start[['Sl']] <- 1 else Sl <- slope
        if (is.na(range)) {
            ndt <- dt
            if (is.na(sill)) S <- max(dt$y) else S <- sill
        } else {
            ndt <- dt[dt$x<range,]
            if (is.na(sill)) S <- mean(dt$y[dt$x >= range]) else S <- sill
        }
 
        if (length(start) > 0) {
            # fit linear with points < range
            fit <- nls(y ~ N + Sl * x,
                       data=ndt, control=ctrl,
                       start=start)
        }

        gv$model$sill <- S
        cf <- coef(fit)
        if ('N' %in% names(cf)) 
            gv$model$nugget <- cf[['N']] else gv$model$nugget <- N
        if ('Sl' %in% names(cf)) 
            gv$model$slope <- cf[['Sl']] else gv$model$slope <- Sl
        # estimate range by the model
        gv$model$range <- (gv$model$sill - gv$model$nugget) / gv$model$slope
    }

    class(gv) <- 'gv'

    # add residuals to gv
    pred <- predict(gv)
    residuals <- dt$y - pred
    gv$residuals <- residuals

    gv
}

