fill.demogdata <-
function(data, series = names(data$rate)[1],
                           method=c('perks', 'interpolate', 'mspline'), ...){
    method <- match.arg(method)
    if (class(data) != "demogdata" || data$type != "mortality")
        stop("Not mortality data")
    mx <- data$rate[[series]]
    m.ind <- bool(abs(mx) > 1e-09, na=F)
    warn <- function(nr) {
        if (nr) warning(paste(' A total of', nr, 
                              '0/NA central mortality rates are re-estimated by the',
                              mark(method, view=F), 'method.'), call.=F)
    }
    method.lab <- NULL
    if (method == 'perks'){
        method.lab <- c(method.lab, method)
        # initial number of 0/NA indexes (m.ind is updated after each perks fitting)
        tmp <- sum(!m.ind)
        # determine the years containing 0 or NA rates:
        y.ind <- not(apply(m.ind, 2, all))
        for(y in as.seq(y.ind)){
            x.ind <- m.ind[,y]
            x <- data$age[x.ind]
            m <- mx[x.ind,y]
            # generalized perks model fit (non-linear least square):
            fit <- try(nls(m ~ perks.d(x, a, b, p, 40),  start = c('a'=0.7, 'b'=5.5, 'p'=0.1), ...), silent=T)
            if (class(fit) != 'nls'){
                cat(' - Generalized Perks failed for year ', names(y.ind)[y], ':\n ', fit)
                cat(' Trying simple Perks model ... ')
                fit <- try(nls(m ~ perks.d(x, a=1, b, p, 40),  start = c('b'=6.5, 'p'=0.11),
                               control=list(maxiter=100), ...), silent=T)
            }
            if (class(fit) == 'nls'){
                cat('year', names(y.ind)[y], 'fitted OK. Re-estimated:', sum(!x.ind), 'data cells.\n')
                mx[!x.ind,y] <- predict(fit, list(x=data$age[!x.ind]))
                # reset the index in case some years fail to fit
                m.ind[,y] <- T
            }
            else{
                y.ind[y] <- 'fail'
                cat('\t this failed too:\n', fit)
            }
        }
        warn(tmp - sum(!m.ind)) # total number of Perks corrections
        if (sum(!m.ind)){ 
            method <- 'interpolate'
            y.ind <- ifelse(y.ind == 'fail', T, F)
            cat('\n Failed years:', names(y.ind)[y.ind],
                '\n are interpolated instead.\n', fill=T)
        } 
    }
    if (method == 'mspline'){
        method.lab <- c(method.lab, method)
        sm.data <- smooth.demogdata(data)
        mx[!m.ind] <- sm.data$rate[[series]][!m.ind]
    }
    if (method == 'interpolate'){
        method.lab <- c(method.lab, method)
        mx[!m.ind] <- 0
        x.ind <- not(apply(m.ind, 1, all))
        for(x in as.seq(x.ind))
            mx[x,] <- fill.zero(mx[x,])
    }
    data$rate[[series]] <- mx
    method.lab <- paste(method.lab, collapse=' + ')
    data$label <- paste(data$label, ' (Corrected: ', method.lab, ')', sep='')
    warn(sum(!m.ind)) # total number of mspline or interpolate corrections
    data
}
