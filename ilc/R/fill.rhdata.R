fill.rhdata <-
function(data, method=c('mspline', 'interpolate', 'perks'), ...){
    # check data type:
    if (class(data) != "rhdata") stop("Not \"rhdata\" class mortality data object!")
    # rearrange data into a 2 dimensional format:
    mu <- matrix(data$mu, length(data$age), byrow=F)      
    pop <- matrix(data$pop, length(data$age), byrow=F)
    dims <- dim(data$mu); dimns <- dimnames(data$mu)
    ddata <- demogdata(mu, pop, data$age, seq(prod(dims[-1])), 'mortality', 
                       data$label, data$name)
    method <- match.arg(method)
    fdata <- try(fill.demogdata(ddata, method=method, ...), silent=TRUE)
    if (class(fdata) != 'demogdata'){
        if (method != 'interpolate') { # retry with interpolation
            fdata <- try(fill.demogdata(ddata, method='interpolate', ...), silent=T)
            if (class(fdata) != 'demogdata'){ 
                warning(paste('Both', mark(method,F), 'and \"interpolate\" methods failed.'))
                fdata <- ddata
            }
            else warning(paste('Method', mark(method,F), 'had to be replaced by \"interpolate\".')) 
        }
    } 
    data$mu <- array(fdata$rate[[1]], dim=dims, dimnames=dimns)
    # data$deaths <- data$mu*data$pop
    data
}
