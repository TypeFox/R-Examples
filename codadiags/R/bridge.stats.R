E.bridge <- function(bridge) {
    return(abs(mean(bridge)))
}

E.studentbridge <- function(studentbridge) {
    return(E.bridge(studentbridge))
}

E.brownianbridge <- function(brownianbridge) {
    return(E.bridge(brownianbridge))
}

E.loglikbridge <- function(loglikbridge) {
    return(E.bridge(loglikbridge))
}

var.bridge <- function(bridge) {
    return(var(bridge))
}

var.studentbridge <- function(studentbridge) {
    return(var.bridge(studentbridge))
}

var.brownianbridge <- function(brownianbridge) {
    return(var.bridge(brownianbridge))
}

var.loglikbridge <- function(loglikbridge) {
    return(var.bridge(loglikbridge))
}

autocov.bridge <- function(bridge) {
    Nb = length(bridge)
    return(cov(bridge[1:(Nb-1)],bridge[2:Nb]))
}

autocov.studentbridge <- function(studentbridge) {
    return(autocov.bridge(studentbridge))
}

autocov.brownianbridge <- function(brownianbridge) {
    return(autocov.bridge(brownianbridge))
}

autocov.loglikbridge <- function(loglikbridge) {
    return(autocov.bridge(loglikbridge))
}

loglik_mean.brownianbridge <- function(brownianbridge) {
    Nbb = length(brownianbridge)
    N = Nbb-1
    seq.n = seq(from=1,to=N-1,by=1)
    return(-mean(brownianbridge[2:N]^2 / (seq.n/N * (1-seq.n/N))))
}

loglik_mean.studentbridge <- function(studentbridge) {
    Nsb = length(studentbridge)
    return(-mean(log(dt(studentbridge,df=Nsb-2))))
}

loglik_extremum.brownianbridge <- function(brownianbridge) {
    Nbb = length(brownianbridge)
    n.ext = which.max(abs(brownianbridge))
    return(brownianbridge[n.ext]^2 / ((n.ext-1)/(Nbb-1) * (1-(n.ext-1)/(Nbb-1))))
}

loglik_extremum.studentbridge <- function(studentbridge) {
    Nsb = length(studentbridge)
    return(-log(dt(max(abs(studentbridge)),df=Nsb-2)))
}

extremum.bridge <- function(bridge) {
    return(max(abs(bridge)))
}

extremum.studentbridge <- function(studentbridge) {
    return(extremum.bridge(studentbridge))
}

extremum.brownianbridge <- function(brownianbridge) {
    return(extremum.bridge(brownianbridge))
}

extremum.loglikbridge <- function(loglikbridge) {
    return(extremum.bridge(loglikbridge))
}

ratio_extremum.bridge <- function(bridge) {
    max.bridge = max(bridge)
    mmin.bridge = -min(bridge)
    return(max(max.bridge/mmin.bridge,mmin.bridge/max.bridge))
}

ratio_extremum.studentbridge <- function(studentbridge) {
    return(ratio_extremum.bridge(studentbridge))
}

ratio_extremum.brownianbridge <- function(brownianbridge) {
    return(ratio_extremum.bridge(brownianbridge))
}

ratio_extremum.loglikbridge <- function(loglikbridge) {
    return(ratio_extremum.bridge(loglikbridge))
}

ratio_loglik_extremum.brownianbridge <- function(brownianbridge) {
    Nbb = length(brownianbridge)
    N = Nbb-1
    n.max.bridge = which.max(brownianbridge)-1
    n.min.bridge = which.min(brownianbridge)-1
    ll.max.bridge = brownianbridge[n.max.bridge+1]^2 / (n.max.bridge/N * (1-n.max.bridge/N))
    ll.min.bridge = brownianbridge[n.min.bridge+1]^2 / (n.min.bridge/N * (1-n.min.bridge/N))
    return(max(ll.max.bridge/ll.min.bridge,ll.min.bridge/ll.max.bridge))
}

# loglik_maximum.brownianbridge <- function(brownianbridge) {
#     Nbb = length(brownianbridge)
#     N = Nbb-1
#     n.max.bridge = which.max(brownianbridge)-1
#     #n.min.bridge = which.min(brownianbridge)-1
#     ll.max.bridge = brownianbridge[n.max.bridge+1]^2 / (n.max.bridge/N * (1-n.max.bridge/N))
#     #ll.min.bridge = brownianbridge[n.min.bridge+1]^2 / (n.min.bridge/N * (1-n.min.bridge/N))
#     #return(max(ll.max.bridge/ll.min.bridge,ll.min.bridge/ll.max.bridge))
#     return(ll.max.bridge)
# }

ratio_loglik_extremum.studentbridge <- function(studentbridge) {
    Nsb = length(studentbridge)
    ll.max.bridge = log(dt(max(studentbridge),df=Nsb-2))
    ll.min.bridge = log(dt(min(studentbridge),df=Nsb-2))
    return(max(ll.max.bridge/ll.min.bridge,ll.min.bridge/ll.max.bridge))
}
