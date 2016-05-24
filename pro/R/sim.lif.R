#' @title  Simulate optogenetic stimulation on a leaky-integrate-fire neuron
#'
#' @description
#' Simulate various kinds of neural measures (e.g. membrane potentials and spikes) from a LIF neuron.
#'
#' @export
#' 
#' @param n Number of time bins. The total time is \code{n} times \code{bin}.
#' @param I Input stimulus vector of  length \code{n}.
#' @param C Membrane capacitance of the simulated neuron.
#' @param R  Membrane resistance  of the simulated neuron.
#' @param Vth  Membrane potential threshold for spiking.
#' @param V0  Membrane potential reset value after spiking.
#' @param bin Time length for each time bin. Default 5 millisecond.
#' @param dt  Time length for each simulation step. Default 0.05 millisecond.
#' @return  a \code{list} of simulated neural spikes, optogenetic light flashes,  and simulation parameters.
#' @examples
#' n<- 500
#' set.seed(100)
#' re <- sim.lif(n, rbinom(n, 1, 0.14), 7, 3)
sim.lif <- function(n, I, C, R, Vth=1, V0=0, bin=5, dt=0.05) {
    ## bin is 5ms, dt is 0.05ms
    ## I[t] = 1, meaning from [t-bin, t] is activated
    T <- n*bin
    nstep <- T/dt
    
    
    ## dt is time scale, only matter when plotting
    V <- rep(0, nstep)
    S <- rep(0, nstep)
    
    stime <- NULL
    for (i in 2:nstep) {
        Ibin <- ceiling((i-1)*dt/bin)
        V[i] <- V[i-1]+(I[Ibin]-V[i-1]/R)*dt/C
        ## if spike or not
        if (V[i]>Vth) {
            S[i] <- 1
            V[i] <- V0
        }
    }

    binsize <- bin/dt
    ind <- ceiling(seq(nstep)/binsize)
    sbin <- sapply(1:n, function(x) sum(S[ind==x]))
    stime <- seq(n)[sbin>=1]

    ## time for V on the fractional counts of bins
    ts <-  (1:nstep)*dt/bin
    ## input sequence, the same length as ts
    Is <- rep(I, each = bin/dt) 
    
    return(list(V=V, S=S, stime=stime, sbin=sbin, I=I, bin=bin, dt=dt,
                C=C, R=R, Vth=Vth, V0=V0, ts=ts, Is = Is))
}


