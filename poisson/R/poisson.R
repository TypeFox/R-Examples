# Functions for simulations and calculations involving Poisson processes.

hpp.event.times <- function(rate, num.events, num.sims=1, t0=0) {
  # Simulate homogeneous Poisson process event times.
  #
  # Get the n consecutive event times of an homogeneous poisson process with given rate.
  # Note: the rate parameter is often referred to as lambda.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :num.sims: number of simulated paths to create
  #   :t0: start time
  #
  # Returns: 
  #   a numeric vector of length num.events if num.sims=1
  #   else, a num.events * num.sims matrix [num.events+1 is prepend.zero=T]

  if(num.sims==1) {
    return(t0 + cumsum(rexp(n=num.events, rate=rate)))
  }
  else {
    x = matrix(rexp(n=num.events*num.sims, rate=rate), num.events)
    return(t0 + apply(x, 2, cumsum))
  }
}

hpp.sim <- function(rate, num.events, num.sims=1, t0=0, prepend.t0=T) {
  # Simulate homogeneous Poisson process(es).
  #
  # Get the n consecutive event times of an homogeneous poisson process with given rate.
  # Note: the rate parameter is often referred to as lambda.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :num.sims: number of simulated paths to create
  #   :t0: start time
  #   :prepend.t0: TRUE to include t0 at the start of the process
  #
  # Returns: 
  #   a numeric vector of length num.events if num.sims=1
  #   else, a num.events * num.sims matrix [num.events+1 is prepend.zero=T]
  
  if(num.sims==1) {
    x = hpp.event.times(rate, num.events, num.sims=1, t0=t0)
    if(prepend.t0)
      return(c(t0, x))
    else
      return(x)
  }
  else {
    x = hpp.event.times(rate, num.events, num.sims=num.sims, t0=t0)
    if(prepend.t0)
      return(rbind(rep(t0, num.sims), x))
    else
      return(x)
  }
}

hpp.mean <- function(rate, t0=0, t1=1, num.points=100, maximum=NULL) {
  # Calculate the expected value of an homogeneous Poisson process at points in time.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :t0: start time
  #   :t1: end time
  #   :num.points: number of points between t0 and t1 to use in estimating mean
  #   :maximum: the optional maximum value that the process should take
  
  times = seq(t0, t1, length.out=num.points)
  x = rate*times
  if(is.null(maximum))
    return(x)
  else
    return(pmin(x, maximum))
}

hpp.mean.event.times <- function(rate, num.events) {
  # Calculate the expected event times of an homogeneous Poisson process.
  # 
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: observe mean event times at this many points
  
  return(seq(from=1/rate, length.out=num.events, by=1/rate))
}

hpp.plot <- function(rate, num.events, num.sims=100, t0=0, t1=NULL, 
                     num.points=100, quantiles=c(0.025, 0.975), ...) {
  # Plot num.events simulated homogeneous Poisson processes, plus
  # the mean and quantiles
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :num.sims: number of simulated paths to plot
  #   :t0: plot start time
  #   :t1: plot end time
  #   :num.points: number of points to use in estimating mean and quantile processes
  #   :quantiles: plot these quantile processes
  
  x = hpp.sim(rate, num.events, num.sims)
  if(is.null(t1))
    t1 = 1.1*max(x)
  plotprocesses(x, xlim=c(t0, t1), ...)
  # Expected value process
  x.bar = hpp.mean(rate, t0, t1, num.points, maximum=num.events)
  lines(seq(t0, t1, length.out=num.points), x.bar, col='darkorange1', lwd=2)
  # Quantile processes
  x.q = t(apply(x, 1, function(x) quantile(x, probs=quantiles)))
  plotprocesses(x.q, col='red', lwd=2, lty=3, add=T)
  return(list(x=x, x.bar=x.bar, x.q=x.q))
}


nhpp.event.times <- function(rate, num.events, prob.func, num.sims=1, t0=0) {
  # Simulate non-homogeneous Poisson process event times.
  #
  # Get the n consecutive event times of a non-homogeneous poisson process.
  # Events are simulated using an homogeneous process with rate,
  # and an event at time t is admitted with probability prob.func(t).
  # This method, called 'thinning' by Lewis & Shedler (1978) is described in
  # Simulation of Non-Homogeneous Poisson Processes by Thinning
  # The rate parameter of an homogeneous process is often called lambda.
  #
  # Params
  #   :rate: the rate at which events occur in the equivalent homogeneous 
  #           Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1
  #   :num.sims: number of simulated paths to create
  #   :t0: the reference start time of all events
  #
  # Returns: 
  #   a numeric vector of length num.events if num.sims=1
  #   else, a num.events * num.sims matrix [num.events+1 is prepend.zero=T]
  
  n <- num.events
  if(n<=0) return(c())
  if(num.sims==1) {
    # As a first attempt, take 2n homogeneous events
    times = hpp.event.times(rate, 2*n, t0=t0)
    # and accept each event with time-varying probability, determined by prob.func
    prob.accept = prob.func(times)
    rands = runif(2*n)
    accept = rands < prob.accept
    if(sum(accept) >= n) {
      # We have enough events to accept
      returnable = times[accept][1:n]
      return(returnable)
    }
    else {
      # We need more events
      extra = hpp.event.times(rate, n-sum(accept), num.sims=1, t0=max(times))
      returnable = c(times[accept], extra)
      return(returnable)
    }
  }
  else {
    f = function(x) nhpp.event.times(rate, num.events, prob.func, num.sims=1, t0=t0)
    return(sapply(1:num.sims, f))
  }
}

nhpp.sim <- function(rate, num.events, prob.func, num.sims=1, t0=0, prepend.t0=T) {
  # Simulate non-homogeneous Poisson process(es).
  #
  # Get the n consecutive event times of a non-homogeneous poisson process.
  # Events are simulated using an homogeneous process with rate,
  # and an event at time t is admitted with probability prob.func(t).
  # This method, called 'thinning' by Lewis & Shedler (1978) is described in
  # Simulation of Non-Homogeneous Poisson Processes by Thinning
  # The rate parameter of an homogeneous process is often called lambda.
  #
  # Params
  #   :rate: the rate at which events occur in the equivalent homogeneous 
  #           Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1
  #   :num.sims: number of simulated paths to create
  #   :t0: the reference start time of all events
  #   :prepend.t0: T to include t0 at the start of the process  
  #
  # Returns: 
  #   a numeric vector of length num.events if num.sims=1
  #   else, a num.events * num.sims matrix [num.events+1 is prepend.zero=T]
  
  if(num.sims==1) {
    x = nhpp.event.times(rate, num.events, prob.func=prob.func, num.sims=1, t0=t0)
    if(prepend.t0)
      return(c(t0, x))
    else
      return(x)
  }
  else {
    x = nhpp.event.times(rate, num.events, prob.func=prob.func, num.sims=num.sims, t0=t0)
    if(prepend.t0)
      return(rbind(rep(t0, num.sims), x))
    else
      return(x)
  }
}

# A slower but easier to understand (e.g. non-recursive) version of nhpp.sim
nhpp.sim.slow <- function(rate, num.events, prob.func, num.sims=1, t0=0, prepend.t0=T) {
  # Simulate a non-homogeneous Poisson process.
  #
  # Get the n consecutive event times of a non-homogeneous poisson process.
  # Events are simulated using an homogeneous process with rate,
  # and an event at time t is admitted with probability prob.func(t).
  # This method, called 'thinning' by Lewis & Shedler (1978) is described in
  # Simulation of Non-Homogeneous Poisson Processes by Thinning
  # The rate parameter of an homogeneous process is often called lambda.
  #
  # Params
  #   :rate: the rate at which events occur in the equivalent homogeneous 
  #           Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1
  #   :num.sims: number of simulated paths to create
  #   :t0: the reference start time of all events
  #   :prepend.t0: T to include t0 at the start of the process  
  #
  # Returns: 
  #   a numeric vector of length num.events if num.sims=1
  #   else, a num.events * num.sims matrix
  
  if(num.sims==1) {
    i=0; t=t0
    times = numeric(length=num.events)
    while(i < num.events) {
      tau = rexp(n=1, rate=rate)
      t.new = t + tau
      if(runif(1) < prob.func(t.new)) {
        i = i + 1
        times[i] = t.new
      }
      t = t.new
    }
    if(prepend.t0)
      return(c(t0, times))
    else
      return(times)
  }
  else {
    f = function(x) nhpp.sim.slow(rate, num.events, prob.func)
    return(sapply(1:num.sims, f))
  }
}

nhpp.mean <- function(rate, prob.func, t0=0, t1=1, num.points=100, 
                      maximum=NULL) {
  # Calculate the expected value of a non-homogeneous Poisson process at points in time.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1  
  #   :t0: start time
  #   :t1: end time
  #   :num.points: number of points between t0 and t1 to use in estimating mean
  #   :maximum: the optional maximum value that the process should take
  
  f <- function(x) rate * prob.func(x)
  times = seq(t0, t1, length.out=num.points)
  y = sapply(times, function(x) integrate(f, lower=0, upper=x)$value)
  if(is.null(maximum))
    return(y)
  else
    return(pmin(y, maximum))
}

nhpp.mean.event.times <- function(rate, num.events, prob.func, max.time=1000) {
  # Calculate the expected event times of a non-homogeneous Poisson process.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: calculate the event times for this many events 
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1  
  #   :max.time: the maximum time needs to be specified to act as upper bound in the 
  #               solver routine. I would like to remove the need to set this variable
  #               but for now I use the arbitrarily high value of 1000...
  
  t0=0
  times = c()
  f <- function(x) rate * prob.func(x)
  for(i in 1:num.events) {
    t1 = uniroot(function(x) integrate(f, lower=t0, upper=x)$value - 1, 
                 lower=t0, upper=max.time)$root
    times = c(times, t1)
    t0 = t1
  }
  return(times)
}

nhpp.plot <- function(rate, num.events, prob.func, num.sims=100,
                              t0=0, t1=NULL, num.points=100, 
                              quantiles=c(0.025, 0.975), ...) {
  # Plot num.events simulated homogeneous Poisson process, plus
  # the mean and quantiles
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :prob.func: function that takes time as sole argument and returns value between 0 and 1
  #   :num.sims: number of simulated paths to plot
  #   :t0: plot start time
  #   :t1: plot end time
  #   :num.points: number of points to use in estimating mean and quantile processes
  #   :quantiles: plot these quantile processes
  
  x = nhpp.sim(rate, num.events, prob.func, num.sims)
  if(is.null(t1))
    t1 = 1.1*max(x)
  plotprocesses(x, xlim=c(t0, t1), ...)
  # Expected value process
  x.bar = nhpp.mean(rate, prob.func, t0, t1, num.points, maximum=num.events)
  lines(seq(t0, t1, length.out=num.points), x.bar, col='darkorange1', lwd=2)
  # Quantile processes
  x.q = t(apply(x, 1, function(x) quantile(x, probs=quantiles)))
  plotprocesses(x.q, col='red', lwd=2, lty=3, add=T)
  return(list(x=x, x.bar=x.bar, x.q=x.q))
}



# Inference
nhpp.lik <- function(x, T1, rate, prob.func) {
  # TODO: explain
  lambda.func = function(x) rate * prob.func(x)
  expectation = integrate(lambda.func, lower=0, upper=T1)$value
  return(sum(log(lambda.func(x))) - expectation)
}

hpp.lik <- function(x, T1, rate) {
  # Get the likelihood of a rate parameter at a specific time for observed HPP event times.
  # Params
  #   :x: a vector of HPP event times
  #   :T1: Calculate likelihood at this time
  #   :rate: the putative HPP event rate
  return(nhpp.lik(x, T1, rate, prob.func=function(x) rep(1, length(x))))
}

hpp.mle <- function(x, T1) {
  # Get the maximum-likelihood rate parameter for given HPP event times.
  # Params
  #   :x: a vector of HPP event times
  #   :T1: Calculate MLE at this time
  return(length(x) / T1)
}

nhpp.mle <- function(x, T1, prob.func, max.val) {
  # TODO: explain
  opt.f = function(z) nhpp.lik(x, T1, rate=z, prob.func=prob.func)
  return(optimize(opt.f, interval=c(0, max.val), maximum=T)$maximum)
}



# Plotting
plotprocesses = function(x, y=NULL, xlab='t (years)', ylab='N', type='l', lty=2, col='cadetblue3', xlim=c(0, 1.1*max(x)), 
                          lwd=0.5, add=F, ...) {
  # TODO: explain
  if(is.null(y))
    y = 0:(dim(x)[1]-1)
  matplot(x, y, xlab=xlab, ylab=ylab, type=type, lty=lty, col=col, xlim=xlim, lwd=lwd, add=add, ...)
}


# Object-oriented approach
setClass("PoissonProcessScenario",
         representation(x="matrix", x.bar="numeric", x.bar.index="numeric", x.q="matrix"),
         prototype(x=matrix(), x.bar=numeric(), x.bar.index=numeric(), x.q=matrix())
)

setMethod("plot" , "PoissonProcessScenario",
          function(x, plot.mean=T, plot.quantiles=T, ...) {
            # Plot simulated processes
            plotprocesses(x@x, ...)
            # Plot mean process
            if(plot.mean)
              lines(x@x.bar.index, x@x.bar, col='darkorange1', lwd=2)
            # Plot quantile processes
            if(plot.quantiles)
              plotprocesses(x@x.q, col='red', lwd=2, lty=3, add=T)
          }
)
setMethod ("show", "PoissonProcessScenario",
           function(object) {
             cat(dim(object@x), "\n")
             cat("\n")
           }
)
hpp.scenario <- function(rate, num.events, num.sims=100, t0=0, t1=NULL, 
                         num.points=100, quantiles=c(0.025, 0.975), ...) {
  # Simulate an homogeneous Poisson process scenario, with sample paths, 
  # expected value process, and quantile processes.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :num.sims: number of simulated paths to plot
  #   :t0: start time of processes
  #   :t1: end time of mean process, inferred automiatically when NULL (default)
  #   :num.points: number of points to use in estimating mean process
  #   :quantiles: calculate these quantile processes
  
  # Simulated processes
  x = hpp.sim(rate, num.events, num.sims, t0=t0, ...)
  if(is.null(t1))
    t1 = 1.1*max(x)
  # Expected value process
  x.bar = hpp.mean(rate, t0=t0, t1=t1, num.points=num.points, maximum=num.events)
  # Quantile processes
  x.q = t(apply(x, 1, function(x) quantile(x, probs=quantiles)))
  return(
    new("PoissonProcessScenario", x=x, x.bar=x.bar, 
        x.bar.index=seq(t0, t1, length.out=num.points), x.q=x.q)
  )
}

nhpp.scenario <- function(rate, num.events, prob.func, num.sims=100, t0=0, t1=NULL, 
                         num.points=100, quantiles=c(0.025, 0.975), ...) {
  # Simulate a non-homogeneous Poisson process scenario, with sample paths, 
  # expected value process, and quantile processes.
  #
  # Params
  #   :rate: the rate at which events occur in the Poisson process, aka lambda
  #   :num.events: number of event times to simulate in each process
  #   :num.sims: number of simulated paths to plot
  #   :t0: start time of processes
  #   :t1: end time of mean process, inferred automiatically when NULL (default)
  #   :num.points: number of points to use in estimating mean process
  #   :quantiles: calculate these quantile processes
  
  # Simulated processes
  x = nhpp.sim(rate, num.events, prob.func, num.sims=num.sims, t0=t0, ...)
  if(is.null(t1))
    t1 = 1.1*max(x)
  # Expected value process
  x.bar = nhpp.mean(rate, prob.func, t0=t0, t1=t1, num.points=num.points, maximum=num.events)
  # Quantile processes
  x.q = t(apply(x, 1, function(x) quantile(x, probs=quantiles)))
  return(
    new("PoissonProcessScenario", x=x, x.bar=x.bar, 
        x.bar.index=seq(t0, t1, length.out=num.points), x.q=x.q)
  )
}

