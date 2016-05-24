MatchCoefDPZ = function(PZ, dt, N, niter = 50000, burn = 0, sigfac = 1, fh = 0.25/dt, k = 0.001, verbose = TRUE){
#  require(signal)

  #### subfunctions:
  PZ2Mod = function(PZ){
    np = length(PZ$poles)
    nz = length(PZ$zeros)
    cpoles = NULL
    czeros = NULL
    zerogroup = NULL
    polegroup = NULL
    poles = NULL
    zeros = NULL
    if(length(poles) > 0){
      poles = PZ$poles[order(abs(PZ$poles))]
    }
    if(length(zeros) > 0){
      zeros = PZ$zeros[order(abs(PZ$zeros))]
    }
    i = 1; while(i <= np){
      if(abs(Im(PZ$poles[i])) < abs(Re(PZ$poles[i] * 1e-6))){ # small numerical errors can cause small imaginary values to appear--make sure they actually are small
        cpoles = c(cpoles, Re(PZ$poles[i]))
        polegroup = c(polegroup, i)
        i = i + 1
      }else if(Re(PZ$poles[i]) == Re(PZ$poles[i+1]) && Im(PZ$poles[i]) == -Im(PZ$poles[i+1])){
        cpoles = c(cpoles, Re(PZ$poles[i]), -abs(Im(PZ$poles[i]))) # make the imaginary part negative so it works with the 'logprior' function in the inversion.
        polegroup = c(polegroup, i, i)
        i = i + 2
      }else{
        stop('All poles must either be real or paired complex conjugates')
      }
    }
    i = 1; while(i <= nz){
      if(PZ$zeros[i] ==  Re(PZ$zeros[i])){
        czeros = c(czeros, Re(PZ$zeros[i]))
        zerogroup = c(zerogroup, i)
        i = i + 1
      }else if(Re(PZ$zeros[i]) == Re(PZ$zeros[i+1]) && Im(PZ$zeros[i]) == -Im(PZ$zeros[i+1])){
        czeros = c(czeros, Re(PZ$zeros[i]), -abs(Im(PZ$zeros[i])))
        zerogroup = c(zerogroup, i, i)
        i = i + 2
      }else{
        stop('All zeros must either be real or paired complex conjugates')
      }
    }
    return(list(v = c(czeros, cpoles), np = np, nz = nz, group = c(zerogroup, polegroup)))
  }

##################
  Mod2PZ = function(mod, group, np, nz){
    poles = NULL
    zeros = NULL
    i = 1; while(i <= nz){
    # real root:
      if(i == nz || group[i] != group[i+1]){
        zeros = c(zeros, mod[i])
        i = i + 1
      }else{ # complex root
        zeros = c(zeros, mod[i] + 1i * mod[i+1], mod[i] - 1i * mod[i+1])
        i = i + 2
      }
    }
    
    i = 1 + nz; while(i <= (nz + np)){
    # real root:
      if(i == (np + nz) || group[i] != group[i+1]){
        poles = c(poles, mod[i])
        i = i + 1
      }else{ # complex root
        poles = c(poles, mod[i] + 1i * mod[i+1], mod[i] - 1i * mod[i+1])
        i = i + 2
      }
    }
    return(list(poles = poles, zeros = zeros, np = length(poles), nz = length(zeros)))
  }
##### done with subfunctions
  
  f = (1:N - 1)/(N*dt); f = f - 1/dt * (f > 1/(2*dt))
  trueresp = abs(PZ2Resp(PZ, f, FALSE))
  M = PZ2Mod(PZ)
  guess = M$v
  group = M$group
  np = M$np
  nz = M$nz

  sigma = abs(guess)/sigfac

  # find the index of the normalization frequency (10 * low corner)
  fnorm = CalcCorners(PZ, PLOT = FALSE)[1] * 10
  ampnorm = 1
  if(length(fnorm) == 0){ # just in case there are no 6dB corners--use the peak instead
    fnorm = f[which.max(abs(PZ2Resp(PZ, f, PLOT = FALSE)))]
    ampnorm = max(abs(PZ2Resp(PZ, f, PLOT = FALSE)))
  }
  if(fnorm > max(f)/2){
    fnorm = max(f)/10
  }

  wnorm = 2
  
  logprior = function(x){
    # prior knowledge: any stable response must have poles with negative real parts.
    # some elements of x correspond to imaginary parts, but those are +/- pairs, so force them to be negative.
    if(any(x > 0)){
      -Inf
    }else{
      0
    }
  }
    
  
  CalcMisfit = function(mod, trueresp, np, nz, group, N, dt, freq, wnorm, ampnorm, k, optgain = TRUE){
    trueresp = trueresp/abs(trueresp[wnorm]) # normalize trueresp with respect to wnorm
    # convert mod to poles/zeros
    DPZ = Mod2PZ(mod, group, np, nz)
    coef = PZ2Coef(DPZ, dt)
    a = coef$a
    b = coef$b
    # inputs coefficients of difference equation; returns misfit
    f = (1:N - 1)/(N*dt); f = f - 1/dt * (f > 1/(2*dt))
    impulse = 1:N == 1
    testresp = abs(fft(filter(b, a, impulse)))
    w1 = which(f > 0 & f <= (fh))
    w2 = which(f > 0 & f > fh)  
    if(optgain){
        testresp = testresp * sum(abs(trueresp[w1] * testresp[w1]))/sum(abs(testresp[w1])^2) # this normalization minimizes rms misfit between normalized testresp and trueresp
    }else{
        testresp = testresp/abs(testresp[wnorm]) * ampnorm # normalize testresp with respect to wnorm
    }
    out = -0.5 * sum((log(testresp[w1]) - log(trueresp[w1]))^2/f[w1]) # divide by f[w] to make sure that all frequencies are weighted equally over a log scale--so each decade contributes the same amount
    out = out + k * (-0.5 * sum((log(testresp[w2]) - log(trueresp[w2]))^2/f[w2])) # this term is to try to match high frequencies with secondary importance--k should be low to stop misfit from accumulating for f < fh
    if(is.na(out)){
      out = -Inf
    }
    out
  }

  loglikelihood = function(coef, ...){
    CalcMisfit(coef, ...)
  }

  loglikelihood = CalcMisfit
  
  generate = function(x, sigma){
    w = ceiling(runif(1) * length(x))
    x[w] = x[w] + rnorm(1, 0, sigma[w])
    return(x)
  }
  
  logproposal = function(x1, x2, sigma){
    -0.5 * sum(((x1) - (x2))^2/(sigma+1e-12)^2)
  }
  
  inv = Metropolis(loglikelihood, sigma, guess, niter, generate, logproposal, logprior, burn, save_int = 1, trueresp = trueresp, np = np, nz = nz, group = group, N = N, dt = dt, wnorm = 2, ampnorm = ampnorm, k = k, verbose = verbose)

  DPZ = Mod2PZ(inv$best$mod, group, np, nz)
  coef = PZ2Coef(DPZ, dt)
  DPZ$Knorm = 1
  a = coef$a
  b = coef$b
  f = (1:N - 1)/(N*dt); f = f - 1/dt * (f > 1/(2*dt))
  impulse = 1:N == 1
  testresp = abs(fft(filter(b, a, impulse)))
  w = which(f > 0 & f < fh)
  DPZ$Knorm = sum(abs(trueresp[w] * testresp[w]))/sum(abs(testresp[w]^2)) * DPZ$Knorm/PZ$Sense

  DPZ$Sense = PZ$Sense
  DPZ$dt = dt
  coef = PZ2Coef(DPZ, dt)
  a = coef$a
  b = coef$b

  DPZ = PZ
  DPZ$dt = dt
  DPZ$Zpg = Zpg(pole = 1/polyroot(a), zero = 1/polyroot(b), gain = PZ$Sense)

  impulse = 1:N == 1
  testresp = abs(fft(filter(DPZ$Zpg, impulse)))
  w = which(f > 0 & f < fh)
#  browser()
#  DPZ$Zpg$gain = DPZ$Zpg$gain * sum(abs(log(abs(trueresp[w])) * log(abs(testresp[w]))))/sum(abs(log(abs(testresp[w])^2)))
  DPZ$Zpg$gain = DPZ$Zpg$gain * exp(mean(log(trueresp[w]/testresp[w])))
  
  f = (1:N - 1)/(N*dt); f = f - 1/dt * (f > 1/(2*dt))
  impulse = 1:N == 1
  testresp = abs(fft(filter(DPZ$Zpg, impulse)))
  DPZ$fmax = f[which(exp(abs(log(abs(testresp/trueresp)))) > 1.01)[2]] # [2] instead of [1] because 1 corresponds to f = 0, which typically has infinite misfit
#  browser()

  w = which(f > 0 & f < 0.5/dt)
  plot(f[w], trueresp[w], type = 'l', log = 'xy')
  lines(f[w], testresp[w], col = 'red')

  w = which(f > 0 & f < fh)
  grms = exp(mean(log(testresp[w]/trueresp[w])^2 / f[w]) / mean(1/f[w]))
  abline(v = range(f[w]))
  return(list(b = b, a = a, analogresp = trueresp, digitalresp = testresp, inv = inv, error = grms, DPZ = DPZ))

}
