getint <- function(rescount, N, cp, epsilon){
  cupraw <- getcup(N, N-rescount[1], (1-cp)/2)
  clowraw <- getclow(N, rescount[2], (1-cp)/2)
  cup <- cupraw/(1-epsilon)
  clow <- (clowraw-epsilon)/(1-epsilon)
  c(max(0,clow),min(1,cup))
}

getcup <- function(n,x,gamma){
  if (x==n)  1
  else
    uniroot(function(cup)pbinom(x,n,cup)-gamma, c(0,1), tol = .Machine$double.eps^0.25/10000)$root
}

getclow <- function(n,x,gamma){
  if (x==0)  0
  else
    uniroot(function(clow)1-pbinom(x-1,n,clow)-gamma, c(0,1), tol = .Machine$double.eps^0.25/10000)$root
}

maxdelta = function(betapilot, N, cp, epsilon){
  if (N <= 2){stop("N must be larger or equal to 3"); return(NA)}
  lastnones = N-2; lastlastnones = 0
  while (abs(lastnones - lastlastnones)>1.5){
    currentnones = floor((lastlastnones + lastnones)/2)
    betaplus = getint(c(N - 2 - currentnones, currentnones),N, cp,epsilon)[2]
    if (betaplus <= betapilot){
      lastlastnones = max(lastnones, lastlastnones)
    }
    else{
      lastlastnones = min(lastnones, lastlastnones)
    }
    lastnones = currentnones
  }
  delta1 = abs(diff(pmin(getint(c(N - 2 - lastnones, lastnones),N,cp,epsilon), betapilot)))
  delta2 = abs(diff(pmin(getint(c(N - 2 - lastlastnones, lastlastnones),N,cp,epsilon), betapilot)))
  max(delta1, delta2)
}


estop <-function(N, beta, delta, cp, epsilon){
  nzerosmax <- floor(N*(1-beta)); nonesmax = ceiling(N*beta);
  lenmin <- abs(diff(getint(rescount = c(nzerosmax, nonesmax), N=N, cp=cp, epsilon=epsilon)))
  if (lenmin > delta){return(NA)}
  lenmax <- abs(diff(getint(rescount = c(0, 1), N=N, cp=cp, epsilon=epsilon)))
  if (lenmax <= delta){return(1)}
  
  lastlastt <- 1;
  lastt <- N
  while (abs(lastt - lastlastt)>1.5){
    currentt <- floor((lastlastt + lastt)/2)
    nzeros <- floor(currentt*(1-beta))
    nones <- ceiling(currentt*beta)
    len <- abs(diff(getint(c(nzeros, nones), N=N, cp=cp, epsilon=epsilon)))
    if (len <= delta){
      lastlastt <- min(lastt, lastlastt)
    }
    else{
      lastlastt <- max(lastt, lastlastt)
    }
    lastt <- currentt
  }
  max(lastt, lastlastt)
}


getFakt <- function(nsteps,alg){
  Lakt <-getL(alg, nsteps)
  Uakt <-getU(alg, nsteps)
  porig <- (alg@internal$porig)
  stepfun((Lakt+1):(Uakt-1),c(0,cumsum(porig)/sum(porig)))
}

getprobs <- function(Ss, eta, xi, nsteps, alg){
  n = length(Ss)
  pvalones <- sapply(1:n, function(i){
    minones = i; minzeros = i;
    Fakt = getFakt(nsteps, alg)
    n = length(Ss)
    Ssort = sort(Ss)
    Smax = Ssort[1:(n-minones+1)]
    Tminones = sum(Fakt(Smax) >= 1-eta)
    pvalones = 1-pbinom(Tminones -1, n, eta)
    Smin = Ssort[minzeros:n]
    Tminzeros = sum(Fakt(Smin) <= eta)
    pvalzeros = 1-pbinom(Tminzeros -1, n, eta)
    c(pvalones, pvalzeros)
  }

                     )
}




testhyp <- function(rescount, rescounthyp, Ss, eta, xi, nsteps, alg, reports){
  n <- rescounthyp[2]-rescount[2] ## hoping to find at least n ones
  m <- rescounthyp[1]-rescount[1] ## hoping to find at least m ones

  r <- length(Ss); 
  Ssort <- sort(Ss)
  Fakt <- getFakt(nsteps, alg)

  if (n==0){pvalones <- -1}
  else{
    Tforones <- sum(Fakt(Ssort[n:r])<= eta)
    pvalones <- 1-pbinom(Tforones -1, r-n+1, eta)
  }

  if (m==0){pvalzeros <- -1}
  else{
    Tforzeros <- sum(Fakt(Ssort[1:(r-m+1)]) >= 1- eta)
    pvalzeros <- 1-pbinom(Tforzeros -1, r-m+1, eta)
  }
  ##if (reports){print(paste("p-value H+:", pvalones, "p-value H-:", pvalzeros, "xi/2:", xi/2))}
  if (pvalones < xi/2 & pvalzeros < xi/2){
  ##  if (reports){print("Test accepted")}
    1
  }
  else {0}
}


setClass("mcpres", representation(N="numeric", effort="numeric", beta="numeric", rescount="numeric", int="numeric", cp="numeric", taccepted="logical", truncated="logical"))

setMethod("show", signature(object="mcpres"),
          function(object){
            int1 <- floor(object@int[1]*1e4)/1e4
            int2 <- ceiling(object@int[2]*1e4)/1e4
            beta <- min(max(round(object@beta, digits=4), int1), int2)
            cat("Confidence interval (coverage prob: ", object@cp, "): [", int1, ",", int2, "]\n", sep="")
            cat("Estimate of power: ", beta, "\n", sep="")
          }
          )

mkmcpresobj<-function(res){
  r = new("mcpres", N=res$N, effort=res$effort, beta=res$beta, rescount=res$rescount, int=res$int, cp=res$cp, taccepted=res$taccepted, truncated=res$truncated)
}


getpower <- function(genstream, alpha, delta, getintcallwtest, getintcallwotest, N, alg, pilot, options){
  tryCatch({
    undecided <- list()
    rescount <- c(0,0)
    beta <- 1
    effort <- pilot$effort
    int <- pilot$int
    taccepted <- FALSE
    truncated <- FALSE
    finished <- FALSE
    
    makeres <- function(){
      list(N=N, effort=effort, beta=beta, rescount=rescount, int=int, cp=options$cp, taccepted=taccepted, truncated=truncated)
    }

    calcprog <- function(){
      progress <- round(effort/min(pilot$esteffort, options$maxeffort), digits = 2)*100
    }
    
    catprog <- function(){
      if (options$reports){
        nbck <- nchar(options$mess)
        cat(rep("\b", nbck), sep="")
        digits <- nchar(progress)
        if (finished){
          if (resmainloop$taccepted){
            options$mess <<- paste(progress,"%", paste(rep(" ", 18-digits), collapse=""),"[", round(resmainloop$int[1], digits=3),",", round(resmainloop$int[2], digits=3), "]*\n", sep="")
          }
          else {
            options$mess <<- paste(progress,"%", paste(rep(" ", 18-digits), collapse=""),"[", round(resmainloop$int[1], digits=3),",", round(resmainloop$int[2], digits=3), "]\n", sep="")
          }
        }
        else {
          options$mess <<- paste(progress,"%", paste(rep(" ", 18-digits), collapse=""),"[", round(resmainloop$int[1], digits=3),",", round(resmainloop$int[2], digits=3), "]", sep="")
        }
        cat(options$mess)
      }    
    }


    if(!is.null(options$file)){
      cat("CIlow", "CIhigh", "power", "positives", "negatives", "remaining", "effort", sep=",", file=options$file)
      cat("\n", append=TRUE, file=options$file)
    }
    catfile <- function(){
      if (!is.null(options$file)){
        cat(resmainloop$int[1], resmainloop$int[2], resmainloop$beta, resmainloop$rescount[1], resmainloop$rescount[2], resmainloop$N-(resmainloop$rescount[1]+resmainloop$rescount[2]), resmainloop$effort, sep=",", file=options$file, append=TRUE); cat("\n", file=options$file, append=TRUE)
      }
    }

    
    resmainloop <- makeres()
    progress <- calcprog()
    catprog()
    catfile()

    nbatches <- 1
    inc <- options$maxstepsinc
    maxbatch <- options$maxbatch
    maxsteps <- min(options$maxstepsbase, maxbatch)
    nsteps <- maxsteps


    ##start-up
    genericres <- run(alg, function(){},0)
    L <-getL(alg, nsteps);
    replicate(N,{
      if (!finished){
      beta <<- beta - 1/N
      res <- genericres
      res@gen <- genstream()
      res <- cont(res,maxsteps)
      effort <<- effort+res@steps
      if(is.na(res@p.value)) {
        undecided[[length(undecided)+1]] <<- res
        beta <<- beta + (res@pos <= res@steps*alpha)/N
      }
      else {
        if (res@p.value<=alpha){
          rescount[2] <<- rescount[2]+1; beta <<- beta + 1/N
        }
        else{
          rescount[1] <<- rescount[1]+1; 
        }
        int <<- getintcallwotest(rescount)
        len <- diff(int)
        mid <- (int[1]+int[2])/2
        deltareq <- options$deltamid(mid)
        resmainloop <<- makeres()
        progressnew <- calcprog()
        if (progressnew > progress){
          progress <<- progressnew
          catprog()
        }
        if (len <= deltareq){
          finished <<- TRUE
        }
      }
      if (effort > options$maxeffort){
        finished <<-TRUE
        truncated <<-TRUE
      }
    }})
    
    if (finished){cat("\n");return(mkmcpresobj(resmainloop))}

    if (nbatches%%options$batchesbetweentests == 0){
      Ss = c()
      for (i in undecided){
        Ss = append(Ss, i@pos)
      }
      int <- getintcallwtest(rescount, Ss, nsteps, alg)
      len <- diff(int)
      mid <- (int[1]+int[2])/2
      deltareq <- options$deltamid(mid)

      if (len <= deltareq){
        taccepted <- TRUE
        finished <- TRUE
        resmainloop <- makeres()
        catprog()
        catfile()
        return(mkmcpresobj(resmainloop))
      }
    }
    
    resmainloop <- makeres()
    progress <- calcprog()
    catprog()
    catfile()


    ## Then:
    while (len > deltareq){
      nbatches <- nbatches + 1
      maxsteps <- min(ceiling(inc*maxsteps), maxbatch)
      nsteps <- nsteps+maxsteps

      pos <- 1
      for (j in 1:length(undecided)){
        i <- undecided[[pos]]
        effort <- effort - i@steps
        beta <- beta - (i@pos <= i@steps*alpha)/N
        L <-getL(alg, nsteps)        
        aktpath <- cont(i,maxsteps)
        effort <- effort + aktpath@steps
        
        
        if (is.na(aktpath@p.value)){
          beta <- beta + (aktpath@pos <= aktpath@steps*alpha)/N
          undecided[[pos]] <- aktpath
          pos <- pos+1
        }
        else{
          undecided[[pos]] <-c()
          ##don't update pos
          if (aktpath@p.value<=alpha){
            rescount[2] <- rescount[2]+1; beta <- beta + 1/N
          }
          else{
            rescount[1] <- rescount[1]+1; 
          }
          int <- getintcallwotest(rescount)
          len <- diff(int)
          mid <- (int[1]+int[2])/2
          deltareq <- options$deltamid(mid)
          resmainloop <- makeres()
          progressnew <- calcprog()
          if (progressnew > progress){
            progress <- progressnew
            catprog()
          }


          if (len <= deltareq){
            finished <- TRUE;
            break
          }
        }
        if (effort > options$maxeffort){
          finished <-TRUE
          truncated <-TRUE
          break
        }
      }      
      resmainloop <- makeres()

      catprog()
      catfile()


      if (finished){cat("\n");return(mkmcpresobj(resmainloop))}

      if (nbatches%%options$batchesbetweentests == 0){
        Ss = c()
        for (i in undecided){
          Ss = append(Ss, i@pos)
        }
        int <- getintcallwtest(rescount, Ss, nsteps, alg)
        len <- diff(int)
        mid <- (int[1]+int[2])/2
        deltareq <- options$deltamid(mid)

        if (len <= deltareq){
          taccepted <- TRUE
          finished <- TRUE
          resmainloop <- makeres()
          catprog()
          catfile()
          return(mkmcpresobj(resmainloop))
        }
      }
    }

  }, interrupt = function(cond){resmainloop$truncated <- TRUE; mkmcpresobj(resmainloop)})}



getminN = function(delta, cp, epsilon, minN = 3, maxN = 1e9){
  if (epsilon/(1-epsilon) >= delta){stop("epsilon too large vs Delta")}
  nzerosmax = floor((maxN-2)/2); nonesmax = ceiling((maxN-2)/2)
  len <- abs(diff(getint(rescount = c(nzerosmax, nonesmax), N=maxN, cp=cp, epsilon=epsilon)))
  if (len > delta){stop("maxN not large enough")}
  nzerosmin = floor((minN-2)/2); nonesmin = ceiling((minN-2)/2)
  len <- abs(diff(getint(rescount = c(nzerosmin, nonesmin), N=minN, cp=cp, epsilon=epsilon)))
  if (len <= delta){return(minN)}
  
  lastN = maxN; lastlastN = minN
  while (abs(lastN - lastlastN)>1.5){
    currentN = floor((lastlastN + lastN)/2)
    len = abs(diff(getint(rescount = c(floor((currentN - 2)/2), ceiling((currentN - 2)/2)), N = currentN, cp = cp, epsilon = epsilon)))
    if (len <= delta){
      lastlastN = min(lastN, lastlastN)
    }
    else{
      lastlastN = max(lastN, lastlastN)
    }
    lastN = currentN
  }
  max(lastN, lastlastN)
}

getlowpowerminN = function(delta, cp, betapilot, epsilon, minN = 3, maxN = 1e9){
                                        #remember that (1-cp) for betapilot + (1-cp) for main run must = 1-total cp.
                                        #i.e. here use, e.g., cp = 1-(1-cp total cp)/2
  if (betapilot <= delta){print("Nothing to do - we should never get here"); return(0)}
  lenmax = maxdelta(betapilot = betapilot, N = minN, cp = cp, epsilon=epsilon); lenmin = maxdelta(betapilot, N = maxN, cp, epsilon)
  if (lenmax <= delta){return(minN)}
  else if (lenmin > delta){print("Error: maxN is not large enough"); return(NA)}
  lastN = maxN; lastlastN = minN
  while (abs(lastN - lastlastN)>1.5){
    currentN = floor((lastlastN + lastN)/2)
    len = maxdelta(betapilot=betapilot, N = currentN, cp = cp, epsilon=epsilon)
    if (len <= delta){
      lastlastN = min(lastN, lastlastN)
    }
    else{
      lastlastN = max(lastN, lastlastN)
    }
    lastN = currentN
  }
  max(lastN, lastlastN)
}

checkN <- function(N, deltamid, getintcall, tol = 1e-10){
  nzeros = sample(0:(N-2), size = N-1, replace = FALSE)
  for (nzero in nzeros){
    int = getintcall(c(nzero, N-2-nzero), N)
    mid = (int[1] + int[2])/2
    if (diff(int) > deltamid(mid)+tol){
      return(FALSE)
    }
  }
  return(TRUE)
}

getminNfunctional <- function(deltamid, getintcall){
  N <- 3
  while (!checkN(N = N, deltamid = deltamid, getintcall = getintcall)){
    N <- 2*N
    if (N > 1e9){warning("N needs to be larger than 1e9!")}
  }
  if (N == 3) return(3)
  lastN <- N
  lastlastN <- N/2
  while(abs(lastN - lastlastN)>1.5){
    currentN = floor((lastlastN + lastN)/2)
    if (checkN(N = currentN, deltamid = deltamid, getintcall = getintcall)){
      lastlastN <- min(lastN, lastlastN)
    }
    else {
      lastlastN <- max(lastN, lastlastN)
    }
    lastN = currentN    
  }
  max(lastN, lastlastN)
}



getpilot <-function(genresamp, alpha, epsilon, n=1000, maxsteps=1000){
  pilot <- list()
  pilot$n <- n; pilot$maxsteps = maxsteps; pilot$alpha = alpha; pilot$epsilon = epsilon
  alg = getalgprecomp(level = alpha, epsilon = epsilon)
  pilot$res <- replicate(n, {
    r <- run(alg,genresamp(), maxsteps=maxsteps)
    if (is.na(r@p.value)){
      c(as.integer(r@pos <=  maxsteps*alpha), NA)
    }
    else{
      c(as.integer(r@p.value <= alpha), r@steps)
    }
  })
  pilot
}

eeffort <- function(N, delta, cp, pilot, maxsurv = 1e30){
  pearly <- sum(pilot$res[2,]>0, na.rm = TRUE)/pilot$n
  if (pearly == 0){eearly = 0}
  else {eearly <- mean(pilot$res[2,], na.rm = TRUE)} 
  beta <- mean(pilot$res[1,])
  k <- N - estop(N = N, beta = beta, delta = delta, cp = cp, epsilon = pilot$epsilon)
  if (is.na(k)){
                                        #print("N not large enough");
    return(Inf)}
  if (pearly ==1||k> N*(1-pearly)){
    ts <- sort(pilot$res[2,])
    index = (N-k)/N * pilot$n;
    if (ceiling(index)==index){return((N-k)*mean(ts[1:index]) + k*ts[index])}
    else {
      w = ceiling(index) - index
      t = w*ts[floor(index)] + (1-w)*ts[ceiling(index)]
      return((N-k)*mean(c(t, ts) + k*t))
    }
  }
  surv <- function(x){sqrt(log(x)/x) /(sqrt(log(pilot$maxsteps)/pilot$maxsteps))}
  if (k < 2){
                                        #print("Need more than N-2 to terminate, returning Inf");
    return(Inf)}
  if (surv(maxsurv) >= k/(N*(1-pearly))){
                                        #print("Effort high, returning Inf");
    return(Inf)}
  else{
    q <- uniroot(function(x){surv(x)-k/(N*(1-pearly))}, c(pilot$maxsteps, maxsurv))$root
    elate <- pilot$maxsteps + integrate(surv, lower = pilot$maxsteps, upper = q)$value
    result <- N*(eearly * pearly + elate*(1-pearly))
  }
}

mkdeltamid <- function(mindelta=0.02, maxdelta=0.1, llim=0.05, rlim=0.95){
  if (mindelta > maxdelta){stop("mindelta must be smaller than maxdelta")}
  if (llim > rlim){stop("llim must be smaller than rlim")}
  rlim <- min(1,rlim)
  llim <- max(0,llim)
  Vectorize(function(mid){    
    if (mid <= llim){mindelta}
    else if (mid >= rlim){mindelta}
    else {min(max(mindelta, min(2*(mid - llim), 2*(rlim -mid))), maxdelta)}
  })
}

mcp <- function(genstream, alpha=0.05, delta="adaptive", cp=0.99, maxeffort = Inf, options = list()){
  ##processing options
  if (is.null(options$reports)){options$reports <- TRUE}
  else{options$reports <- FALSE}
  if (is.numeric(delta)){options$deltamid <-function(mid){delta}; options$epsilon = delta/100; fixeddelta <- TRUE}
  else if (delta == "adaptive"){
    if (is.null(options$deltamid)){    
      options$deltamid <- mkdeltamid()
      options$epsilon = 0.0001
      fixeddelta <- FALSE
    }
    else {
      if (is.null(options$epsilon)){stop("if using non-default adaptive Delta must specify options$epsilon = lowest point of options$deltamid / 100")}
      fixeddelta <- FALSE
    }
  }
  else {stop("Delta must either be a positive numeric or \"adaptive\"")}

  if (is.null(options$gammapilotprop)){options$gammapilotprop = 0.1}
  if (is.null(options$gammatestprop)){options$gammatestprop = 0.1}
  if (is.null(options$eta)){options$eta = 0.05}
  if (is.null(options$spendgammatest)){
    options$spendgammatest <- function(callno){
      callno/(callno + 10)
    }
  }
  if (is.null(options$maxstepsbase)){options$maxstepsbase = 500} 
  if (is.null(options$maxstepsinc)){options$maxstepsinc = 1.5}
  if (is.null(options$maxbatch)){options$maxbatch = 200000}
  if (is.null(options$batchesbetweentests)){options$batchesbetweentests = 1}
  options$maxeffort <- maxeffort
  options$cp <- cp


  if (is.null(options$pilotn)){options$pilotn <- 1000}
  if (is.null(options$pilotmaxsteps)){options$pilotmaxsteps <- 1000}
  if (options$reports){
    cat("Pilot with ", options$pilotn, " streams stopped after ", options$pilotmaxsteps, " steps... ", sep="")
  }

  pilot <- getpilot(genstream, alpha=alpha, epsilon = options$epsilon, n = options$pilotn, maxsteps = options$pilotmaxsteps)
  if (options$reports){cat("OK\n")}
  res <- pilot$res[, !is.na(pilot$res[2,])]; rescount = c(sum(res[1,] == 0), sum(res[1,]==1))
  pilot$int = getint(rescount=rescount, N=pilot$n, cp=1-(1-cp)*options$gammapilotprop, options$epsilon)
  pilotmid <- (pilot$int[1]+pilot$int[2])/2
  pilotdeltareq <- options$deltamid(pilotmid)
  pilot$effort <- pilot$maxsteps*(sum(is.na(pilot$res[2,]))) + sum(pilot$res[2,], na.rm= TRUE)
  if (diff(pilot$int)<=pilotdeltareq | pilot$effort > maxeffort){
    cat("Est. progress      Conf. Int.\n", sep="")
    cat(100,"%               [", round(pilot$int[1], digits=3),",", round(pilot$int[2], digits=3), "]\n", sep="")

    result <- list()
    result$int <- pilot$int
    result$beta <- sum(pilot$res[1,] == 1)/pilot$n
    result$effort <-pilot$maxsteps*(sum(is.na(pilot$res[2,]))) + sum(pilot$res[2,], na.rm= TRUE)
    result$N <- pilot$n
    result$rescount <- rescount
    result$cp <- options$cp
    result$taccepted <- FALSE
    if (pilot$effort > maxeffort){result$truncated <- TRUE} else {result$truncated <- FALSE}
    return(mkmcpresobj(result))
  }


  if (options$report){cat("Choosing number of streams... ")}
  if (fixeddelta){
    deltapred<-delta
    if (pilot$int[2] < 0.5){    
      minN <- getlowpowerminN(delta = delta, cp=1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), betapilot = pilot$int[2], epsilon=options$epsilon)
      Ntried <- floor(seq(minN, 20*minN, length.out = 500))
      N <- Ntried[which.min(sapply(Ntried, function(n){eeffort(N = n, delta = delta, 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), pilot = pilot)}))]
      if (options$reports){cat(N, "\n")}
    }
    else if (pilot$int[1]>0.5){
      minN <- getlowpowerminN(delta=delta, cp=1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), betapilot = 1-pilot$int[1], options$epsilon)
      Ntried <- floor(seq(minN, 20*minN, length.out = 500))
      N <- Ntried[which.min(sapply(Ntried, function(n){eeffort(N = n, delta = delta, 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), pilot = pilot)}))]
    }
    else {
      minN <- getminN(delta, 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), options$epsilon)
      Ntried <- floor(seq(minN, 20*minN, length.out = 500))
      N <- Ntried[which.min(sapply(Ntried, function(n){eeffort(N = n, delta = delta, cp = 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), pilot = pilot)}))]
    }
  }
  else{#more involved computations needed because almost no properties of Delta known
    getintcall <- function(rescount, N){
      mainint <- getint(rescount = rescount, N = N, cp = 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), options$epsilon)
      int <- c(max(mainint[1], pilot$int[1]), min(mainint[2], pilot$int[2]))
    }
    minN <- getminNfunctional(deltamid=options$deltamid, getintcall= getintcall)
    ##if (options$reports){print(paste("The minimum N: ", minN, sep = ""))}
    Ntried <- floor(seq(minN, 20*minN, length.out = 500))
    beta <- (sum(pilot$res[1,]) + 1)/(pilot$n+2)
    deltapred <- options$deltamid(beta)
    ##if (options$reports){print(paste("We reckon the midpoint will be: ", beta))}
    ##if (options$reports){print(paste("The delta we predict we'll need: ", deltapred, sep = ""))}
    N <- Ntried[which.min(sapply(Ntried, function(n){eeffort(N = n, delta = deltapred, cp = 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), pilot = pilot)}))]
    while(TRUE){
      if (checkN(N = N, deltamid = options$deltamid, getintcall = getintcall)){break}
      N <- max(N+1, ceiling(N*1.01))
    }
    if (options$reports){cat(N, "\n")}
  }

  pilot$esteffort <- pilot$effort + eeffort(N = N, delta = deltapred, cp = 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), pilot = pilot)
  if (options$report){
    progress <- round(pilot$effort/min(pilot$esteffort, options$maxeffort), digits = 2)*100
    digits <- nchar(progress)
    cat("Est. progress      Conf. Int.\n", sep="");
    options$mess <- paste(progress,"%", paste(rep(" ", 18-digits), collapse=""),"[", round(pilot$int[1], digits=3),",", round(pilot$int[2], digits=3), "]", sep="")
    cat(options$mess)
  }

  
  ##Define functions to compute interval:
  getintcallwotest <- function(rescount){
      mainint <- getint(rescount = rescount, N = N, cp = 1-(1-cp)*(1-options$gammapilotprop-options$gammatestprop), options$epsilon)
      int <- c(max(mainint[1], pilot$int[1]), min(mainint[2], pilot$int[2]))
  }  
  callno <- 1
  getintcallwtest <- function(rescount, Ss, nsteps, alg){
    deltahyp <- 1
    deltareq <- 0
    k <- -1
    while (deltahyp > deltareq){#determine what hypotheses to test
      k <- k+1
      if (rescount[1]<=rescount[2]){
        rescounthyp <- rescount+c(floor(k/2), ceiling(k/2))
      }
      else {rescounthyp <- rescount+c(ceiling(k/2), floor(k/2))}
      int <- getintcallwotest(rescounthyp)
      mid <- (int[1] + int[2])/2
      deltareq <- options$deltamid(mid) 
      deltahyp <- diff(int)
    }
    if (!testhyp(rescount=rescount, rescounthyp = rescounthyp, Ss=Ss, eta=options$eta, xi = (1-cp)*options$gammatestprop*(options$spendgammatest(callno+1)-options$spendgammatest(callno)), nsteps = nsteps, alg=alg, reports = options$reports)){int <-getintcallwotest(rescount)
      }
    callno <<- callno + 1
    int
  }

  ##start main loop
  alg <- getalgprecomp(level = alpha, epsilon = options$epsilon)
  result <- getpower(genstream=genstream, alpha=alpha, delta=delta, getintcallwtest=getintcallwtest, getintcallwotest=getintcallwotest, N=N, alg=alg, pilot=pilot, options = options)
}


