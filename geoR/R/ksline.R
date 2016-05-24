"ksline" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            locations, borders = NULL, 
            cov.model = "matern",
            cov.pars = stop("covariance parameters (sigmasq and phi) needed"), 
            kappa = 0.5, nugget = 0, micro.scale = 0,
            lambda = 1, m0 = "ok", nwin = "full", 
            n.samples.backtransform = 500, 
            trend = 1, d = 2, ktedata = NULL, ktelocations = NULL,
            aniso.pars = NULL,  signal = FALSE,  dist.epsilon = 1e-10,
            messages) 
{
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  ##
  ## selecting locations inside the borders 
  ##
  locations <- .check.locations(locations)
  ## Checking for 1D prediction 
  if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
    krige1d <- TRUE
  else krige1d <- FALSE
  if(!is.null(borders)){
    locations <- locations.inside(locations, borders, bound=TRUE)
    if(messages.screen)
      cat("ksline: results will be returned only for prediction locations inside the borders\n")
  }
  call.fc <- match.call()
  cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
  if(abs(lambda-1) > 0.001) {
    if(messages.screen)
      cat("ksline: Box-Cox data transformation performed.\n")
    if(abs(lambda) < 0.001)
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  coords <- as.matrix(coords)
  locations <- as.matrix(locations)
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  if(!is.null(ktedata) & !is.null(ktelocations) & m0 != "kte"){
    cat("ksline: external variable (covariate) provided. Kriging ste to KTE\n")
    m0 <- "kte"
  }
  ##
  ## anisotropy correction
  ##
  if(!is.null(aniso.pars)) {
    if(length(aniso.pars) != 2 | mode(aniso.pars) != "numeric")
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(messages.screen)
      cat("ksline: anisotropy correction performed\n")
    coords.c <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
    locations.c <- coords.aniso(coords = locations, aniso.pars = aniso.pars)
  }
  else {
    coords.c <- coords
    locations.c <- locations
  }
  ## 2. Preparing KTE matrices #####
  ##  
  if(m0 == "kte") {
    ktedata <- as.matrix(ktedata)
    ktelocations <- as.matrix(ktelocations)
    dimnames(ktedata) <- list(NULL, NULL)
    dimnames(ktelocations) <- list(NULL, NULL)
  }
  n <- length(data)
  ni <- length(locations[, 1])
  tausq <- nugget
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  if(nwin == "full") {
    est <- rep(0, ni)
    dif <- rep(0, ni)
    kvar <- rep(0, ni)
    sumw <- rep(0, ni)
    wofmean <- rep(0, ni)
    iv <- varcov.spatial(coords = coords.c, cov.model = cov.model, 
                         kappa = kappa, nugget = nugget, cov.pars = cov.pars, 
                         inv = TRUE, det = FALSE, func.inv = "cholesky")$inverse
    av <- mean(data)
    sd <- sqrt(var(data))
    one <- rep(1, n)
    tone <- t(one)
    toneiv <- crossprod(one, iv)
    den <- solve(toneiv %*% one)
    ml <- den %*% toneiv %*% data
    kmsd <- sqrt(den)
    means <- c(average = av, stdev = sd, kmean = ml, kmsd = kmsd)
    if(m0 != "kt") {
      mktlocations <- "Constant trend"
      beta <- ml
    }
    else {
      mktlocations <- rep(0, ni)
      if(m0 == "kt" & trend == 1) {
        if(d == 1) {
          xmat <- cbind(rep(1, n), coords[, 2])
          xmati <- cbind(rep(1, ni), locations[, 2])
        }
        else {
          xmat <- cbind(rep(1, n), coords[, 1], coords[, 2])
          xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                               2])
        }
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkt <- xmat %*% beta
      }
      if(m0 == "kt" & trend == 2) {
        if(d == 1) {
          xmat <- cbind(rep(1, n), coords[, 2], (coords[, 2])^2)
          xmati <- cbind(rep(1, ni), locations[, 2], (locations[, 
                                                                2])^2)
        }
        else {
          xmat <- cbind(rep(1, n), coords[, 1], coords[, 2], 
                        (coords[, 1])^2, (coords[, 2])^2, coords[, 1] * coords[, 
                                                                               2])
          xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                               2], (locations[, 1])^2, (locations[, 2])^2, locations[, 
                                                                                                                     1] * locations[, 2])
        }
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkt <- xmat %*% beta
      }
    }
    if(m0 != "kte") 
      mktelocations <- "No external trend"
    else {
      if(m0 == "kte") {
        mktelocations <- rep(0, ni)
        xmat <- cbind(rep(1, n), ktedata)
        xmati <- cbind(rep(1, ni), ktelocations)
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkte <- xmat %*% beta
      }
    }
    for (i in 1:ni) {
      if(messages.screen) {
        if(ni < 11) 
          cat(paste("ksline: kriging location: ", i, "out of", 
                    ni, "\n"))
        else {
          if(ni < 101 & (i%%10 == 1)) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if(ni > 100 & i%%100 == 1) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if(i == ni) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
        }
      }
      coords0 <- cbind((coords.c[, 1] - locations.c[i, 1]), (coords.c[, 2] - locations.c[i, 
                                                                                 2]))
      dm0 <- sqrt(coords0[, 1]^2 + coords0[, 2]^2)
      v0 <- cov.spatial(obj = dm0, cov.model = cov.model, 
                        kappa = kappa, cov.pars = cov.pars)
      v0[dm0 < dist.epsilon] <- micro.scale + sigmasq
      tv0 <- t(v0)
      v0iv <- crossprod(v0, iv)
      v0ivv0 <- v0iv %*% v0
      skw <- crossprod(v0,iv)
      wofmean[i] <- 1 - sum(skw)
      ##
      ## 4.2.1 Simple kriging with known mean
      ##
      if(mode(m0) == "numeric") {
        dif[i] <- skw %*% (data - m0)
        est[i] <- m0 + dif[i]
        if(signal == TRUE) 
          kvar[i] <- sigmasq - v0ivv0
        else kvar[i] <- tausq + sigmasq - v0ivv0
        sumw[i] <- sum(skw)
      }
      ##
      ## 4.2.2 Simple kriging with data average mean
      ##
      if(m0 == "av") {
        dif[i] <- skw %*% (data - av)
        est[i] <- av + dif[i]
        if(signal == TRUE) 
          kvar[i] <- sigmasq - v0ivv0
        else kvar[i] <- tausq + sigmasq - v0ivv0
        sumw[i] <- sum(((tone/n) + skw - ((skw %*% one %*% 
                                           tone)/n)))
      }
      ##
      ## 4.2.3 Ordinary kriging (or SK with G.L.S. mean)
      ##
      if(m0 == "ok") {
        dif[i] <- skw %*% (data - ml)
        est[i] <- ml + dif[i]
        redu <- as.vector(1 - toneiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 + (redu %*%
                                         den %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% den %*% redu)
        sumw[i] <- sum((den %*% one + tv0 - v0iv %*% 
                        one %*% den %*% tone) %*% iv)
      }
      ##
      ## 4.2.4 Universal Kriging (or Kriging with trend model) 
      ##
      if(m0 == "kt") {
        dif[i] <- skw %*% (data - mkt)
        est[i] <- xmati[i,  ] %*% beta + dif[i]
        redu <- as.vector(xmati[i,  ]) - as.vector(
                                                   txiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 + (redu %*%
                                         iviv %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% iviv %*% redu)
        sumw[i] <- sum(skw + xmati[i,  ] %*% iviv %*%
                       txiv - skw %*% xmat %*% iviv %*% txiv)
        mktlocations[i] <- xmati[i,  ] %*% beta
      }
      ##
      ## 4.2.5 Kriging with external trend 
      ##
      if(m0 == "kte") {
        dif[i] <- skw %*% (data - mkte)
        est[i] <- xmati[i,  ] %*% beta + dif[i]
        redu <- as.vector(xmati[i,  ]) - as.vector(
                                                   txiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 - (redu %*%
                                         iviv %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% iviv %*% redu)
        sumw[i] <- sum(skw + xmati[i,  ] %*% iviv %*%
                       txiv - skw %*% xmat %*% iviv %*% txiv)
        mktelocations[i] <- xmati[i,  ] %*% beta
      }
      NULL
    }
    message <- "Kriging performed using global neighbourhood"
    if(messages.screen) 
      cat(paste(message,"\n"))
    results <- list(predict = est, krige.var = kvar, dif = dif, summary = means, 
                    ktrend = mktlocations, ktetrend = mktelocations, beta = beta, 
                    wofmean = wofmean)
  }
  else {
    nwin <- min(n, nwin)
    avwin <- rep(0, ni)
    sdwin <- rep(0, ni)
    mlwin <- rep(0, ni)
    kmsdwin <- rep(0, ni)
    estwin <- rep(0, ni)
    difwin <- rep(0, ni)
    kvarwin <- rep(0, ni)
    sumwwin <- rep(0, ni)
    wofmean <- rep(0, ni)
    if(m0 != "kt") 
      mkt <- "Constant position trend"
    else mkt <- rep(0, ni)
    if(m0 != "kte") 
      mkte <- "No external trend"
    else mkte <- rep(0, ni)
    if(m0 != "kt" & m0 != "kte") 
      betawin <- "No polynomial or external trend"
    if(m0 == "kt") {
      if(trend == 1) {
        if(d == 1) 
          xmati <- cbind(rep(1, ni), locations[, 2])
        else xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                                  2])
      }
      if(trend == 2) {
        if(d == 1) 
          xmati <- cbind(rep(1, ni), locations[, 2], locations[, 
                                                               2]^2)
        else xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                                  2], (locations[, 1])^2, (locations[, 2])^2, locations[, 1] * 
                            locations[, 2])
      }
      betawin <- matrix(0, nrow = (ncol(xmati) * ni), ncol = ncol(xmati))
    }
    if(m0 == "kte") {
      xmati <- cbind(rep(1, ni), ktelocations)
      if(is.vector(ktedata) == TRUE) 
        betawin <- matrix(0, nrow = (2 * ni), ncol = 2)
      else betawin <- matrix(0, nrow = ((ncol(ktedata) + 
                                         1) * ni), ncol = (ncol(ktedata) + 1))
    }
    for (i in 1:ni) {
      temp.win <- .ksline.aux.1(coords = coords, coords.c = coords.c,
                               data = data, n = n,
                               locations = locations[i,  ],
                               locations.c = locations.c[i,  ],
                               cov.pars = cov.pars, nugget = nugget,
                               cov.model = cov.model, kappa = kappa, m0 = m0,
                               nwin = nwin, trend = trend, d = d, ktedata = 
                               ktedata, ktelocations = ktelocations,
                               micro.scale = micro.scale, 
                               location.number = i, xmati = xmati[i,  ],
                               mkte = NULL, mkt = NULL, betawin = NULL,
                               signal = signal, dist.epsilon = dist.epsilon)
      avwin[i] <- temp.win$avwin
      sdwin[i] <- temp.win$sdwin
      mlwin[i] <- temp.win$mlwin
      kmsdwin[i] <- temp.win$kmsdwin
      estwin[i] <- temp.win$estwin
      difwin[i] <- temp.win$difwin
      kvarwin[i] <- temp.win$kvarwin
      sumwwin[i] <- temp.win$sumwwin
      wofmean[i] <- temp.win$wofmean
      if(m0 == "kt") 
        mkt[i] <- temp.win$mkt
      if(m0 == "kte") 
        mkte[i] <- temp.win$mkte
      if(m0 == "kt" | m0 == "kte") 
        betawin[i, ] <- temp.win$betawin
      if(messages.screen) {
        if(ni < 11) 
          cat(paste("ksline: kriging location: ", i, "out of", 
                    ni, "\n"))
        else {
          if(ni < 101 & (i%%10 == 1)) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if(ni > 100 & i%%100 == 1) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if(i == ni) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
        }
      }
    }
    message <- "kriging performed in moving neighbourhood"
    if(messages.screen) 
      cat(paste(message,"\n"))
    results <- list(predict = estwin, krige.var = kvarwin, dif = difwin, 
                    avtrend = avwin, sd = sdwin, oktrend = mlwin, oksd = kmsdwin, 
                    ktrend = mkt, ktetrend = mkte, beta = betawin, sumw = sumwwin, 
                    wofmean = wofmean)
  }  
  if(abs(lambda - 1) > 0.001) {
    if(messages.screen)
      cat("Back-transforming the predicted mean and variance.\n")
    if(abs(lambda) < 0.001) {
      predict.transf <- results$predict
      results$predict <- exp(predict.transf) - 0.5 * results$krige.var
      results$krige.var <- (exp(2 * predict.transf - results$krige.var)) * (exp(results$krige.var) - 1)
    }
    if(lambda > 0.001) {
      if(messages.screen)
        cat("Back-transformation by simulating from the normal predictive distribution\n")
      ##ap.warn <- options()$warn
      ##options(warn = -1)
      temp.data <- suppressWarnings(matrix(rnorm(ni * n.samples.backtransform,
                                mean = results$predict,
                                sd = sqrt(results$krige.var)),
                          nrow = ni))
      ##options(warn = ap.warn)
      temp.data[(results$krige.var == 0),  ] <- results$predict[(results$krige.var == 0)]
      temp.data[temp.data < -1/lambda] <- -1/lambda     
      temp.data <- ((temp.data * lambda) + 1)^(1/lambda)
###      temp.data[is.na(temp.data)] <- Inf
      results$predict <- as.vector(rowMeans(temp.data))
      results$krige.var <- as.vector(apply(temp.data, 1, var))
    }
    if(lambda < -0.001) {
      cat("Resulting distribution has no mean for lambda < 0 - back transformation not performed\n"
          )
    }
  }
  results$locations <- locations
  results$message <- message
  results$call <- call.fc
  attr(results, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  attr(results, "prediction.locations") <- call.fc$locations
  if(!is.null(call.fc$borders))
    attr(results, "borders") <- call.fc$borders
  oldClass(results) <- c("kriging")
  return(invisible(results))
}

".ksline.aux.1" <-
  function (coords, coords.c, data, n, locations, locations.c, cov.pars,
            nugget, cov.model, kappa, 
            m0, nwin, trend, d, ktedata, ktelocations, mbased,
            micro.scale, location.number, 
            xmati, mkte, mkt, betawin, signal, dist.epsilon) 
{
  i <- location.number
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  tausq <- nugget
  coords0 <- cbind((coords.c[, 1] - locations.c[1]), (coords.c[, 2] -
                                                      locations.c[2]))
  dm0 <- sqrt(coords0[, 1]^2 + coords0[, 2]^2)
  coordswin <- coords[order(dm0)[1:nwin],  ]
  coordswin.c <- coords.c[order(dm0)[1:nwin],  ]
  datawin <- data[order(dm0)[1:nwin]]
  ivwin <- varcov.spatial(coords = coordswin.c, cov.model = cov.model,
                          kappa = kappa, nugget = nugget, cov.pars = cov.pars, inv = TRUE,
                          det = FALSE, func.inv = "cholesky", only.decomposition = FALSE)$inverse
  avwin <- mean(datawin)
  sdwin <- sqrt(var(datawin))
  onewin <- rep(1, nwin)
  toneivwin <- crossprod(onewin, ivwin)
  denwin <- solve(toneivwin %*% onewin)
  mlwin <- denwin %*% toneivwin %*% datawin
  kmsdwin <- sqrt(denwin)
  coords0win <- cbind((coordswin[, 1] - locations[1]), (coordswin[, 2] -
                                                        locations[2]))
  coords0win.c <- cbind((coordswin.c[, 1] - locations.c[1]), (coordswin.c[
                                                                          , 2] - locations.c[2]))
  dm0win <- sqrt(coords0win.c[, 1]^2 + coords0win.c[, 2]^2)
  v0win <- cov.spatial(obj = dm0win, cov.model = cov.model, kappa = kappa,
                       cov.pars = cov.pars)
  v0win[dm0win < dist.epsilon] <- micro.scale + sigmasq
  skwwin <- crossprod(v0win, ivwin)
  wofmean <- 1 - sum(skwwin)
  if(m0 == "kt" & trend == 1) {
    if(d == 1)
      xmatwin <- cbind(rep(1, nwin), coordswin[, 2])
    else xmatwin <- cbind(rep(1, nwin), coordswin[, 1], coordswin[
                                                                  , 2])
    txivwin <- crossprod(xmatwin, ivwin)
    ivivwin <- solve(txivwin %*% xmatwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    mktwin <- xmatwin %*% betawin
    mkt <- xmati %*% betawin
  }
  if(m0 == "kt" & trend == 2) {
    if(d == 1)
      xmatwin <- cbind(rep(1, nwin), coordswin[, 2], (
                                                      coordswin[, 2])^2)
    else xmatwin <- cbind(rep(1, nwin), coordswin[, 1], (coordswin[
                                                                   , 1])^2, coordswin[, 2], (coordswin[, 2])^
                          2, coordswin[, 1] * coordswin[, 2])
    xmatwin.cent <- xmatwin
    xmatwin.cent[, 2] <- xmatwin.cent[, 2] - mean(xmatwin[, 2])
    xmatwin.cent[, 3] <- xmatwin.cent[, 3] - mean(xmatwin[, 3])
    ivivwin <- solve(crossprod(xmatwin.cent, ivwin) %*% 
                     xmatwin.cent)
    txivwin <- crossprod(xmatwin.cent, ivwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    betawin <- mean(datawin) - crossprod(betawin, c(0, mean(xmatwin[
                                                                    , 2]), mean(xmatwin[, 3])))
    mktwin <- xmatwin %*% betawin
    mkt <- xmati %*% betawin
  }
  if(m0 == "kte") {
    if(is.vector(ktedata))
      ktedatawin <- ktedata[order(dm0)[1:nwin]]
    else ktedatawin <- ktedata[order(dm0)[1:nwin],  ]
    xmatwin <- cbind(rep(1, nwin), ktedatawin)
    ivivwin <- solve(crossprod(xmatwin, ivwin) %*% xmatwin)
    txivwin <- crossprod(xmatwin, ivwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    mktewin <- xmatwin %*% betawin
    mkte <- xmati %*% betawin
  }
  ##
  ##  Simple kriging with know mean
  ##
  if(mode(m0) == "numeric") {
    difwin <- skwwin %*% (data - m0)
#    estwin <- m0win + difwin
    estwin <- m0 + difwin
    if(signal)
      kvarwin <- sigmasq - crossprod(v0win, ivwin) %*% v0win
    else kvarwin <- tausq + sigmasq - crossprod(v0win, ivwin) %*%
      v0win
    sumwwin <- sum(skwwin)
  }
  ##
  ## 4.2.2 Simple kriging with data average mean
  ##
  if(m0 == "av") {
    difwin <- skwwin %*% (datawin - avwin)
    estwin <- avwin + difwin
    if(signal)
      kvarwin <- sigmasq - crossprod(v0win, ivwin) %*% v0win
    else kvarwin <- tausq + sigmasq - crossprod(v0win, ivwin) %*%
      v0win
    sumwwin <- sum(((t(onewin)/nwin) + skwwin - ((skwwin %*% onewin %*%
                                                  t(onewin))/n)))
  }
  ##
  ## Ordinary kriging (or SK with G.L.S. mean)
  ##
  if(m0 == "ok") {
    difwin <- skwwin %*% (datawin - mlwin)
    estwin <- mlwin + difwin
    redu <- as.vector(1 - toneivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - v0win %*% ivwin %*% v0win + (
                                                        redu %*% denwin %*% redu)
    else kvarwin <- tausq + sigmasq - v0win %*% ivwin %*% v0win +
      (redu %*% denwin %*% redu)
    sumwwin <- sum((denwin %*% onewin + t(v0win) - crossprod(v0win,
                                                             ivwin) %*% onewin %*% denwin %*% t(onewin)) %*% ivwin)
  }
  ##
  ## Universal Kriging (or Kriging with trend model) 
  ##
  if(m0 == "kt") {
    difwin <- skwwin %*% (datawin - mktwin)
    estwin <- mkt + difwin
    xmati <- as.vector(xmati)
    redu <- as.vector(xmati) - as.vector(txivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    else kvarwin <- tausq + sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    sumwwin <- sum(skwwin + xmati %*% ivivwin %*% txivwin - skwwin %*%
                   xmatwin %*% ivivwin %*% txivwin)
  }
  ##
  ## Kriging with external trend 
  ##
  if(m0 == "kte") {
    difwin <- skwwin %*% (datawin - mktewin)
    estwin <- mkte + difwin
    xmati <- as.vector(xmati)
    redu <- as.vector(xmati) - as.vector(txivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    else kvarwin <- tausq + sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    sumwwin <- sum(skwwin + xmati %*% ivivwin %*% txivwin - skwwin %*%
                   xmatwin %*% ivivwin %*% txivwin)
  }
  ##
  ##
  ##  
  results <- list(avwin = avwin, sdwin = sdwin, mlwin = mlwin, kmsdwin = 
                  kmsdwin, wofmean = wofmean, betawin = betawin, mkt = mkt, mkte
                  = mkte, difwin = difwin, estwin = estwin, kvarwin = kvarwin,
                  sumwwin = sumwwin)
  return(results)
}












