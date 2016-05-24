ftide <- function(x, dto, hcn = hc60, astlon = c("task","cartwright"),
      nodal = TRUE, smsl = FALSE, span = 0.75, degree = 1, ...)
{
  astlon <- match.arg(astlon)
  if(missing(x) || !is.numeric(x)) 
    stop("'x' must be a numeric vector or ts object")
  if(length(x) != length(dto))
    stop("'x' and 'dto' must have the same length")
  if(any(is.na(dto))) {
    stop("missing values in time/date object")
  }
  
  # match either naming scheme
  namei <- match(hcn, harmonics$name)
  namei2 <- match(hcn, harmonics$sname)
  namei[is.na(namei)] <- namei2[is.na(namei)]
  
  # try all upper case
  if(any(is.na(namei))) {
    nameui <- match(hcn, toupper(harmonics$name))
    nameui2 <- match(hcn, toupper(harmonics$sname))
    nameui[is.na(nameui)] <- nameui2[is.na(nameui)]
    namei[is.na(namei)] <- nameui[is.na(namei)]
  }
  
  # try all lower case
  if(any(is.na(namei))) {
    nameui <- match(hcn, tolower(harmonics$name))
    nameui2 <- match(hcn, tolower(harmonics$sname))
    nameui[is.na(nameui)] <- nameui2[is.na(nameui)]
    namei[is.na(namei)] <- nameui[is.na(namei)]
  }
  
  # try full greek letters
  # also try V for nu swap
  if(any(is.na(namei))) {
    harname <- harmonics$name
    harsname <- harmonics$sname
    harname <- sub("sig","sigma",harname)
    harname <- sub("lam","lambda",harname)
    harname <- sub("the","theta",harname)
    harname <- sub("gam","gamma",harname)
    harname <- sub("alp","alpha",harname)
    harname <- sub("del","delta",harname)
    harname <- sub("nu","XXX",harname)
    harname <- sub("V","nu",harname)
    harname <- sub("XXX","V",harname)
    harname <- sub("D","lambda",harname) 
    harsname <- sub("sig","sigma",harsname)
    harsname <- sub("lam","lambda",harsname)
    harsname <- sub("the","theta",harsname)
    harsname <- sub("gam","gamma",harsname)
    harsname <- sub("alp","alpha",harsname)
    harsname <- sub("del","delta",harsname)
    harsname <- sub("nu","XXX",harsname)
    harsname <- sub("V","nu",harsname)
    harsname <- sub("XXX","V",harsname)
    harsname <- sub("D","lambda",harsname)
  }
  
  # try all upper case with full greek
  if(any(is.na(namei))) {
    nameui <- match(hcn, toupper(harname))
    nameui2 <- match(hcn, toupper(harsname))
    nameui[is.na(nameui)] <- nameui2[is.na(nameui)]
    namei[is.na(namei)] <- nameui[is.na(namei)]
  }
  
  # try all lower case with full greek
  if(any(is.na(namei))) {
    nameui <- match(hcn, tolower(harname))
    nameui2 <- match(hcn, tolower(harsname))
    nameui[is.na(nameui)] <- nameui2[is.na(nameui)]
    namei[is.na(namei)] <- nameui[is.na(namei)]
  }
  
  # try special cases
  if(any(is.na(namei))) {
    hcn[hcn == "RHO" | hcn == "rho"] <- "rho1" 
    hcn[hcn == "lamO1" | hcn == "LAMO1"] <- "1D.1o1"
    hcn[hcn == "lamdaO1" | hcn == "LAMDAO1"] <- "1D.1o1"
    hcn[hcn == "lamda2" | hcn == "LAMDA2"] <- "lam2"
    hcn[hcn == "Msf"] <- "MSf"
    hcn[hcn == "2N4"] <- "N4"
    namei2 <- match(hcn, harmonics$name)
    namei[is.na(namei)] <- namei2[is.na(namei)]
  }
  
  # not recognized
  if(any(is.na(namei))) 
    stop(paste(c("names not recognized",hcn[is.na(namei)]),collapse=" "))
  hcn <- harmonics$name[namei]
  
  nc <- length(hcn)
  nt <- length(x)
  speed <- harmonics$speed[namei]
  pcon <- harmonics$phi[namei]
  dood <- cbind(harmonics$i1[namei],
                harmonics$i2[namei]-harmonics$i1[namei],
                harmonics$i3[namei]+harmonics$i1[namei],
                harmonics$i4[namei],harmonics$i5[namei],harmonics$i6[namei])
  
  # similar speeds check
  spddiff <- diff(sort(speed))
  if(any(spddiff < 1e-05))
    stop("two components at similar speeds detected")
  
  # mean sea level
  z0 <- mean(x, na.rm = TRUE)
  x <- x - z0
  
  # time zone issues
  times <- as.POSIXct(dto, tz = "UTC")
  attr(times, "tzone") <- "UTC"
  
  # origin and hours since origin
  orig <- times[length(times) %/% 2]
  hours <- as.numeric(difftime(times, orig, units = "hours"))
  
  # smooth mean sea level
  if(smsl) {
    z1 <- loess(x ~ hours, na.action = na.exclude, 
      trace.hat = "approximate", span = span, degree = degree)
    x <- x - fitted(z1)
    z0 <- z0 + fitted(z1)
  }
  else z1 <- NULL
  
  # days of period, lambdas and day indicator
  days <- seq(as.Date(times[1]), as.Date(times[nt]), 1)
  nd <- length(days)
  lamb <- lambdas(days, astlon = astlon)
  dind <- as.numeric(floor(julian(times, origin = as.Date(times[1])) + 1))
  
  # reference phase
  t0 <- which(days == as.Date(orig))
  vn0 <- (as.numeric(dood[,-1] %*% lamb[,t0,drop=FALSE]) + pcon) %% 360
  
  # nodal adjustments
  if(nodal) {
    ndl <- nodal_adj(lamb[3,], lamb[4,], lamb[5,])
    fn <- ndl$fn[hcn,,drop=FALSE]
    un <- ndl$un[hcn,,drop=FALSE]
  } else {
    fn <- matrix(1, nrow = nc, ncol = nd)
    un <- matrix(0, nrow = nc, ncol = nd)
  }
  
  # design matrix for lm fit
  xmat <- matrix(NA, nrow = nt, ncol = 2*nc)
  colnames(xmat) <- c(paste0(hcn, "_S"), paste0(hcn, "_C"))
  for(i in 1:nc) {
    ff <- fn[i,dind]
    targ <- pi/180 * (speed[i]*hours + un[i,dind] + vn0[i])
    xmat[,i] <- ff*sin(targ)
    xmat[,i+nc] <- ff*cos(targ)
  }

  xmat <- as.data.frame(xmat)
  colnames(xmat) <- c(paste0(hcn, "_S"), paste0(hcn, "_C"))
  lmfit <- lm(x ~ -1 + ., data = xmat, na.action = na.exclude, ...)
  cfmat <- matrix(coef(lmfit), nrow = nc)
  dimnames(cfmat) <- list(hcn, c("sine", "cosine"))
  
  apmat <- cbind(sqrt(rowSums(cfmat^2)), (180/pi * atan2(cfmat[,1],cfmat[,2])) %% 360)
  apmat <- apmat[sort.list(apmat[,1], decreasing=TRUE),]
  colnames(apmat) <- c("amplitude", "phase")
  
  fval <- NA
  if(all(c("M2","S2","K1","O1") %in% hcn)) {
    amp <- apmat[,1]
    MSL <- mean(z0, na.rm = TRUE)
    # form factor #0.25-3 mixed #3+ diurnal #0.5
    fval <- as.numeric((amp["K1"] + amp["O1"])/(amp["M2"] + amp["S2"]))
    
    # diurnal and mixed
    MHHW <- MSL + (amp["M2"] + amp["K1"] + amp["O1"])
    MLHW <- MSL + abs(amp["M2"] - amp["K1"] - amp["O1"])
    MHLW <- MSL - abs(amp["M2"] - amp["K1"] - amp["O1"])
    MLLW <- MSL - (amp["M2"] + amp["K1"] + amp["O1"])
    
    # semi-diurnal
    MHWS <- MSL + (amp["M2"] + amp["S2"])
    MHWN <- MSL + abs(amp["M2"] - amp["S2"])
    MLWN <- MSL - abs(amp["M2"] - amp["S2"])
    MLWS <- MSL - (amp["M2"] + amp["S2"])
    
    lmfit$features1 <- c(MLWS, MLWN, MSL, MHWN, MHWS)
    names(lmfit$features1) <- c("MLWS", "MLWN", "MSL", "MHWN", "MHWS")
    lmfit$features2 <- c(MLLW, MHLW, MSL, MLHW, MHHW)
    names(lmfit$features2) <- c("MLLW", "MHLW", "MSL", "MLHW", "MHHW")
  }
  
  lmfit$msl <- z0
  lmfit$cfmat <- cfmat
  lmfit$apmat <- cbind(apmat, cfmat[match(rownames(apmat), rownames(cfmat)),])
  lmfit$fval <- fval
  lmfit$nodal <- nodal
  lmfit$astlon <- astlon
  lmfit$orig <- orig
  lmfit$vn0 <- vn0
  names(lmfit$vn0) <- hcn
  lmfit$lobj <- z1
  class(lmfit) <- c("tide", "lm")
  lmfit
}

plagtz <- function(plag, tzd, indegree = TRUE, outdegree = TRUE) 
{
  hcn <- names(plag)
  namei <- match(hcn, harmonics$name)
  if(any(is.na(namei))) 
    stop("names in 'plag' are not recognized")
  speed <- harmonics$speed[namei]
  if(!indegree) plag <- plag * 180/pi
  plag <- (plag + speed*tzd) %% 360
  if(!outdegree) plag <- plag * pi/180
  plag
}

coef.tide <- function(object, hc = FALSE, mat = hc, utc = 0, ...) 
{
  if(!hc && !mat) {
    out <- object$coefficients
    nms <- rownames(object$cfmat)
    nms <- paste0(rep(nms,2), rep(c("_S","_C"), each = length(nms)))
    names(out) <- nms
  }
  if(!hc && mat) out <- object$cfmat
  if(hc && mat) {
    out <- object$apmat
    if(utc != 0) out[,"phase"] <- plagtz(out[,"phase"], utc)
  }
  if(hc && !mat) {
    nms <- rownames(object$apmat)
    nms <- paste0(rep(nms,2), rep(c("_A","_P"), each = length(nms)))
    out <- object$apmat[,1:2]
    if(utc != 0) out[,"phase"] <- plagtz(out[,"phase"], utc)
    out <- as.numeric(out)
    names(out) <- nms
  }
  out
}

print.tide <- function(x, digits = max(3L, getOption("digits") - 3L), ...) 
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if(!is.na(x$fval)) {
    cat("Form Factor:", x$fval, "\n\n")
    if(x$fval < 0.5) {
      cat("Features (Semi-Diurnal):\n")
      print.default(x$features1)
    } else {
      cat("Features (Diurnal or Mixed):\n")
      print.default(x$features2)
    }
  }
  cat("\n")
  if (length(coef(x))) {
    cat("Harmonics:\n")
    print.default(format(coef(x, hc = TRUE), digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No harmonics\n")
  cat("\n")
  invisible(x)
}
  
predict.tide <- function(object, from, to, by = NULL, split = FALSE, which = NULL, msl = !split, ...)
{
  hc <- coef(object, hc = TRUE)[,1:2]
  astlon <- object$astlon
  nodal <- object$nodal
  
  if(!is.null(which)) {
    if(!is.character(which))
      stop("'which' should be a character vector")
    if(!all(which %in% rownames(hc)))
      stop("components not included in harmonics object")
    hc <- hc[which,,drop=FALSE]
  }
  ampl <- hc[,1]
  plag <- hc[,2]
  hcn <- rownames(hc)
  nc <- length(hcn)
  namei <- match(hcn, harmonics$name)
  
  # times vector
  t1 <- as.POSIXct(from, tz = "UTC")
  t2 <- as.POSIXct(to, tz = "UTC")
  if(t1 >= t2) stop("'to' must be later than 'from'")
  if(is.null(by)) by <- as.numeric(difftime(t2,t1,units="hours")/999)
  times <- seq(t1, t2, by = by * 3600)
  attr(times, "tzone") <- "UTC"
  nt <- length(times)
  
  # mean sea level
  if(!is.null(object$lobj)) {
    nhours <- as.numeric(difftime(times, object$orig, units = "hours"))
    z0 <- predict(object$lobj, newdata = nhours) + mean(object$msl, na.rm=TRUE)
  } else z0 <- object$msl
  
  # speed, pcon, dood, origin and hours since origin
  speed <- harmonics$speed[namei]
  pcon <- harmonics$phi[namei]
  dood <- cbind(harmonics$i1[namei],
                harmonics$i2[namei]-harmonics$i1[namei],
                harmonics$i3[namei]+harmonics$i1[namei],
                harmonics$i4[namei],harmonics$i5[namei],harmonics$i6[namei])
  hours <- as.numeric(difftime(times, object$orig, units = "hours"))
  
  # days of period, lambdas and day indicator
  days <- seq(as.Date(times[1]), as.Date(times[nt]), 1)
  nd <- length(days)
  lamb <- lambdas(days, astlon = astlon)
  dind <- as.numeric(floor(julian(times, origin = as.Date(times[1])) + 1))
  
  # reference phase
  vn0 <- object$vn0[hcn]
  
  # nodal adjustments
  if(nodal) {
    ndl <- nodal_adj(lamb[3,], lamb[4,], lamb[5,])
    fn <- ndl$fn[hcn,,drop=FALSE]
    un <- ndl$un[hcn,,drop=FALSE]
  } else {
    fn <- matrix(1, nrow = nc, ncol = nd)
    un <- matrix(0, nrow = nc, ncol = nd)
  }
  
  predmat <- matrix(NA, nc, nt)
  rownames(predmat) <- hcn
  for (j in 1:nt) {
    for (i in 1:nc) {
      predmat[i,j] <- fn[i,dind[j]] * ampl[i] * 
        cos(pi/180 * (speed[i] * hours[j] + un[i,dind[j]] + vn0[i] - plag[i]))
    }
  }
  
  preds <- colSums(predmat)
  if(msl) {
    preds <- preds + z0
    predmat <- predmat + matrix(z0, nrow=nc, ncol=nt, byrow=TRUE)
  }
  
  if(split) preds <- predmat 
  preds
}
  
plot.tide <- function(x, from, to, by = NULL, split = FALSE, which = NULL, msl = !split, 
  ask = split && dev.interactive(), main = NULL, xlab = "Times", ylab = "Level",...)
{
  # times and preds
  t1 <- as.POSIXct(from, tz = "UTC")
  t2 <- as.POSIXct(to, tz = "UTC")
  if(t1 >= t2) stop("'to' must be later than 'from'")
  if(is.null(by)) by <- as.numeric(difftime(t2,t1,units="hours")/999)
  times <- seq(t1, t2, by = by * 3600)
  preds <- predict(x, from, to, by = by, split = split, which = which, msl = msl) 
  
  opar <- par(ask = ask)
  if(split) {
    hcn <- rownames(preds)
    for(i in 1:length(hcn)) {
      if(is.null(main)) maini <- hcn[i]
      plot(times, preds[i,], main = maini, xlab = xlab, ylab = ylab, type = "l", ...)
    }
  }
  else {
    plot(times, preds, main = main, xlab = xlab, ylab = ylab, type = "l", ...)
  }
  par(opar)
  invisible(list(times = times, preds = preds))
}

