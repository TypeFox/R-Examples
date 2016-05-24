
gdrate <- function(input, pval, plots) {

  # Function for given model and dataset
  gdX <- function(input1, i) {
    tit <- paste("ID=", unique(input1$name), sep = "")
    dset <- input1[order(input1$date), ]
    f <- as.matrix(dset$size)
    f = f/f[1]
    time <- as.matrix(dset$date)
    time = (time - time[1]) + 1
    jdta <- data.frame(cbind(time, f))
    colnames(jdta) <- c("time", "f")
    v <- subset(foo, foo$IDmodel == i)

    try({
      outgd <- nlsLM(eval(parse(text = paste(v$model))), data = jdta, start = eval(parse(text = paste(v$start))),
                     control = nls.lm.control(maxiter = 1000, maxfev = 1000, factor = 0.01,
                                              ftol = sqrt(.Machine$double.eps),
                                              ptol = sqrt(.Machine$double.eps)),
                     lower = eval(parse(text = paste(v$lb))),
                     upper = eval(parse(text = paste(v$ub))))
    }, silent = TRUE)

    if(!exists("outgd") && i == 4 ){
      try({
        outgd2 <- stats::nls(eval(parse(text = paste(v$model))), data = jdta,
                      start = eval(parse(text = paste(v$start))),
                      algorithm = 'port',
                      control = stats::nls.control(maxiter = 1000, warnOnly = FALSE, minFactor = .000001),
                      lower = eval(parse(text = paste(v$lb))),
                      upper = eval(parse(text = paste(v$ub))))
      }, silent = TRUE)
      return(outgd2)
    } else {
      return(outgd)
    }
  }

  # Function to prepare user input data for modeling
  inputprep <- function(input1) {

    if (is.null(input1)) {
      stop("input argument missing")
    } else {
      try({
        input <- input1[stats::complete.cases(input1), ]
      }, silent = TRUE)
      if (dim(input)[1] < 1) {
        stop("input contains no non-missing data")
      } else {
        if (!c('name') %in% colnames(input) | !c('size') %in% colnames(input) | !c('date') %in% colnames(input)) {
          stop("please rename columns as described in help page")
        }  else {
          if (!is.numeric(input[,c(1)]) | !is.numeric(input[,c(2)]) | !is.numeric(input[,c(3)])) {
            stop("all input data must be numeric")
          } else {
            input2 <- input
          }
        }
      }
    }

    name <- unique(input2[, "name"])
    lasti <- as.numeric(length(name))

    # return number of row entries, evaluation data points by patient
    la <- length(unique(input1$name))
    idall <- rep(1:la, 1)
    id <- cbind(idall, unique(input1$name))
    colnames(id) <- c("idall", "name")
    input0 <- merge(id, input1, by = "name", all = TRUE)
    input00 <- input0[!is.na(input0$name) & !is.na(input0$date), ]

    countent <- function(i) {
      # number of entries by patient
      r <- subset(input0, input0$idall == i)
      nrow <- dim(r)[1]
      name <- as.numeric(paste(unique(r$name)))

      # number of datapoints (ie size and date) by patient
      r2 <- subset(input00, input00$idall == i)
      nc <- dim(r2)[1]
      if (nc > 0) {
        numcyc <- nc
        # number of unique f values for a patient to identify error data
        i3 <- as.data.frame(unique(r2[, c("name", "size")]))
        nunique <- dim(i3)[1]
        nout <- cbind(numcyc, nunique)
      } else {
        nout <- cbind(numcyc = NaN, nunique = NaN)
      }
      cout <- cbind(name, nrow, nout)
    }
    ninfo <- data.frame(do.call("rbind", sapply(1:la, countent, simplify = FALSE)))
    IDm <- rep(1:lasti, 1)
    nameid <- data.frame(cbind(name, IDm))
    x3 <- merge(nameid, input2, by = "name")

    countD <- function(i) {
      r <- subset(x3, x3$IDm == i)

      # time and size initial and last
      r1 <- r[order(r$date), ]
      dr1 <- dim(r1)[1]
      r$time0 <- r1[1, 3]
      r$size0 <- r1[1, 4]
      r$timelast <- r1[dr1, 3]
      r$sizelast <- r1[dr1, 4]

      # time and size at nadir
      r2 <- r[order(r$size), ]
      r$timemin <- r2[1, 3]
      r$sizemin <- r2[1, 4]

      # t and f calc for model input
      r$t <- r$date - r$time0
      # if size initial =0
      if (r1[1, 4] == 0) {
        r$f <- (r$size + 0.5)/(r$size0 + 0.5)
      } else {
        r$f <- r$size/r$size0
      }

      r5 <- subset(r, r$f != "NaN" & r$f != "Inf")
      # total time evaluated
      r5$tseq <- r5$timelast - r5$time0
      r5
    }
    x5 <- data.frame(do.call("rbind", sapply(1:lasti, countD, simplify = FALSE)))

    # merge ninfo with input details
    ninfo2 <- merge(x5, ninfo, by = "name", all = TRUE)
    ninfo2$twock <- ifelse((ninfo2$size0 == 0), (ninfo2$sizelast + 0.5)/(ninfo2$size0 +
                                                                           0.5), ninfo2$sizelast/ninfo2$size0)

    # assign type- either analyze or reason for exclusion
    ninfo2$calcfinal <- ifelse((ninfo2$nunique == "NA" | is.na(ninfo2$nunique) |
                                  ninfo2$size0 == 0), "No measurement data", ifelse((ninfo2$numcyc == 1),
                                                                                    "only 1 eval", ifelse((ninfo2$nunique == 1 & ninfo2$numcyc > 2), "error data",
                                                                                                          ifelse((ninfo2$size0 == 0 & ninfo2$sizelast == 0), "error data",
                                                                                                                 ifelse((ninfo2$numcyc == 2 & ninfo2$twock > 0.8 & ninfo2$twock <
                                                                                                                           1.2), "2 evals not 20% diff", ifelse((ninfo2$numcyc == 2 & ninfo2$twock <=
                                                                                                                                                                   0.8 | ninfo2$twock >= 1.2), "analyze", "analyze"))))))


    ex <- subset(ninfo2, ninfo2$calcfinal != "analyze")
    if (dim(ex)[1] > 0) {
      zexc <- unique(ex[, c("IDm", "name", "calcfinal", "numcyc")])
    } else {
      zexc <- data.frame(NaN, NaN, NaN, NaN)
      colnames(zexc) <- c("IDm", "name", "calcfinal", "numcyc")
    }

    ninfo3 <- ninfo2[, c(1:16, 18)]
    z8a <- subset(ninfo3, ninfo3$calcfinal == "analyze")
    z8 <- z8a[, c(2, 1, 3:13, 15)]
    name <- unique(z8$name)
    nl <- length(name)

    # return result if data has only excluded cases, else input for analysis
    if (nl < 1) {
      resultip <- list(excluded = zexc, inputdata = zexc, cont = 0)
    } else {
      IDr <- rep(1:nl, 1)
      z8id <- cbind(name, IDr)
      z9a <- merge(z8, z8id, by = "name", all = TRUE)
      z9 <- z9a[, c(2, 1, 3:ncol(z9a))]
      resultip <- list(excluded = zexc, inputdata = z9, cont = 1)
    }
    return(resultip)
  }

  # Function to plot observed and predicted values for given patient and model
  plotgdX <- function(input1, i) {
    # data
    tit <- paste("ID=", unique(input1$name), sep = "")
    dset <- input1[order(input1$date), ]
    tseq <- as.numeric(paste(unique(dset$tseq)))
    time <- dset$t
    f <- dset$f
    jdta <- data.frame(cbind(time, f))

    # model info
    v <- subset(foo, foo$IDmodel == i)
    ft <- paste(v$fit)
    cc <- paste(v$cc)

    # model given input and i
    outgd <- gdX(input1, i)
    newx <- seq(1, tseq, by = 1)
    dnew <- data.frame(time = newx)
    #prd <- stats::predict(outgd, newdata = dnew)
    #length(seq(1, tseq, by = 1))
    prd <- stats::predict(outgd, newdata = data.frame(time = newx))

    # merge pred with input for calc rmse
    yhat <- cbind(newx, prd)
    colnames(yhat)[1] <- "time"
    h <- merge(jdta, yhat, by = "time")
    h$remove <- ifelse((h$time == 0 & h$prd != 1), 1, 0)
    h2 <- subset(h, h$remove == 0)
    h2$res <- h2$f - h2$prd
    h2$resSq <- h2$res * h2$res
    rmse <- sqrt(mean(h2$resSq))

    # plot
    #par(mar = c(6.5, 4.5, 1, 1.5))
    graphics::par(mar = c(6.5, 4.5, 1, 1.5))
    graphics::plot(f ~ time, data = jdta, frame = FALSE, col = "red", cex = 1.3, cex.axis = 1.4,
                   cex.lab = 1.6, pch = 19, xlab = "Days", ylab = "Tumor Q/Q0", main = tit)
    #lines(newx, prd, col = cc, lty = 1, lwd = 3)
    graphics::lines(newx, prd, col = cc, lty = 1, lwd = 3)

    lp <- ifelse((ft == "dx"), "topright", "topleft")
    #legend(lp, ft, col = cc, bty = "n", lty = c(1), lwd = 3, cex = 1.2)
    graphics::legend(lp, ft, col = cc, bty = "n", lty = c(1), lwd = 3, cex = 1.2)


    # observed values
    graphics::points(f ~ time, data = jdta, pch = 21, col = c("black"), bg = "red", lwd = 1.2,
                     cex = 1.5)
    # return(rmse)
  }

  # Function to compare models and return selected fit with estimates or not fit
  checkconv2 <- function() {

    # empty mod row if lm LT np or model not fit
    xfit <- function(fit, name00, iMod, lm) {
      isconv <- "NA"
      stopMessage <- "NA"
      stopcode <- "NA"
      zout <- cbind(fit, iMod, name00, stopcode, stopMessage, isconv)
      laout <- cbind(name00, sigp = NaN, np = NaN, LL = NaN, AIC = NaN, AICc = NaN,
                     lm)
      zaout <- merge(zout, laout, by = "name00")
      zaout
    }

    # empty coef tab if lm LT np or model not fit
    xcof <- function(fit, name00) {
      parameter <- "NA"
      modelnames <- fit
      name00 <- name00
      q0 <- data.frame(cbind(Estimate = NaN, Std..Error = NaN, t.value = NaN,
                             Pr...t.. = NaN))
      q1 <- cbind(name00, parameter, modelnames)
      q <- cbind(q0, q1)
      q
    }

    # by patient
    fid4 <- function(k) {
      input1a <- subset(c, c$ID4 == k)

      # by model given patient k
      fmod <- function(i) {
        name00 <- as.numeric(paste(unique(input1a$name)))
        fit <- paste(foo[foo$IDmodel == i, c(2)])
        iMod <- as.numeric(paste(foo[foo$IDmodel == i, c(1)]))

        # number of parameters in model
        np <- as.numeric(paste(foo[foo$IDmodel == i, c(9)]))

        # number of measurement values
        lm <- length(input1a$size)

        # cond exe if n measurement values LE n parameters
        if (lm <= np) {
          name00 <- as.numeric(paste(unique(input1a$name)))
          zaout0 <- xfit(fit = fit, name00 = name00, iMod = iMod, lm = lm)
          zcof <- xcof(fit = fit, name00 = name00)
          zaout <- merge(zaout0, zcof, by = "name00", all = TRUE)
          zaout
        } else {
          try({
            outgd <- gdX(input1a, i)
          }, silent = TRUE)

          if (exists("outgd")) {
            z <- outgd$convInfo
            stopcode <- z$stopCode
            stopMessage <- z$stopMessage
            isconv <- z$isConv
            zout <- cbind(fit, iMod, name00, stopcode, stopMessage, isconv)

            if (isconv == "TRUE") {
              LL <- stats::logLik(outgd)
              AIC <- as.numeric(paste(-2 * LL + 2 * np))
              AICc <- as.numeric(paste(AIC + 2 * np * (np + 1)/(lm - np -
                                                                  1)))
              q <- data.frame(summary(outgd)$coefficients)
              q$name00 <- name00
              q$parameter <- row.names(q)
              q$modelnames <- fit
              q$sig <- ifelse((q$Pr...t.. < pval), 1, 0)
              sigp <- sum(q$sig)
              laout <- data.frame(cbind(name00, sigp, np, LL, AIC, AICc,
                                        lm))
              zaout0 <- merge(zout, laout, by = "name00")
              zcof <- data.frame(q[, c(1:7)])
              zaout <- merge(zaout0, zcof, by = "name00", all = TRUE)
              zaout
            } else {
              zout <- cbind(fit, iMod, name00, stopcode, stopMessage, isconv)
              laout <- data.frame(cbind(name00, sigp = NaN, np, LL = NaN,
                                        AIC = NaN, AICc = NaN, lm))
              zaout0 <- merge(zout, laout, by = "name00")
              zcof <- xcof(fit = fit, name00 = name00)
              zaout <- merge(zaout0, zcof, by = "name00", all = TRUE)
              zaout
            }
          } else {
            zaout0 <- xfit(fit = fit, name00 = name00, iMod = iMod, lm = lm)
            zcof <- xcof(fit = fit, name00 = name00)
            zaout <- merge(zaout0, zcof, by = "name00", all = TRUE)
            zaout
          }  #end if exists outgd
          zaout
        }  #end if n params Lt n measurements and size0 NE 0
        return(zaout)
      }

      fmd <- do.call("rbind", sapply(1:nm, fmod, simplify = FALSE))

      # selected or not fit
      sigm <- unique(fmd[fmd$sigp == fmd$np & fmd$isconv == "TRUE", 1:10])
      dsigm <- dim(sigm)[1]
      if (dsigm > 0) {
        fmd1 <- sigm[(sigm$AIC == min(sigm$AIC)), ]
        selected <- paste(fmd1$fit)
        plotiMod <- as.numeric(paste(unique(fmd1$iMod)))
        if (plots==TRUE) { plotgdX(input1a, plotiMod) }
        fmd$selected <- selected
        fmd$Analyzed <- "yes"
        fmd$Group <- "included"
        kep <- c("name00", "lm", "Analyzed", "Group", "selected", "iMod",
                 "modelnames", "parameter", "Estimate", "Std..Error", "t.value",
                 "Pr...t..")
        fmd3 <- fmd[, c(kep)]
        fmd3$IDr <- k
        input1a$ploti <- as.numeric(paste(unique(fmd1$iMod)))
        fmd3
      } else {
        fmd$selected <- "not fit"
        fmd$Analyzed <- "yes"
        fmd$Group <- "excluded"
        kep <- c("name00", "lm", "Analyzed", "Group", "selected", "iMod",
                 "modelnames", "parameter", "Estimate", "Std..Error", "t.value",
                 "Pr...t..")
        fmd3 <- fmd[, c(kep)]
        fmd3$IDr <- k
        fmd3
      }
      colnames(fmd3)[1:2] <- c("name", "Nobs")
      return(fmd3)
    }

    ip <- inputprep(input)
    allinput <- ip$inputdata
    a <- allinput
    name <- unique(a$name)
    ln <- length(name)
    nm <- dim(foo)[1]

    if (ln > 0) {
      ID4 <- rep(1:ln, 1)
      b <- cbind(name, ID4)
      c <- merge(a, b, by = "name")
      conout <- do.call("rbind", sapply(1:ln, fid4, simplify = FALSE))
      conout
    } else {
      zaout <- xfit(fit = "NA", name00 = "name00")
      conout <- zaout
      conout
    }

    allconv <- conout
    pEst <- conout[, c(1, 4, 5, 7:12, 2, 13)]
    colnames(pEst) <- c("name", "type", "selected", "fit", "parameter", "Estimate",
                        "StdError", "t.value", "p.value", "N", "IDr")
    return(pEst)
  }

  # Function to create finalg/d/phi columns from results for des stats
  aggPE <- function(pEst) {
    # pEst is the output data frame from the checkconv2 function

    out3 <- data.frame(pEst)
    ncyc <- unique(out3[, c("name", "N", "type")])

    # parameter names from fit
    fits <- paste(unique(pEst$selected))
    foovars <- foo[(fits %in% foo$fit), ]
    gnames <- paste(foovars[foovars$gvar != "NA", c("gvar")])
    dnames <- paste(foovars[foovars$dvar != "NA", c("dvar")])
    phinames <- paste(foovars[foovars$pvar != "NA", c("pvar")])

    colkep <- c("Estimate", "StdError", "t.value", "p.value", "name", "selected")
    colsg <- c("g", "SEg", "tValueg", "pValueg", "name", "selectedFit")
    colsd <- c("d", "SEd", "tValued", "pValued", "name", "selectedFit")
    colsp <- c("phi", "SEphi", "tValuephi", "pValuephi", "name", "selectedFit")

    # gs
    out <- subset(out3, out3$parameter %in% gnames)
    if (dim(out)[1] > 0) {
      out1g <- unique(out[, colkep])
      colnames(out1g) <- colsg
      outg <- data.frame(out1g)
    } else {
      out1g <- matrix(c(NaN, NaN, NaN, NaN, NaN, NaN), nrow = 1, ncol = 6)
      colnames(out1g) <- colsg
      outg <- data.frame(out1g)
    }

    # ds
    out <- subset(out3, out3$parameter %in% dnames)
    if (dim(out)[1] > 0) {
      out1d <- unique(out[, colkep])
      colnames(out1d) <- colsd
      outd <- data.frame(out1d)
    } else {
      out1d <- matrix(c(NaN, NaN, NaN, NaN, NaN, NaN), nrow = 1, ncol = 6)
      colnames(out1d) <- colsd
      outd <- data.frame(out1d)
    }

    # phis
    out <- subset(out3, out3$parameter %in% phinames)
    if (dim(out)[1] > 0) {
      out1phi <- unique(out[, colkep])
      colnames(out1phi) <- colsp
      outphi <- data.frame(out1phi)
    } else {
      out1phi <- matrix(c(NaN, NaN, NaN, NaN, NaN, NaN), nrow = 1, ncol = 6)
      colnames(out1phi) <- colsp
      outphi <- data.frame(out1phi)
    }

    out6 <- merge(outg, outd, by = c("name", "selectedFit"), all = TRUE)
    out7 <- merge(out6, outphi, by = c("name", "selectedFit"), all = TRUE)
    out7$finalg <- ifelse((out7$selectedFit == "gdphi"), (out7$g * sqrt((1 -
                                                                           out7$phi))), out7$g)
    out7$finald <- out7$d
    out7$finalphi <- out7$phi
    out8 <- merge(out7, ncyc, by = "name")
    out9 <- out8[, c("name", "N", "type", "selectedFit", "finalg", "finald",
                     "finalphi")]
    colnames(out9)[5:7] <- c("g", "d", "phi")
    return(out9)
  }

  # Function to create mod table
  initf <- function() {
    IDmod <- rep(1:4, 1)
    fit <- c("gd", "dx", "gx", "gdphi")
    model <- c("f~exp(-d*time)+exp(g*(time))-1", "f~exp(-dx*time)", "f~exp(gx*time)",
               "f~(1-p)*exp(gt*time)+p*exp(-dt*time)")
    start <- c("list(g=0.00511,d=0.00511)", "list(dx=0.00511)", "list(gx=0.00511)",
               "list(p=.9,gt=.00511,dt=.00511)")
    cc <- c("blue", "purple3", "navy", "midnightblue")
    lb <- c("c(0,0)", "c(0)", "c(0)", "c(0,0,0)")
    ub <- c("c(1,1)", "c(1)", "c(1)", "c(1,1,1)")
    foo <- data.frame(cbind(IDmod, fit, model, start, cc, lb, ub))
    foo$IDmodel <- as.numeric(paste(foo$IDmod))
    foo$K <- c(2, 1, 1, 3)
    foo$gvar <- c("g", "NA", "gx", "gt")
    foo$dvar <- c("d", "dx", "NA", "dt")
    foo$pvar <- c("NA", "NA", "NA", "p")
    return(foo)
  }

  # Function to return result list
  generateresults <- function() {

    # Function to generate outlist1 for gdrate fx return
    genoutlist <- function(xx) {
      y <- data.frame(stats::aggregate(xx$name ~ xx$calcfinal, data = xx, length))
      colnames(y) <- c("Type", "N")
      y$Percentage <- round((y$N/sum(y$N)), digits = 2) * 100
      y$Group <- ifelse((y$Type %in% paste(foo$fit)), "included", "excluded")
      olt <- c(paste(foo$fit), "not fit")
      y$Analyzed <- ifelse((y$Type %in% olt), "yes", "no")
      y <- y[, c(4, 5, 1:3)]
      y <- data.frame(y[order(y$Group), ])
      rownames(y) <- NULL
      return(y)
    }

    # data prep
    ip <- inputprep(input)

    # where all input data is non-analyzed excluded cases
    if (ip$cont == 0) {
      resultsexc1only <- function() {
        ol <- data.frame(ip$excluded)
        outlist1 <- genoutlist(ol)
        ol2 <- ol[, c("name", "numcyc", "calcfinal")]
        colnames(ol2)[2:3] <- c("N", "selectedFit")
        OutputData <- cbind(ol2, g = NaN, d = NaN, phi = NaN)
        noa <- "no analyzable cases in input data"
        result <- list(allest = noa, results = OutputData, models = outlist1,
                       sumstats = noa)
        return(result)
      }

      result <- resultsexc1only()
      return(result)

    } else {

      # where analyzable cases
      resultsanalyzable <- function() {
        ip <- inputprep(input)
        allinput <- ip$inputdata

        # number of nonexcluded cases
        lnDSET <- length(unique(allinput$name))

        # model and optionally plot selected or return excluded
        allconv0 <- checkconv2()
        allconv <- allconv0[(allconv0$selected == allconv0$fit | allconv0$selected ==
                               "not fit"), ]
        kep <- c("name", "selected", "IDr", "N")
        res00 <- unique(allconv[, c(kep)])

        # create sets by inclusion and analysis status
        setex <- function() {
          all <- res00
          all$calcfinal <- all$selected
          allm <- all[, c("name", "calcfinal", "N")]

          # create empty dataset for if no cases
          naempty <- function() {
            name <- NaN
            dnone <- data.frame(cbind(name, calcfinal = NA, N = NaN))
            dnone$calcfinal <- as.character(dnone$calcfinal)
          }

          # included cases
          incl <- subset(allm, allm$calcfinal != "not fit")
          inc1 <- dim(incl)[1]
          if (inc1 > 0) {
            inc2 <- incl
            inc2$name <- as.numeric(as.character(inc2$name))
          } else {
            inc2 <- naempty()
          }

          # excluded cases analyzed
          notfit <- subset(allm, allm$calcfinal == "not fit")
          nf <- dim(notfit)[1]
          if (nf > 0) {
            nf2 <- notfit
            nf2$name <- as.numeric(as.character(nf2$name))
          } else {
            nf2 <- naempty()
          }

          # excluded cases not analyzed
          excl <- ip$excluded
          colnames(excl)[4] <- "N"
          exclm <- excl[stats::complete.cases(excl), c("name", "calcfinal", "N")]
          excl1 <- dim(exclm)[1]
          if (excl1 > 0) {
            excl2 <- exclm
          } else {
            excl2 <- naempty()
          }

          # merge and return list
          all1 <- data.frame(rbind(inc2, nf2, excl2))
          all2 <- all1[(all1$calcfinal != "NA"), ]
          row.names(all2) <- NULL
          combo <- list(inc1 = inc1, nf = nf, excl1 = excl1, all2 = all2)
          return(combo)
        }
        combos <- setex()
        inc1 <- combos$inc1
        excl1 <- combos$excl1
        nf <- combos$nf
        all2 <- data.frame(combos$all2)

        # summary by fit/exclusion reason
        outlist1 <- genoutlist(all2)

        # plots and output estimates where analyzed cases by inc exc status
        if (inc1 > 0) {
          # where included cases create table of estimates for selected fits of analyzed
          pEst <- allconv0[allconv0$selected == allconv0$fit, ]

          # generate dset with finalg/d/phi columns
          res <- aggPE(pEst = pEst)

          # descriptive stats by variable
          retAll <- function(vr, vc) {
            vals <- stats::na.omit(vr)
            lv <- length(vals)
            if (lv > 0) {
              N <- length(vals)
              Median <- round(stats::median(vals), digits = 6)
              IQR <- paste("(", round(stats::quantile(vals)[2], digits = 6), ", ",
                           round(stats::quantile(vals)[4], digits = 6), ")", sep = "")
              Mean <- round(mean(vals), digits = 6)
              SD <- round(stats::sd(vals), digits = 6)
              Parameter <- vc
              out1 <- data.frame(cbind(Parameter, N, Median, IQR, Mean, SD))
            } else {
              Parameter <- vc
              nv <- NaN
              out1 <- data.frame(cbind(Parameter, N = nv, Median = nv, IQR = nv,
                                       Mean = nv, SD = nv))
            }
            return(out1)
          }
          g <- retAll(res$g, "g")
          d <- retAll(res$d, "d")
          phi <- retAll(res$phi, "phi")
          outlist2 <- data.frame(rbind(g, d, phi))

          # merge results table with non analyzed excluded if existing
          if (excl1 > 0) {
            # excluded cases not analyzed
            excl <- ip$excluded
            excl$type <- "excluded"
            aid <- c("name", "numcyc", "type", "calcfinal")
            exclm <- excl[!is.na(excl$name), aid]
            colnames(exclm)[2:4] <- c("N", "type", "selectedFit")
            res$name <- as.numeric(as.character((res$name)))
            rc <- merge(res, exclm, by = c("name", "N", "type", "selectedFit"),
                        all = TRUE)
          } else {
            rc <- res
          }

          if (nf > 0) {
            notfit <- allconv0[(allconv0$selected == "not fit"), ]
            colnames(notfit)[3] <- c("selectedFit")
            kep <- c("name", "N", "type", "selectedFit")
            notfit <- unique(notfit[, c(kep)])
            rescalc <- merge(rc, notfit, by = c("name", "N", "type", "selectedFit"),
                             all = TRUE)
          } else {
            rescalc <- rc
          }

        } else {
          # where no included cases (ie nf must be GT 0)
          notfit <- allconv0[(allconv0$selected == "not fit"), ]
          colnames(notfit)[3] <- c("selectedFit")
          kep <- c("name", "N", "type", "selectedFit")
          notfit <- unique(notfit[, c(kep)])
          notfit$g <- NaN
          notfit$d <- NaN
          notfit$phi <- NaN
          rescalc <- notfit
          colnames(rescalc) <- c("name", "N", "type", "selectedFit", "g",
                                 "d", "phi")
          outlist2 <- "no estimates when zero included cases"
          allconv0 <- "no estimates when zero included cases"
        }  #end combos inc1 GT 0 else

        # list to output
        result <- list(allest = allconv0, results = rescalc, models = outlist1,
                       sumstats = outlist2)
        return(result)
      }

      result <- resultsanalyzable()
      return(result)
    }
  }

  foo <- initf()
  result <- generateresults()
  return(result)
}
