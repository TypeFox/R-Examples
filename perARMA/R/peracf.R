peracf<-function (x, T, tau, missval, datastr, ...)
{
    peracf_full <- function(x, T, tau, missval, datastr, prttaus,
        plottaus, cialpha, typeci, typerho, pchci, pchrho, colci,
        colrho) {
        nx = length(x)
        pmean1 <- matrix(0, nx, 1)
        pmean <- matrix(0, T, 1)
        numper = floor(nx/T)
        numtau = length(tau)
        B <- matrix(0, T, numtau)
        Rho <- matrix(0, T, numtau)
        nsamp <- matrix(0, T, numtau)
        if (is.nan(missval)) {
            missisnan = 1
            imissx = x[(is.nan(x))]
        }
        else {
            missisnan = 0
            imissx = x[x == missval]
            x[imissx] = NaN
        }
        for (t in 1:T) {
            index <- seq(t, nx, T)
            z = x[index]
            zomit = na.omit(z)
            pmean[t] = mean(zomit)
            pmean1[index] = pmean[t]
        }
        xd = x - t(pmean1)
        for (t in 1:T) {
            indt = seq(t, nx, T)
            for (k in 1:numtau) {
                lag = tau[k]
                indtt = seq(t + lag, nx, T)
                if (t + lag) {
                  mfirst = 1
                }
                else {
                  mfirst = 2 - fix((t + lag)/T)
                }
                indt = indt[mfirst:length(indt)]
                indtt = indtt[mfirst:length(indtt)]
                lent = length(indt)
                lentt = length(indtt)
                len_to_use = min(lent, lentt)
                indt = indt[1:len_to_use]
                indtt = indtt[1:len_to_use]
                xx = xd[indtt] * xd[indt]
                xxomit = na.omit(xx)
                Bcorr = sum(xxomit)/len_to_use
                nsamp[t, k] = sum(!is.nan(xx))
                xd_indt_omit = na.omit(xd[indt] * xd[indt])
                sigt = sqrt(sum(xd_indt_omit)/len_to_use)
                xd_indtt_omit = na.omit(xd[indtt] * xd[indtt])
                sigtt = sqrt(sum(xd_indtt_omit)/len_to_use)
                B[t, k] = Bcorr
                Rho[t, k] = Bcorr/(sigt * sigtt)
            }
        }
        prttaus = intersect(prttaus, tau)
        plottaus = intersect(plottaus, tau)
        saveequalpv = NULL
        savezeropv = NULL
        if (!length(union(prttaus, plottaus)) == 0) {
            for (k in 1:numtau) {
                lag = tau[k]
                if (is.element(lag, union(prttaus, plottaus)) &
                  !lag == 0) {
                  lower <- matrix(NaN, T, 1)
                  upper <- matrix(NaN, T, 1)
                  inotnan = !is.nan(B[, k])
                  ibigenuf = as.logical(nsamp[, k] > 3)
                  iuse = as.logical(inotnan * ibigenuf)
                  rhoci <- matrix(0, length(iuse), 2)
                  for (z in 1:length(iuse)) {
                    rho <- rhoci(Rho[z, k], nsamp[z, k], cialpha)
                    lower = as.matrix(rho$lower)
                    upper = as.matrix(rho$upper)
                    rhoci[z, 1] = lower
                    rhoci[z, 2] = upper
                  }
                  rct <- rho.constant.test(Rho[iuse, k], nsamp[iuse,
                    k])
                  equalpv <- rct$pv
                  equalchi2 <- rct$chi2
                  ngood <- rct$ngood
                  rzt <- rho.zero.test(Rho[iuse, k], nsamp[iuse,
                    k])
                  zeropv <- rzt$pv
                  mzscore <- rzt$mzscore
                  ngood <- rzt$ngood

                  cat(paste("\n"))


                 cat(paste("lag=", tau[k],"\n"))
                 detail <- matrix(c(equalchi2,mzscore,equalpv,zeropv),ncol=2)
                 colnames(detail) <- c("test", "pv")
                 rownames(detail)<-c("rho(t+lag,t)=rho(lag), chi2 =", "rho(t+lag,t)=0, mzscore = ")
                 print(detail)
                 saveequalpv=c(saveequalpv,equalpv)
                 savezeropv=c(savezeropv,zeropv)


                  if (is.element(lag, prttaus)) {

                 detail <- matrix(c(Rho[, k],rhoci[, 1],rhoci[, 2],nsamp[, k]),ncol=4)
                 colnames(detail) <- c("rho(t,lag)", "lower", "upper", "nsamp")
                 row.names(detail)<-paste("t=",seq(1,T), sep="")
                 print(detail)
                    }


                  if (is.element(lag, plottaus)) {
                    dev.set(which = 1)
                    matplot(seq(1, T), rhoci, xlab = "time t",
                      ylab = "rho(t,lag)", type = typeci, lwd = 1,
                      lty = 4, col = colci, pch = pchci, new = TRUE)
                    points(seq(1, T), Rho[, k], type = typerho,
                      lwd = 2, lty = 4, col = colrho, pch = pchrho,
                      new = TRUE)
                    legend("bottom", c("coefficients", "confidence intervals"),
                      fill = c(colrho, colci), ncol = 2, title = paste("Tests results: equalpv= ",
                        equalpv, "; zeropv= ", zeropv))
                    title(main = paste(" Correlation coefficients rho(t,lag) for:",
                      datastr, " for lag=", tau[k]))
                  }
             }
}

        }

        nbf = length(saveequalpv)
        if (nbf) {
            nmsaveequalpv = nbf * min(saveequalpv)
            if (nmsaveequalpv >= 1) {
                nmsaveequalpv = 1
            }
            nmsavezeropv = nbf * min(savezeropv)
            if (nmsavezeropv >= 1) {
                nmsavezeropv = 1
            }
            cat(paste("\n"))
            cat(paste("least equalpv bonferroni corrected for",
                nbf, " lags tried:", nmsaveequalpv, "\n"))
            cat(paste("least equalpv bonferroni corrected for",
                nbf, " lags tried:", nmsavezeropv, "\n"))
        }
        result = list(B = B, Rho = Rho, nsamp = nsamp)
        class(result) = "peracf"
        result
    }
    L <- modifyList(list(prttaus = seq(1, T/2), plottaus = seq(1,
        T/2), cialpha = 0.05, typeci = "b", typerho = "b", pchci = 10,
        pchrho = 15, colci = "blue", colrho = "red"), list(x = x,
        T = T, tau = tau, missval = missval, datastr = datastr,
        ...))
    do.call(peracf_full, L)
}
