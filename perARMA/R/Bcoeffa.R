Bcoeffa<-function (x, T, tau, missval, datastr, ...)
{
    Bcoeffa_full <- function(x, T, tau, missval, datastr, printflg,
        meth) {
        nout = floor((T + 2)/2)
        nx = length(x)
        pmean1 <- matrix(0, nx, 1)
        pmean <- matrix(0, T, 1)
        numtau = length(tau)
        Bkhat <- matrix(0, nout, numtau)
        Rhokhat <- matrix(0, nout, numtau)
        ahat <- matrix(0, T, numtau)
        n1 <- matrix(1, T, numtau)
        n2 <- matrix(1, T, numtau)
        fratio <- matrix(1, T, numtau)
        pvec <- matrix(1, T, numtau)
        nsamp <- matrix(0, 1, numtau)
        if (is.nan(missval)) {
            missisnan = 1
            imissx = x[(is.nan(x))]
        }
        else {
            missisnan = 0
            imissx = x[x == missval]
        }
        nmissx = length(imissx)
        for (i in 1:T) {
            index = seq(i, nx, T)
            z = x[index]
            if (missisnan) {
                igood = which(!is.nan(z))
                imiss = which(is.nan(z))
            }
            else {
                igood = which((z != missval))
                imiss = which((z == missval))
            }
            z = z[igood]
            pmean[i] = mean(z)
            x[index[imiss]] = pmean[i]
            pmean1[index] = pmean[i]
        }
        xd = x - pmean1
        xd[imissx] = 0
        for (k in 1:numtau) {
            lag = tau[k]
            alag = abs(lag)
            indt = seq(1, (nx - alag))
            indtt = seq((1 + alag), nx)
            corlen = length(indt)
            ncorper = floor(corlen/T)
            ncor = ncorper * T
            nsamp[k] = ncor
            indt = indt[seq(1, ncor)]
            indtt = indtt[seq(1, ncor)]
            xx = xd[indt] * xd[indtt]
            xxt = fft(xx)/ncor
            xxt2 = abs(xxt)^2
            kindex = seq(1, floor((ncor + 2 * ncorper)/2), ncorper)
            nk = length(kindex)
            Bkhat[, k] = xxt[kindex]
            if (lag < 0) {
                phase_correction = exp(i * 2 * pi * (kindex -
                  1) * lag/ncor)
                Bkhat[, k] = Bkhat[, k] * phase_correction
            }
            ahat[1, k] = Bkhat[1, k]
            kupper = floor((T - 1)/2)
            evenind = seq(2, 2 * kupper, 2)
            oddind = seq(3, (2 * kupper + 1), 2)
            ahat[evenind, k] = 2 * Re(Bkhat[(2:(kupper + 1)),
                k])
            ahat[oddind, k] = -2 * Im(Bkhat[(2:(kupper + 1)),
                k])
            if (T%/%2) {
                ahat[T, k] = Re(Bkhat[(T/2) + 1, k])
            }
            if (meth == 0) {
                kgood = setdiff(seq(1, floor(ncor/2)), kindex)
                n1[, k] = 1
                for (kk in 1:T) {
                  n2[kk, k] = 2 * length(kgood)
                  backav = sum(xxt2[kgood])/n2[kk, k]
                  num = (abs(ahat[kk, k])^2)/n1[kk, k]
                  fratio[kk, k] = num/backav
                  pvec[kk, k] = 1 - pf(fratio[kk, k], n1[kk,
                    k], n2[kk, k])
                }
            }
            else {
                if (meth >= ncorper) {
                  meth = 2 * floor((ncorper - 1)/2)
                }
                half = meth/2
                n[, k] = 1
                for (kk in 1:T) {
                  if (kk == 1) {
                    kgood = seq(2, (half + 1))
                  }
                  else {
                    if (T%%2 == 0 & kk == T) {
                      center = (T/2) * ncorper + 1
                      kgood = seq((center - half), (center -
                        1))
                    }
                    else {
                      center = (floor(kk/2)) * ncorper + 1
                      kgood = c(seq((center - half), (center -
                        1)), seq((center + 1), (center + half)))
                    }
                  }
                  n2[kk, k] = 2 * length(kgood)
                  backav = sum(xxt2[kgood])/n2[kk, k]
                  num = (abs(ahat[kk, k])^2)/n1[kk, k]
                  fratio[kk, k] = num/backav
                  pvec[kk, k] = 1 - pf(fratio[kk, k], n1[kk,
                    k], n2[kk, k])
                }
            }
        }
        if (printflg) {
            for (k in 1:numtau) {
                cat(paste("\n"))
                cat(paste("Bcoeffs in ab form for", datastr,"lag=", tau[k], "\n"))

          detail <- matrix(c(Re(ahat[, k]), Re(n1[, k]), Re(n2[,k]), Re(fratio[, k]),Re( pvec[, k])),ncol=5)
          colnames(detail) <- c( "ahat_k", "n1", "n2", "Fratio", "pv")
          row.names(detail)<-paste("k=",seq(0,T-1), sep="")
          print(detail)
            }
        }
        result = list(ahat = ahat, nsamp = nsamp, n1 = n1, n2 = n2,
            fratio = fratio, pvec = pvec)
        class(result) = "Bcoeffa"
        result
    }
    L <- modifyList(list(printflg = 1, meth = 0), list(x = x,
        T = T, tau = tau, missval = missval, datastr = datastr,
        ...))
    do.call(Bcoeffa_full, L)
}
