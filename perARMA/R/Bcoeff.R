Bcoeff<-function (x, T, tau, missval, datastr, ...)
{
    Bcoeff_full <- function(x, T, tau, missval, datastr, printflg,
        meth) {
        nout = floor((T + 2)/2)
        nx = length(x)
        pmean1 <- matrix(0, nx, 1)
        pmean <- matrix(0, T, 1)
        numtau = length(tau)
        Bkhat <- matrix(0, nout, numtau)
        Rhokhat <- matrix(0, nout, numtau)
        n1 <- matrix(1, nout, numtau)
        n2 <- matrix(1, nout, numtau)
        fratio <- matrix(1, nout, numtau)
        pvec <- matrix(1, nout, numtau)
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
        for (k in 1:numtau) {
            lag = tau[k]
            alag = abs(lag)
            indt = seq(1, (nx - alag))
            indtt = seq((1 + alag), nx)
            corlen = length(indt)
            ncorper = floor(corlen/T)
            ncor = ncorper * T
            nsamp[k] = ncor
            indt = indt[1:ncor]
            indtt = indtt[1:ncor]
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
            if (meth == 0) {
                kgood = setdiff(seq(1, floor(ncor/2)), kindex)
                for (kk in 1:length(kindex)) {
                  if (kk == 1) {
                    n1[kk, k] = 1
                  }
                  else {
                    n1[kk, k] = 2
                  }
                  n2[kk, k] = 2 * length(kgood)
                  backav = sum(xxt2[kgood])/n2[kk, k]
                  num = (abs(Bkhat[kk, k])^2)/n1[kk, k]
                  fratio[kk, k] = num/backav
                  pvec[kk, k] = 1 - pf(fratio[kk, k], n1[kk,
                    k], n2[kk, k])
                }
            }
            else {
                if (meth >= ncorper) {
                  meth = 2 * floor((ncorper - 1)/2)
                }
                half = floor(meth/2)
                for (kk in 1:length(kindex)) {
                  if (kk == 1) {
                    kgood = seq(2, (half + 1))
                    n1[kk, k] = 1
                  }
                  else {
                    if (T%%2 == 0 & kk == length(kindex)) {
                      center = kk * ncorper
                      kgood = seq((center - half), (center -
                        1))
                      n1[kk, k] = 1
                    }
                    else {
                      center = kk * ncorper
                      kgood = c(seq((center - half), (center -
                        1)), seq((center + 1), (center + half)))
                      n1[kk, k] = 2
                    }
                  }
                  n2[kk, k] = 2 * length(kgood)
                  backav = sum(xxt2[kgood])/n2[kk, k]
                  num = (abs(Bkhat[kk, k])^2)/n1[kk, k]
                  fratio[kk, k] = num/backav
                  pvec[kk, k] = 1 - pf(fratio[kk, k], n1[kk,
                    k], n2[kk, k])
                }
            }
        }
        if (printflg) {
            for (k in 1:numtau) {
                 cat(paste("\n"))
                 cat(paste("Bcoeffs for", datastr, "lag=", tau[k],"\n"))

                detail <- matrix(c(Re(Bkhat[,k]), Im(Bkhat[,k]), n1[,k], n2[,k], fratio[,k],pvec[,k]),ncol=6)
                colnames(detail) <- c( "reB_k"," imB_k ", "n1", "n2", "Fratio", "pv")
                row.names(detail)<-paste("k=",seq(0,(length(kindex)-1)), sep="")
                print(detail)
           }
        }
        result = list(Bkhat = Bkhat, nsamp = nsamp, n1 = n1,
            n2 = n2, fratio = fratio, pvec = pvec)
        class(result) = "Bcoeff"
        result
    }
    L <- modifyList(list(printflg = 1, meth = 0), list(x = x,
        T = T, tau = tau, missval = missval, datastr = datastr,
        ...))
    do.call(Bcoeff_full, L)
}


