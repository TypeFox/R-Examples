persigest<-function (x, T, alpha, missval, datastr, ...)
{
    persigest_full <- function(x, T, alpha, missval, datastr,
        typeci, typepstd, pchci, pchpstd, colci, colpstd, pp) {
        nx = length(x)
        xind = seq(1, nx)
        nper = floor(nx/T)
        nxact = nper * T
        nrem = nx - nxact
        if (!is.nan(datastr)) {
            cat(paste("found ", nper, " periods of length ",
                T, " with remainder of ", nrem, "\n"))
        }
        if (is.nan(missval)) {
            missisnan = 1
            imissx = which(is.nan(x))
        }
        else {
            missisnan = 0
            imissx = which(x == missval)
        }
        nmissx = length(imissx)
        pmean1 <- matrix(0, nx, 1)
        pstd1 <- matrix(0, nx, 1)
        pmean <- matrix(0, T, 1)
        pstd <- matrix(0, T, 1)
        ny <- matrix(0, T, 1)
        psci <- matrix(0, T, 2)
        X <- matrix()

        vimiss<-c()
        for (i in 1:T) {
            index = seq(i, nx, T)
            z = x[index]
            if (missisnan) {
                igood = which(!is.nan(z))
                imiss = which(is.nan(z))
            }
            else {
                igood = which(z != missval)
                imiss = which(z == missval)
            }
            z = z[igood]
            ny[i] = length(z)
            pmean[i] = mean(z)
            pstd[i] = sd(z)
            x[index[imiss]] = pmean[i]
            pmean1[index] = pmean[i]
            pstd1[index] = pstd[i]
            chi20 = qchisq(c(alpha/2, 1 - alpha/2), ny[i] - 1)
            psci[i, 1] = pstd[i] * sqrt((ny[i] - 1)/chi20[2])
            psci[i, 2] = pstd[i] * sqrt((ny[i] - 1)/chi20[1])
            a <- matrix(1, length(z), 1)
            ai <- a * i
            b <- list(t(z), t(ai))

           vimiss[i]=length(imiss)
        }

         if (pp)
         {  detail <- matrix(c(ny,vimiss,pstd,psci[,1], psci[,2]),ncol=5)
            colnames(detail) <- c(" ngood", "nmiss"," pstd ", "lower", "upper")
            row.names(detail)<-paste("i=",seq(1,T), sep="")
            print(detail) }

        htest <- bartlett.test(b)
        pspv <- htest[3]
        xd = t(x) - pmean1
        if (pp) {
            matplot(psci, xlab = "seasons", ylab = "std", type = typeci,
                lwd = 1, lty = 1, col = colci, pch = pchci)
            points(pstd1, type = typepstd, lwd = 1, lty = 1,
                col = colpstd, pch = pchpstd)
            title(main = (paste("Periodic standard deviations: ",
                "No. periods =", nper, " alpha =", alpha)), sub = (paste("Bartlett p-value:",
                pspv)))
            legend("bottomright", c(expression(std), expression(confidence_intervals)),
                fill = c(colpstd, colci), ncol = 2, title = "legend")
        }
        xn = xd/pstd1
        result = list(xn = xn, pstd = pstd, psci = psci, pspv = pspv)
        class(result) = "persigest"
        result
    }
    L <- modifyList(list(typeci = "o", typepstd = "b", pchci = 10,
        pchpstd = 15, colci = "red", colpstd = "blue", pp = 1),
        list(x = x, T = T, alpha = alpha, missval = missval,
            datastr = datastr, ...))
    do.call(persigest_full, L)
}
