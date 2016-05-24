hwmi <- function(yTref, Tref, yTemp, Temp) {
    
    t3daymax <- NULL
    timet3daymax <- NULL
                                                                                              
    
    Tref <- cut366thDay(yTref, (yTref + 31), Tref)
    
    yendTemp <- (trunc(length(Temp) / 365) - 1)
    
    Temp <- cut366thDay(yTemp, (yTemp + yendTemp), Temp)
    
    threshold <- mvThreshold(yTref, (yTref+31), (yTref+1), (yTref+30), Tref, 90, 15)

    for (ij in seq((365 + 1), length(Tref), 365)) {
    
        tnday <- consecDayMax(2, 364, Tref[ ij:(ij + 364) ], 3)
        t3daymax <- c(t3daymax, tnday[2])
        timet3daymax <- c(timet3daymax, tnday[1])
    }
    
    
    bw <- bw.SJ(t3daymax)
    PDF <- density(t3daymax,bw)
    pdfx <- PDF$x
    pdfy <- PDF$y
    
    nyears <- as.integer(length(Temp) / 365)
    hwmi <- hwmiFun(nyears, Temp, threshold, pdfx, pdfy)

    return(list(hwmi = hwmi, thr = threshold, pdfx = pdfx, pdfy = pdfy))

} # end of 'hwmi' function.



cut366thDay <- function(t1, t2, Ti) {

    timey <- c(t1:t2)
    ltime <- length(timey)

    t4 <- timey / 4
    t4int <- as.integer(t4)

    tdiff <- (t4 - t4int)
    tindex <- which(tdiff == 0)

    ltindex <- length(tindex)

    if (ltindex > 0) {

        dayindex <- NULL
        dayi <- NULL

        k <- 1
        j <- 0

        for (i in tindex) {

            index <- (365 * i + k)
            dayindex <- c(dayindex, index)

            if (k == 1) dayi <- c(dayi, c(1:(dayindex[k] - 1)))
            else if (k > 1) dayi <- c(dayi, c((dayindex[k - 1] + 1):(dayindex[k] - 1)))

            k <- k + 1

        } # end of for 'i' loop.

        if (tindex[ ltindex ] < ltime) dayi <- c(dayi, c((dayindex[k - 1] + 1):length(Ti)))

        T365 <- Ti[ dayi ]

    } # end of if 'ltindex > 0' stmts.

    if (ltindex == 0) T365 <- Ti

    if ((as.integer(length(Ti) / 365) - length(Ti) / 365) == 0)  T365 <- Ti

    return(T365)

} # end of 'cut366thDay' function.



mvThreshold <- function(t1, t2, t1ref, t2ref, Ti, thresh.lev, nday) {

    timey <- c(t1:t2)
    lt <- length(timey) * 365
    refper.ind <- which(timey == t1ref)
    refper.ind1 <- ((refper.ind - 1) * 365 + 1)
    refper.ind2 <- ((refper.ind + 29) * 365)

    T2 <- Ti[(refper.ind1 - nday):(refper.ind2 + nday)]

    timeyref <- c(t1ref:t2ref)

    Thresh <- c(1:365)

    for (i in (1 + nday):(365 + nday)) {

        # si parte da 3 e si finisce a 367 perche' abbiamo aggiunto 2 giorni
        # a destra e a sinistra nella serie di temperature T2 quindi T2 ha
        # lunghezza 365*30+4 gg

        i30 <- seq(i, (29 * 365 + i), 365)

        t5days<-T2[i30]
        for (j in 1:nday) t5days <- c(t5days, T2[i30 - j], T2[i30 + j])

        a <- quantile(t5days, probs = seq(0, 1, 0.1), na.rm = TRUE, names = FALSE, type = 8)
        b <- quantile(t5days, probs = seq(0, 1, 0.01), na.rm = TRUE, names = FALSE, type = 8)

        if (thresh.lev == 10) Thresh[(i - nday)] <- a[2]

        if (thresh.lev == 90) Thresh[(i - nday)] <- a[10]

        if (thresh.lev == 5) Thresh[(i - nday)] <- b[6]

        if (thresh.lev == 75) Thresh[(i - nday)] <- b[76]

        if (thresh.lev == 50) Thresh[(i - nday)] <- b[51]

        if (thresh.lev == 25) Thresh[(i - nday)] <- b[26]

        if (thresh.lev == 95) Thresh[(i - nday)] <- b[96]

        if (thresh.lev == 2) Thresh[(i - nday)] <- b[3]

        if (thresh.lev == 98) Thresh[(i - nday)] <- b[99]

    } # end of for 'i' loop.

    return(Thresh)

}



consecDayMax <- function(day1, day2, Ti, nday) {

    lday <- (day2 - day1 + 1)
    step <- trunc((nday - 1) / 2)
    tnday <- NULL

    for (i in day1:day2) tnday <- c(tnday, sum(Ti[ (i - step):(i + step) ]))
    tmax <- max(tnday)
    timetmax <- (which(tnday == tmax) + day1 - 1)

    tndaymax <- c(timetmax, tmax)
    return(tndaymax)

} # end of 'consecDayMax' function.



hwmiFun <- function(nyears, Ti, thrsummer, pdfx, pdfy) {

    hw.scale <- matrix(-111, nyears, 3)
    hw.scale[1:nyears,1:3] <- 0

    f <- approxfun(pdfx, pdfy, yleft = 0, yright = 0)
    k <- 0

    for (year in 1:nyears) {

        Ts <- Ti[(k * 365 + 1):((k + 1) * 365)]
        indTs <- which(is.na(Ts))
        Ts[ indTs ] <- -999.9

        hw <- hwyear(Ts,thrsummer[1:365],3)
        phw3 <- c(1:length(hw[,1]))
        phw3[ 1:length(hw[,1]) ] <- 0

        if (hw[1,1]>0) {

            for (id in 1:length(hw[,1])) {

                hwd <- hw[id, 1]
                hwt <- hw[id, 5]

                if (hwd < 364)  hw3d <- ((round(hwd / 3 + 0.4) * 3))

                if (hwd > 363)  {

                    hw3d <- 363
                    hwd <- 363

                }

                if ((hwd / 3 - as.integer(hwd / 3)) == 0) thw <- Ts[ hwt:(hwt + hwd - 1) ]

                if (0 < (hwd / 3 - as.integer(hwd / 3)) && (hwd / 3 - as.integer(hwd / 3)) < 0.5) {

                    if (hwt > 1 && ((365 - hwt) - hwd) > 2) thw <- Ts[ (hwt - 1):(hwt + hwd) ]
                    if (hwt > 2 && ((365 - hwt) - hwd) < 3) thw <- Ts[ (hwt - 2):(hwt + hwd - 1) ]

                    if (hwt == 2 && 0 < ((365 - hwt) - hwd) && ((365 - hwt) - hwd) < 3) thw <- Ts[ (hwt - 1):(hwt + hwd) ]
                    if (hwt == 2 && ((365 - hwt) - hwd) == 0) thw <- Ts[ 1:120 ]

                    if (hwt == 1) thw <- Ts[ hwt:(hwt + hwd + 1) ]

                }

                if ((hwd / 3 - as.integer(hwd / 3)) > 0.5) {

                    if (hwt > 1 && ((365 - hwt) - hwd) > 2) {

                        thw <- c(Ts[ hwt:(hwt + hwd - 1) ], max(Ts[ (hwt - 1) ], Ts[ hwt:(hwt + hwd) ]))

                    }

                    if (hwt > 2 && ((365 - hwt) - hwd) < 3) thw <- c(Ts[ hwt:(hwt + hwd - 1) ], max(Ts[ (hwt - 1) ], Ts[ hwt - 2 ]))
                    if (hwt == 2 && 0 < ((365 - hwt) - hwd) && ((365 - hwt) - hwd) < 3) {

                        thw <- c(Ts[ hwt:(hwt + hwd - 1) ], max(Ts[ (hwt - 1) ], Ts[ hwt + hwd ]))

                    }
                    if (hwt == 2 && ((365 - hwt) - hwd) == 0) thw <- Ts[ 1:363 ]

                    if (hwt == 1) thw <- c(Ts[ hwt:(hwt + hwd - 1) ], max(Ts[ (hwt + hwd) ], Ts[ (hwt + hwd + 1) ]))

                } # end of if 'hwd' stmts.

                #abbiamo eliminato il sort di thw per calcolare il fattore di probabilita' su tre giorni consecutivi appartenenti alla hwave.
                #sort(thw,decreasing=TRUE)->thw

                for (ik in seq(1, hw3d, 3)) {

                    thw3 <- sum(thw[ ik:(ik + 2) ])
                    fval <- integrate(f, lower = -Inf, upper = thw3, stop.on.error = FALSE)
                    fvalue <- fval$val
                    phw3[ id ] <- (phw3[ id ] + fvalue)

                } # end of for 'ik' loop.

            } # end of for 'id' loop.

        } # end of if 'id in 1:length(hw[,1])' stmts.

        if (max(phw3)>0) {

            indtime <- which(phw3 == max(phw3))

            if (length(indtime) > 1) indtime <- indtime[1]

            hw.time <- hw[indtime, 5]

            hw.duration <- hw[indtime, 1]

            hw.scale[year, 1] <- max(phw3)

            ###duration in position 2
            hw.scale[year, 2] <- hw.duration

            ###duration in position 3
            hw.scale[year, 3] <- hw.time

        } # end of if max phw3 stmts.

        k <- k + 1

    } # end of for 'year' loop.

    return(hw.scale)

} # end of 'hwmiFun' function.



hwyear <- function(Ti, threshold, nd) {

    index <- Ti > threshold

    temp.ind <- index * Ti

    time <- c(1:length(Ti))

    v <- time * index
    z1 <- which(v == 0)

    HWF<-0
    HWtime <- HWD <- HWI <- HWInd <- HWIndT <- HW <- NULL

    if (sum(index) >= 122) {

        HWD <- sum(index)
        HWI <- sum(Ti - threshold)
        HWsort <- sort((Ti - threshold), decreasing = TRUE)
        HWInd <- sum(HWsort[ 1:nd ])
        HWsortT <- sort(Ti, decreasing = TRUE)
        HWIndT <- sum(HWsortT[ 1:nd ])
        HWtime <- 1

    } # end of if 'sum of index >= 122' stmts.

    if (length(z1) == 1) {

        if (z1[1]>nd) {

            HWD0 <- sum(index[ 1:(z1[1] - 1) ])

            HWF <- HWF + 1
            HWtime[HWF] <- 1

            HWD[ HWF ] <- HWD0

            HWI[ HWF ] <- sum(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ])
            HWsort <- sort((Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1)] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ]), decreasing = TRUE)

            HWInd[ HWF ] <- sum(HWsort[ 1:nd ])
            HWsortT <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1)], decreasing = TRUE)
            HWIndT[ HWF ] <- sum(HWsortT[ 1:nd ])

        } # end of if 'z1[1] == 1' stmts.

        if (max(z1) < length(time) && sum(index[ max(z1):length(time) ]) >= nd) {

            HWF <- HWF + 1

            HWD[ HWF ] <- sum(index[ max(z1):length(time) ])

            HWtime[ HWF ] <- max(z1) + 1

            HWI[ HWF ] <- sum(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ])

            HWsort <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)

            HWInd[ HWF ] <- sum(HWsort[ 1:nd ])
            HWsortT <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)
            HWIndT[ HWF ] <- sum(HWsortT[ 1:nd ])

        }

    } # end of if 'length of z1 = 1' stmts.

    if (length(z1)>1) {

        if (z1[1]>nd) {

            HWD0 <- sum(index[ 1:(z1[1] - 1) ])
            HWF <- HWF + 1

            HWtime[ HWF ] <- 1
            HWD[ HWF ] <- HWD0
            HWI[ HWF ] <- sum(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ])

            HWsort <- sort((Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ]), decreasing = TRUE)
            HWInd[ HWF ] <- sum(HWsort[ 1:nd ])

            HWsortT <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)
            HWIndT[ HWF ] <-sum(HWsortT[ 1:nd ])

        } # end of if 'z1[1] > nd' stmts.

        for (i in 1:(length(z1)-1)) {

            HWD0 <- sum(index[ (z1[i] + 1):(z1[i + 1] - 1) ])

            if (HWD0>=nd) {

                HWF <- HWF + 1
                HWtime[ HWF ] <- (z1[i] + 1)

                HWD[ HWF ] <- HWD0

                HWI[ HWF ] <- sum(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ])
                HWsort <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1)], decreasing = TRUE)

                HWInd[ HWF ] <- sum(HWsort[ 1:nd ])
                HWsortT <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)
                HWIndT[ HWF ] <- sum(HWsortT[ 1:nd ])

            } # end of if 'HWD0 >= nd' stmts.

        } # end of for 'i' loop.

        if (max(z1) < length(time) && sum(index[ max(z1):length(time) ]) >= nd) {

            HWF <-HWF+1
            HWD[ HWF ] <- sum(index[ max(z1):length(time) ])

            HWtime[ HWF ] <- max(z1) + 1

            HWI[ HWF ] <- sum(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ])
            HWsort <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ] - threshold[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)

            HWInd[ HWF ] <- sum(HWsort[ 1:nd ])
            HWsortT <- sort(Ti[ (HWtime[HWF]):(HWtime[HWF] + HWD[HWF] - 1) ], decreasing = TRUE)
            HWIndT[ HWF ] <- sum(HWsortT[ 1:nd ])

        }

    } # end of if length of z > 1 stmts.

    HW <- c(HWD, HWI, HWInd, HWIndT, HWtime)

    if (length(HWD)>0) dim(HW)<-c(length(HWD),5)
    else if (length(HWD)==0) {

        HW<-c(-1111,-1111,-1111,-1111,-1111)
        dim(HW)<-c(1,5)

    } # end of if else 'length of HWD' stmts.

    return(HW)

} # end of 'hwyear' function.

