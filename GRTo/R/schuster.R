schuster <-
function (finame, title = "Schuster's diagram", color = c("black", 
    "red", "orange", "green", "cyan", "navy"), hd = FALSE, colidye = NULL, 
    colidmo = NULL, colidda = NULL, colidho = 1, colidmi = 2, 
    colidma = 3, colidz = NULL, utccor = 0, dayt1 = NULL, dayt2 = NULL, 
    ysel1 = NULL, ysel2 = NULL, mosel1 = NULL, mosel2 = NULL, 
    dasel1 = NULL, dasel2 = NULL, magsel1 = NULL, magsel2 = NULL, 
    zsel1 = NULL, zsel2 = NULL, weekday1 = c("mo", "tu", "we", 
        "th", "fr", "sa", "su"), weekday2 = c("mo", "tu", "we", 
        "th", "fr", "sa", "su")) 
{
    color <- match.arg(color)
    if (is.null(dayt1)) 
        dayt1 = 0
    if (is.null(dayt2)) 
        dayt2 = 24
    dt1 = dayt1
    dt2 = dayt2
    loc = read.table(finame, header = hd)
    if (!is.null(colidye)) {
        y = colidye
        year = loc[, y]
        minyear = min(year, na.rm = TRUE)
        maxyear = max(year, na.rm = TRUE)
        cat("Year range", minyear, maxyear, "\n")
        if (minyear < 1900 | maxyear > 2050) 
            stop("Year data is out of range (1900 - 2050).")
        year1 = minyear
        year2 = maxyear
        if (!is.null(ysel1)) {
            year1 = ysel1
            cat(year1, " ")
            if (!is.null(ysel2)) 
                year2 = ysel2
            cat(year2)
            cat("\n")
        }
    }
    if (!is.null(colidmo)) {
        mo = colidmo
        month = loc[, mo]
        minmo = min(month, na.rm = TRUE)
        maxmo = max(month, na.rm = TRUE)
        if (is.null(colidye)) 
            cat("Month range", minmo, maxmo, "\n")
        if (minmo < 1 | maxmo > 12) 
            stop("Month data is out of range (1 - 12).")
        month1 = minmo
        month2 = maxmo
        if (!is.null(mosel1)) {
            month1 = mosel1
            cat(month1, " ")
            if (!is.null(mosel2)) 
                month2 = mosel2
            cat(month2)
            cat("\n")
        }
    }
    if (!is.null(colidda)) {
        da = colidda
        day = loc[, da]
        minda = min(day, na.rm = TRUE)
        maxda = max(day, na.rm = TRUE)
        if (is.null(colidmo)) 
            cat("Day range", minda, maxda, "\n")
        if (minda < 1 | maxda > 31) 
            stop("Month data is out of range (1 - 31).")
        day1 = minda
        day2 = maxda
        if (!is.null(dasel1)) {
            day1 = dasel1
            cat(day1, " ")
            if (!is.null(dasel2)) 
                day2 = dasel2
            cat(day2)
            cat("\n")
        }
    }
    dayofweek = NULL
    if (!is.null(colidda) & !is.null(colidmo) & !is.null(colidye)) {
        da = colidda
        day = loc[, da]
        mo = colidmo
        month = loc[, mo]
        y = colidye
        year = loc[, y]
        formaday = ISOdate(year, month, day)
        dayofweek = format(formaday, "%w")
        dayofweek = ifelse(dayofweek == 0, 7, dayofweek)
        if (!is.null(weekday1)) 
            wk1 = switch(weekday1, mo = 1, tu = 2, we = 3, th = 4, 
                fr = 5, sa = 6, su = 7)
        if (!is.null(weekday2)) 
            wk2 = switch(weekday2, mo = 1, tu = 2, we = 3, th = 4, 
                fr = 5, sa = 6, su = 7)
    }
    if (!is.null(colidma)) {
        m = colidma
        mag = loc[, m]
        minmag = min(mag, na.rm = TRUE)
        maxmag = max(mag, na.rm = TRUE)
        cat("Magnitude range", minmag, maxmag, "\n")
        if (minmag < -3 | maxmag > 9.5) 
            stop("Magnitude data is out of range (-3 - 9.5).")
        mag1 = minmag
        mag2 = maxmag
        if (!is.null(magsel1)) {
            mag1 = magsel1
            cat(mag1, " ")
            if (!is.null(magsel2)) 
                mag2 = magsel2
            cat(mag2)
            cat("\n")
        }
    }
    if (!is.null(colidz)) {
        d = colidz
        z = loc[, d]
        minz = min(z, na.rm = TRUE)
        maxz = max(z, na.rm = TRUE)
        cat("Depth range", minz, maxz, "\n")
        z1 = minz
        z2 = maxz
        if (!is.null(zsel1)) {
            z1 = zsel1
            cat(z1, " ")
            if (!is.null(zsel2)) 
                z2 = zsel2
            cat(z2)
            cat("\n")
        }
    }
    sel1 = TRUE
    if (!is.null(colidye)) 
        sel1 <- (year >= year1 & year <= year2)
    sel2 = TRUE
    if (!is.null(colidmo)) 
        sel2 <- (month >= month1 & month <= month2)
    sel3 = TRUE
    if (!is.null(colidda)) 
        sel3 <- (day >= day1 & day <= day2)
    sel4 = TRUE
    if (!is.null(colidz)) 
        sel4 <- (z >= z1 & z <= z2)
    sel5 = TRUE
    if (!is.null(colidma)) 
        sel5 <- (mag >= mag1 & mag <= mag2)
    sel6 = TRUE
    if (!is.null(dayofweek) & !is.null(weekday1) & !is.null(weekday2)) 
        sel6 <- (dayofweek >= wk1 & dayofweek <= wk2)
    if (!is.null(colidz)) 
        zdef <- z[which(sel1 & sel2 & sel3 & sel4 & sel5 & sel6)]
    if (!is.null(colidma)) 
        magdef <- mag[which(sel1 & sel2 & sel3 & sel4 & sel5 & 
            sel6)]
    if (!is.null(colidye)) 
        yeardef <- year[which(sel1 & sel2 & sel3 & sel4 & sel5 & 
            sel6)]
    if (!is.null(colidmo)) 
        modef <- month[which(sel1 & sel2 & sel3 & sel4 & sel5 & 
            sel6)]
    if (!is.null(colidmo)) 
        daydef <- day[which(sel1 & sel2 & sel3 & sel4 & sel5 & 
            sel6)]
    h = colidho
    htmp = loc[, h]
    total = length(htmp)
    if (min(htmp) < 0 | max(htmp) >= 24) 
        stop("Hour data is out of range (0-23)")
    hh <- htmp[which(sel1 & sel2 & sel3 & sel4 & sel5 & sel6)]
    minu = colidmi
    mtmp = loc[, minu]
    if (min(mtmp) < 0 | max(mtmp) >= 60) 
        stop("Minute data is out of range (0-59)")
    mm <- mtmp[which(sel1 & sel2 & sel3 & sel4 & sel5 & sel6)]
    heure = hh + utccor + mm/60
    heure[which(heure > 24)] = heure[which(heure > 24)] - 24
    heure[which(heure < 0)] = heure[which(heure < 0)] + 24
    htab = data.frame(heure)
    write.table(htab, file = "hour.res", col.names = FALSE, row.names = FALSE)
    heure = heure - dt1
    heure[which(heure > 24)] = heure[which(heure > 24)] - 24
    heure[which(heure < 0)] = heure[which(heure < 0)] + 24
    tint = dt2 - dt1
    if (tint < 0) 
        tint = tint + 24
    angtmp = heure * 360/tint
    print(length(angtmp))
    angle = angtmp[which(angtmp <= 360)]
    x = cos(angle * atan(1)/45)
    y = sin(angle * atan(1)/45)
    coordx = cumsum(x)
    coordy = cumsum(y)
    Rxy = sqrt(coordx * coordx + coordy * coordy)
    Rprob = max(Rxy)
    N = length(x)
    R = sqrt(-N * log(0.01))
    theta = seq(0, 360, 1)
    xcirc = R * sin(theta * atan(1)/45)
    ycirc = R * cos(theta * atan(1)/45)
    coordcirc = data.frame(xcirc, ycirc)
    figname = paste(title, "_schu.png", sep = "")
    png(filename = figname, width = 2400, height = 2715, res = 360)
    cat("\nSchuster's plot in", figname, "\n\n")
    par(mfrow = c(1, 1), pty = "s")
    rscale = max(range(xcirc, x, y))
    limits <- c(-rscale, rscale) * 1.08
    plot(0, 0, ylim = limits, xlim = limits, axes = FALSE, xlab = "", 
        ylab = "", type = "n", asp = 1)
    points(xcirc, ycirc, type = "l", lty = "dotted", lwd = 2)
    xcirc2 = sqrt(-N * log(0.05)) * sin(theta * atan(1)/45)
    ycirc2 = sqrt(-N * log(0.05)) * cos(theta * atan(1)/45)
    points(xcirc2, ycirc2, type = "l", lty = "dotted", lwd = 2)
    walkcolor = color
    points(c(0, coordx[1]), c(0, coordy[1]), type = "l", col = walkcolor, 
        lwd = 2)
    points(coordx, coordy, type = "l", col = walkcolor, lwd = 2)
    ytitl = 1.1 * sqrt(-N * log(0.01))
    titl = paste(title, ", Local time", sep = "")
    text(0, ytitl, titl, cex = 2)
    xtext = 1.1 * sqrt(-N * log(0.01)) * sin(310 * atan(1)/45)
    ytext = 1.1 * sqrt(-N * log(0.01)) * cos(310 * atan(1)/45)
    text(xtext, ytext, label = "99%", cex = 1.3)
    xtext = 0.9 * sqrt(-N * log(0.05)) * sin(290 * atan(1)/45)
    ytext = 0.9 * sqrt(-N * log(0.05)) * cos(290 * atan(1)/45)
    text(xtext, ytext, label = "95%", cex = 1.3, pos = 2, offset = 0.1)
    points(c(0, 0), c(rscale * 0.25, -rscale * 0.25), type = "l")
    points(c(rscale * 0.25, -rscale * 0.3), c(0, 0), type = "l")
    lab1 = paste(dt1, ":00", sep = "")
    text(rscale * 0.25, 0, label = lab1, cex = 0.9, pos = 4, 
        offset = 0.5)
    h = dt1 + 0.25 * tint
    if (h > 24) 
        h = h - 24
    if (h < 0) 
        h = h + 24
    mi = (h - trunc(h)) * 60
    if (!mi) 
        lab2 = paste(trunc(h), ":00", sep = "")
    if (mi) 
        lab2 = paste(trunc(h), ":", mi, sep = "")
    text(0, +rscale * 0.25, label = lab2, cex = 0.9, pos = 3, 
        offset = 0.5)
    h = dt1 + 0.5 * tint
    if (h > 24) 
        h = h - 24
    if (h < 0) 
        h = h + 24
    mi = (h - trunc(h)) * 60
    if (!mi) 
        lab3 = paste(trunc(h), ":00", sep = "")
    if (mi) 
        lab3 = paste(trunc(h), ":", mi, sep = "")
    text(-rscale * 0.25, 0, label = lab3, cex = 0.9, pos = 2, 
        offset = 1)
    h = dt1 + 0.75 * tint
    if (h > 24) 
        h = h - 24
    if (h < 0) 
        h = h + 24
    mi = (h - trunc(h)) * 60
    if (!mi) 
        lab4 = paste(trunc(h), ":00", sep = "")
    if (mi) 
        lab4 = paste(trunc(h), ":", mi, sep = "")
    text(0, -rscale * 0.25, label = lab4, cex = 0.9, pos = 1, 
        offset = 0.5)
    par(xpd = NA)
    leg1 = paste("N=", N, " (total ", total, ")", sep = "")
    if (!is.null(colidye)) 
        leg2 = paste("Year range ", min(yeardef, na.rm = TRUE), 
            "-", max(yeardef, na.rm = TRUE), " (", year1, "-", 
            year2, ")", sep = "")
    if (!is.null(colidmo)) 
        leg3 = paste("Month range ", min(modef, na.rm = TRUE), 
            "-", max(modef, na.rm = TRUE), " (", month1, "-", 
            month2, ")", sep = "")
    if (!is.null(colidda)) 
        leg4 = paste("Day range ", min(daydef, na.rm = TRUE), 
            "-", max(daydef, na.rm = TRUE), " (", day1, "-", 
            day2, ")", sep = "")
    if (!is.null(colidz)) 
        leg5 = paste("Depth range ", min(zdef, na.rm = TRUE), 
            "-", max(zdef, na.rm = TRUE), " km (", z1, "-", z2, 
            ")", sep = "")
    if (!is.null(colidma)) 
        leg6 = paste("Magnitude range ", min(magdef, na.rm = TRUE), 
            "-", max(magdef, na.rm = TRUE), " (", mag1, "-", 
            mag2, ")", sep = "")
    if (!is.null(dayofweek) & !is.null(weekday1) & !is.null(weekday2)) 
        leg7 = paste("Day of week range : ", weekday1, "-", weekday2, 
            sep = "")
    coefini = 1.1
    text(-rscale * 1.3, -rscale * coefini, label = leg1, cex = 1.1, 
        adj = 0)
    dpos = 0.1
    coefleg = coefini + dpos
    if (!is.null(colidye)) {
        text(-rscale * 1.3, -rscale * coefleg, label = leg2, 
            cex = 1.1, adj = 0)
        coefleg = coefleg + dpos
    }
    if (!is.null(colidmo)) {
        text(-rscale * 1.3, -rscale * coefleg, label = leg3, 
            cex = 1.1, adj = 0)
        coefleg = coefleg + dpos
    }
    if (!is.null(colidda)) {
        text(-rscale * 1.3, -rscale * coefleg, label = leg4, 
            cex = 1.1, adj = 0)
        coefleg = coefleg + dpos
    }
    if (!is.null(colidz)) {
        text(rscale * 0.05, -rscale * coefini, label = leg5, 
            cex = 1.1, adj = 0)
        coefleg = coefini + dpos
    }
    if (!is.null(colidma)) {
        text(rscale * 0.05, -rscale * coefleg, label = leg6, 
            cex = 1.1, adj = 0)
        coefleg = coefleg + dpos
    }
    if (!is.null(dayofweek) & !is.null(weekday1) & !is.null(weekday2)) {
        text(rscale * 0.05, -rscale * coefleg, label = leg7, 
            cex = 1.1, adj = 0)
        coefleg = coefleg + dpos
    }
    prob = exp(-Rprob * Rprob/N)
    leg8 = paste("prob = ", signif(prob, 4), sep = "")
    text(rscale * 0.05, -rscale * coefleg, label = leg8, cex = 1.1, 
        adj = 0)
}
