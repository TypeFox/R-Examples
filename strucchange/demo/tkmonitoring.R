if(require("tcltk"))
{
    data(UKDriverDeaths)
    seatbelt <- log10(UKDriverDeaths)
    seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
    colnames(seatbelt) <- c("y", "ylag1", "ylag12")
    seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
    data(GermanM1)
    data(durab)

    data  <- tclVar("M1")
    type <- tclVar("OLS-CUSUM")
    border <- tclVar(1)
    h  <- tclVar(0.5)
    h.sav <- 0.5

    replot <- function(...) {
        h.sav <- hh <- as.numeric(tclvalue(h))
        tp <- tclvalue(type)
	bd <- tclvalue(border)
        dt <- tclvalue(data)
	  switch(dt,
          "UK Seatbelt" = {
	    seat.sub <<- window(seatbelt, start = c(1975,11), end = c(1983,1))
            seat.efp <- efp(y ~ ylag1 + ylag12, data = seat.sub, type = tp, h = hh)
            if(bd > 0 & tp %in% c("OLS-CUSUM", "RE"))
	      bd <- newborder <- function(k) 1.5778*k/seat.efp$nobs
	    else
	      bd <- NULL
	    seat.mefp <- mefp(seat.efp, period = 2, border = bd)
	    seat.sub <<- window(seatbelt, start = c(1975, 11))
            seat.mon <- monitor(seat.mefp, verbose = FALSE)
	    plot(seat.mon)
	  },
	  "M1" = {
            M1 <<- historyM1
	    m1.efp <- efp(dm ~ dy2 + dR + dR1 + dp + ecm.res + season, type = tp, h = hh, data = M1)
            if(bd > 0 & tp %in% c("OLS-CUSUM", "RE"))
	      bd <- newborder <- function(k) 1.5778*k/m1.efp$nobs
	    else
	      bd <- NULL
	    m1.mefp <- mefp(m1.efp, period = 2, border = bd)
	    M1 <<- GermanM1
	    m1.mon <- monitor(m1.mefp, verbose = FALSE)
	    plot(m1.mon)
          },
          "US Durables" = {
	    Durab <<- window(durab, start=1964, end = c(1979, 12))
            durab.efp <- efp(y ~ lag, type = tp, h = hh, data = Durab)
            if(bd > 0 & tp %in% c("OLS-CUSUM", "RE"))
	      bd <- newborder <- function(k) 1.5778*k/durab.efp$nobs
	    else
	      bd <- NULL
            durab.mefp <- mefp(durab.efp, period=2, border = bd)
            Durab <<- window(durab, start=1964)
            durab.mon <- monitor(durab.mefp, verbose = FALSE)
            plot(durab.mon)
	  })
    }


    base <- tktoplevel()
    tkwm.title(base, "Monitoring")

    spec.frm <- tkframe(base, borderwidth = 2)
    left.frm <- tkframe(spec.frm)
    right.frm <- tkframe(spec.frm)

    ## Left frame:
    frame1 <- tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame1, text="Process type"))
    for (i in c("OLS-CUSUM", "OLS-MOSUM", "RE", "ME") ) {
        tmp <- tkradiobutton(frame1, command = replot,
                             text = i, value = i, variable = type)
        tkpack(tmp, anchor="w")
    }

    frame4 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel (frame4, text = "border type"))
    tmp <- tkradiobutton(frame4, command = replot,
                         text = "Chu et al.", value = 0, variable = border)
    tkpack(tmp, anchor="w")
    tmp <- tkradiobutton(frame4, command = replot,
                         text = "Zeileis et al.", value = 1, variable = border)
    tkpack(tmp, anchor="w")


    ## Two right frames:
    frame2 <-tkframe(right.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame2, text="Data set"))
    for (i in c("UK Seatbelt", "M1", "US Durables") ) {
        tmp <- tkradiobutton(frame2, command = replot,
                             text = i, value = i, variable = data)
        tkpack(tmp, anchor="w")
    }

    frame3 <- tkframe(right.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel (frame3, text = "Bandwidth h"))
    for (i in c(0.25, 0.5, 1) ) {
        tmp <- tkradiobutton(frame3, command = replot,
                             text = i, value = i, variable = h)
        tkpack(tmp, anchor="w")
    }

    tkpack(frame1, frame4, fill="x")
    tkpack(frame2, frame3, fill="x")
    tkpack(left.frm, right.frm, side = "left", anchor = "n")

    ## Bottom frame on base:
    q.but <- tkbutton(base, text = "Quit",
                      command = function() tkdestroy(base))

    tkpack(spec.frm, q.but)

    replot()
}
