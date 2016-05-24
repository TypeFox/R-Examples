if(require("tcltk"))
{
    data(Nile)
    data(UKDriverDeaths)
    seatbelt <- log10(UKDriverDeaths)
    seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
    colnames(seatbelt) <- c("y", "ylag1", "ylag12")
    seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
    data(GermanM1)
    data(Grossarl)
    data(durab)

    data  <- tclVar("Nile")
    type <- tclVar("OLS-MOSUM")
    h  <- tclVar(0.15)
    h.sav <- 0.15

    replot <- function(...) {
        h.sav <- hh <- as.numeric(tclvalue(h))
        tp <- tclvalue(type)
        dt <- tclvalue(data)
        if(tp == "data") {
	  switch(dt,
          "Nile" = plot(Nile, ylab = "annual flow", main = "Measurements of the annual flow of the Nile at Ashwan"),
          "UK Seatbelt" = plot(seatbelt[,"y"], ylab = expression(log[10](casualties)), main = "UK seatbelt data"),
          "M1" = plot(GermanM1[,"m"], ylab = "money demand", main = "German M1 money demand"),
          "Grossarl" = plot(Grossarl$fraction, ylab = "fraction of illegitimate births", main = "Illegitimate births in Grossarl"),
          "US Durables" = plot(durab[,"y"], ylab = "productivity in the manufacturing/durables sector", main = "US labor productivity")
          )}
	else if(tp == "F statistics") {
	  switch(dt,
          "Nile" = plot(Fstats(Nile ~ 1), main = "F statistics"),
          "UK Seatbelt" = plot(Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = 0.1), main = "F statistics"),
          "M1" = plot(Fstats(dm ~ dy2 + dR + dR1 + dp + ecm.res + season, data = GermanM1, from = 0.2, to = 0.9), main = "F statistics"),
          "Grossarl" = plot(Fstats(fraction ~ politics, data = Grossarl), main = "F statistics"),
          "US Durables" = plot(Fstats(y ~ lag, data = durab), main = "F statistics")
	  )}
	else {
	  switch(dt,
          "Nile" = plot(efp(Nile ~ 1, type = tp, h = hh)),
          "UK Seatbelt" = plot(efp(y ~ ylag1 + ylag12, type = tp, h = hh, data = seatbelt)),
          "M1" = plot(efp(dm ~ dy2 + dR + dR1 + dp + ecm.res + season, type = tp, h = hh, data = GermanM1)),
          "Grossarl" = plot(efp(fraction ~ politics, type = tp, h = hh, data = Grossarl)),
          "US Durables" = plot(efp(y ~ lag, type = tp, h = hh, data = durab))
	  )}
    }

    replot.maybe <- function(...)
    {
        if((tclvalue(type) %in% c("Rec-MOSUM", "OLS-MOSUM", "ME")) & (as.numeric(tclvalue(h)) != h.sav)) replot()
    }


    base <- tktoplevel()
    tkwm.title(base, "Testing")

    spec.frm <- tkframe(base, borderwidth = 2)
    left.frm <- tkframe(spec.frm)
    right.frm <- tkframe(spec.frm)

    ## Left frame:
    frame1 <- tkframe(left.frm, relief="groove", borderwidth=2)
    tkpack(tklabel(frame1, text="Plot type"))
    for (i in c("data", "Rec-CUSUM", "Rec-MOSUM", "OLS-CUSUM", "OLS-MOSUM",
                 "RE", "ME", "F statistics") ) {
        tmp <- tkradiobutton(frame1, command = replot,
                             text = i, value = i, variable = type)
        tkpack(tmp, anchor="w")
    }

    ## Two right frames:
    frame2 <-tkframe(right.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame2, text="Data set"))
    for (i in c("Nile", "UK Seatbelt", "M1", "Grossarl", "US Durables") ) {
        tmp <- tkradiobutton(frame2, command = replot,
                             text = i, value = i, variable = data)
        tkpack(tmp, anchor="w")
    }

    frame3 <- tkframe(right.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel (frame3, text = "Bandwidth h"))
    tkpack(tkscale(frame3, command = replot.maybe, from = 0.05, to = 0.95,
                   showvalue = TRUE, variable = h,
                   resolution = 0.005, orient = "horiz"))

    tkpack(frame1, fill="x")
    tkpack(frame2, frame3, fill="x")
    tkpack(left.frm, right.frm, side = "left", anchor = "n")

    ## Bottom frame on base:
    q.but <- tkbutton(base, text = "Quit",
                      command = function() tkdestroy(base))

    tkpack(spec.frm, q.but)

    replot()
}
