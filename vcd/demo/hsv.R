if(require("tcltk"))
{
    hue  <- tclVar("hue")
    saturation  <- tclVar("saturation")
    value  <- tclVar("value")
    hue <- tclVar(0)
    hue.sav <- 0
    saturation <- tclVar(1)
    saturation.sav <- 1
    value <- tclVar(1)
    value.sav <- 1

    replot <- function(...) {
        hue.sav <- my.h <- as.numeric(tclvalue(hue))
        saturation.sav <- my.s <- as.numeric(tclvalue(saturation))
        value.sav <- my.v <- as.numeric(tclvalue(value))
	barplot(1, col = hsv(my.h, my.s, my.v), axes = FALSE)
    }

    replot.maybe <- function(...)
    {
        if(!((as.numeric(tclvalue(hue)) == hue.sav) &&
	     (as.numeric(tclvalue(saturation)) == saturation.sav) &&
	     (as.numeric(tclvalue(value)) == value.sav))) replot()
    }


    base <- tktoplevel()
    tkwm.title(base, "HSV Colors")

    spec.frm <- tkframe(base, borderwidth = 2)
    hue.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    saturation.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    value.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)

    tkpack(tklabel(hue.frm, text = "Hue"))
    tkpack(tkscale(hue.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = hue,
                   resolution = 0.01, orient = "horiz"))

    tkpack(tklabel(saturation.frm, text = "Saturation"))
    tkpack(tkscale(saturation.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = saturation,
                   resolution = 0.01, orient = "horiz"))

    tkpack(tklabel(value.frm, text = "Value"))
    tkpack(tkscale(value.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = value,
                   resolution = 0.01, orient = "horiz"))

    tkpack(hue.frm, saturation.frm, value.frm, fill="x")

    ## Bottom frame on base:
    q.but <- tkbutton(base, text = "Quit",
                      command = function() tkdestroy(base))

    tkpack(spec.frm, q.but)

    replot()
}
