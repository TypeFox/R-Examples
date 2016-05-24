if(require("tcltk"))
{
    hue  <- tclVar("hue")
    luminance  <- tclVar("luminance")
    saturation  <- tclVar("saturation")
    hue <- tclVar(0)
    hue.sav <- 0
    luminance <- tclVar(0.5)
    luminance.sav <- 0.5
    saturation <- tclVar(1)
    saturation.sav <- 1

    replot <- function(...) {
        hue.sav <- my.h <- as.numeric(tclvalue(hue))
        saturation.sav <- my.s <- as.numeric(tclvalue(saturation))
        luminance.sav <- my.l <- as.numeric(tclvalue(luminance))
	barplot(1, col = hls(my.h, my.l, my.s), axes = FALSE)
    }

    replot.maybe <- function(...)
    {
        if(!((as.numeric(tclvalue(hue)) == hue.sav) &&
	     (as.numeric(tclvalue(saturation)) == saturation.sav) &&
	     (as.numeric(tclvalue(luminance)) == luminance.sav))) replot()
    }


    base <- tktoplevel()
    tkwm.title(base, "HLS Colors")

    spec.frm <- tkframe(base, borderwidth = 2)
    hue.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    saturation.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    luminance.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)

    tkpack(tklabel(hue.frm, text = "Hue"))
    tkpack(tkscale(hue.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = hue,
                   resolution = 0.01, orient = "horiz"))

    tkpack(tklabel(luminance.frm, text = "Luminance"))
    tkpack(tkscale(luminance.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = luminance,
                   resolution = 0.01, orient = "horiz"))

    tkpack(tklabel(saturation.frm, text = "Saturation"))
    tkpack(tkscale(saturation.frm, command = replot.maybe, from = 0, to = 1,
                   showvalue = TRUE, variable = saturation,
                   resolution = 0.01, orient = "horiz"))

    tkpack(hue.frm, luminance.frm, saturation.frm, fill="x")

    ## Bottom frame on base:
    q.but <- tkbutton(base, text = "Quit",
                      command = function() tkdestroy(base))

    tkpack(spec.frm, q.but)

    replot()
}
