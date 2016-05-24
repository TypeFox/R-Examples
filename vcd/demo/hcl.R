if(require("tcltk"))
{
    hue  <- tclVar("hue")
    chroma  <- tclVar("chroma")
    luminance  <- tclVar("luminance")
    fixup <- tclVar("fixup")
    hue <- tclVar(230)
    hue.sav <- 230
    chroma <- tclVar(55)
    chroma.sav <- 55
    luminance <- tclVar(75)
    luminance.sav <- 75
    fixup <- tclVar(FALSE)

    replot <- function(...) {
        hue.sav <- my.h <- as.numeric(tclvalue(hue))
        chroma.sav <- my.c <- as.numeric(tclvalue(chroma))
        luminance.sav <- my.l <- as.numeric(tclvalue(luminance))
	my.fixup <- as.logical(as.numeric(tclvalue(fixup)))
	barplot(1, col = hcl2hex(my.h, my.c, my.l, fixup = my.fixup), axes = FALSE)
    }

    replot.maybe <- function(...)
    {
        if(!((as.numeric(tclvalue(hue)) == hue.sav) &&
	     (as.numeric(tclvalue(chroma)) == chroma.sav) &&
	     (as.numeric(tclvalue(luminance)) == luminance.sav))) replot()
    }


    base <- tktoplevel()
    tkwm.title(base, "HCL Colors")

    spec.frm <- tkframe(base, borderwidth = 2)
    hue.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    chroma.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    luminance.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
    fixup.frm <- tkframe(spec.frm, relief = "groove", borderwidth = 2)

    tkpack(tklabel(hue.frm, text = "Hue"))
    tkpack(tkscale(hue.frm, command = replot.maybe, from = 0, to = 360,
                   showvalue = TRUE, variable = hue,
                   resolution = 1, orient = "horiz"))

    tkpack(tklabel(chroma.frm, text = "Chroma"))
    tkpack(tkscale(chroma.frm, command = replot.maybe, from = 0, to = 100,
                   showvalue = TRUE, variable = chroma,
                   resolution = 5, orient = "horiz"))

    tkpack(tklabel(luminance.frm, text = "Luminance"))
    tkpack(tkscale(luminance.frm, command = replot.maybe, from = 0, to = 100,
                   showvalue = TRUE, variable = luminance,
                   resolution = 5, orient = "horiz"))
    tkpack(tklabel(fixup.frm, text="Fixup"))
    for (i in c("TRUE", "FALSE") ) {
        tmp <- tkradiobutton(fixup.frm, command = replot,
                             text = i, value = as.logical(i), variable = fixup)
        tkpack(tmp, anchor="w")
    }

    tkpack(hue.frm, chroma.frm, luminance.frm, fixup.frm, fill="x")

    ## Bottom frame on base:
    q.but <- tkbutton(base, text = "Quit",
                      command = function() tkdestroy(base))

    tkpack(spec.frm, q.but)

    replot()
}
