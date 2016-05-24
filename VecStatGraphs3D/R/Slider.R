Slider <- function (value) 
{
    tt <- tktoplevel()
    tkwm.title(tt, "Size: only 5 drawing")
    SliderValue <- tclVar(as.character(value))
    SliderValueLabel <- tklabel(tt, text = as.character(tclvalue(SliderValue)))
    tkgrid(tklabel(tt, text = "Size: "), SliderValueLabel)
    tkconfigure(SliderValueLabel, textvariable = SliderValue)
    slider <- tkscale(tt, from = 0, to = 10, showvalue = TRUE, 
        variable = SliderValue, resolution = 0.1)
    value = 0
    done = 0
    aux_function_slider <- function() {
        val <- as.double(tclvalue(SliderValue))
        tclvalue(done) <- 1
        return(val)
    }
    OK.but <- tkbutton(tt, text = "OK", command = function() value <<- aux_function_slider())
    tkgrid(slider)
    tkgrid(OK.but)
    tkwait.variable(done)
    tkfocus(tt)
    tkdestroy(tt)
    return(value)
}
