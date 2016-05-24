initializeDialogDoE <-
function(title = "", offset = 0, preventCrisp = FALSE)
{
    ## works in principle, but moves a little bit to the right bottom at each refresh
    if ((!preventCrisp) && getRcmdr("crisp.dialogs"))
        tclServiceMode(on = FALSE)
    putRcmdr("topdes2", tktoplevel(borderwidth = 10))
    tkwm.title(topdes2, title)
    position <- if (is.SciViews())
        -1
    else if (exists("topleft2xy", where="RcmdrEnv")) topleft2xy
    else commanderPosition()
    position <- if (any(position < 0))
        "-50+50"
    else paste("+", paste(offset + position, collapse = "+"),
        sep = "")
    tkwm.geometry(topdes2, position)
}
