
library(playwith)

## A tool to draw a line and label it with its equation.
## The annotations are persistent and may be redrawn in other contexts.
## The actions can be reversed with "Undo" menu item or <Ctrl>Z.

## This "callback" is run when the action is activated
## (from toolbar, menu, or keyboard shortcut).
lineEq_handler <- function(widget, playState) {
    ## draw line at drag locations
    lnInfo <- playLineInput(playState)
    if (is.null(lnInfo)) return()
    if (isTRUE(lnInfo$is.click)) return()
    ## annotation call to draw a line
    lnExpr <- with(lnInfo$coords,
                   call("panel.segments",
                        x[1], y[1], x[2], y[2]))
    playAnnotate(playState, lnExpr, space = lnInfo$space)
    ## draw equation at click location
    eqInfo <- playPointInput(playState, prompt =
                          paste("Click to place equation,",
                                "Right-click to cancel."))
    if (is.null(eqInfo)) return()
    grad <- with(lnInfo$coords, diff(y) / diff(x))
    icept <- with(lnInfo$coords, y[1] - grad * x[1])
    ## create the equation as an expression
    eqn <- substitute(expression(y == a * x + b),
                      list(a = signif(grad, 3),
                           b = signif(icept, 3)))
    ## annotation call to draw text
    eqExpr <- with(eqInfo$coords,
                   call("panel.usertext", x[1], y[1], eqn))
    playAnnotate(playState, eqExpr, space = eqInfo$space)
}

lineEqTool <- list("LineEq", "gtk-indent", "Line + Eqn",
                   callback = lineEq_handler)

playwith(plot(1:10), tools = list(lineEqTool))
## Click on tool, drag a line, then click to place equation.
