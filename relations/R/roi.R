ROI_options <-
local({
    options <- list()
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value))
            options[[option]]
        else
            options[[option]] <<- value
    }
})

ROI_solve <-
function(x, solver, control = NULL, ...)
    UseMethod("ROI_solve")

ROI_solve.MILP <-
function(x, solver, control = NULL, ...)
    solve_MILP(x, solver, control)

ROI_solve.MIQP <-
function(x, solver, control = NULL, ...)
    solve_MIQP(x, solver, control)
