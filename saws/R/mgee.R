`mgee` <-
function (formula = formula(data), id = id, data = parent.frame(), 
    subset, na.action, R = NULL, b = NULL, tol = 0.001, maxiter = 25, 
    family = gaussian, corstr = "independence", Mv = 1, silent = TRUE, 
    contrasts = NULL, scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE) 
{
    call <- match.call()
    ## change call to gee, then evaluate
    call[[1]]<-as.name("gee")
    geeOutput<-eval(call)
    ## calculate the U and Omega values
    out<-geeUOmega(geeOutput)
    out
}
