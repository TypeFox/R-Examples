coef.bma <-
function (object, exact = FALSE, order.by.pip = TRUE, include.constant = FALSE, 
    incl.possign = TRUE, std.coefs = FALSE, condi.coef = FALSE, 
    ...) 
{
    estimates.bma(object, exact = exact, order.by.pip = order.by.pip, 
        include.constant = include.constant, incl.possign = incl.possign, 
        std.coefs = std.coefs, condi.coef = condi.coef)
}
