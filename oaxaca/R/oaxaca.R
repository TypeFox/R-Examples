# oaxaca - Blinder-Oaxaca Decomposition
# Author: Marek Hlavac

oaxaca <-
function(formula, data, weights = NULL, R = 100, reg.fun = lm,  ...) {  
  cl <- match.call()
    
  return(.oaxaca.wrap(formula=formula, data=data, weights=weights, R=R, reg.fun=reg.fun, cl=cl, ...))
}

summary.oaxaca <- function(object, ...) {
  return(object)
}

plot.oaxaca <- function(x, decomposition = "threefold", type = "variables",
                        weight = NULL, unexplained.split = FALSE,
                        variables = NULL, components = NULL,
                        component.left = FALSE,
                        component.labels = NULL,
                        variable.labels = NULL,
                        ci = TRUE, ci.level = 0.95, 
                        title = "", xlab = "", ylab = "", 
                        bar.color = NULL, ...) {
  cl <- match.call()
  
  return(.plot.oaxaca(x=x, decomposition=decomposition, type=type,
                      w=weight, unexplained.split=unexplained.split,
                      variables=variables, components=components,
                      component.left=component.left,
                      component.labels=component.labels, variable.labels=variable.labels, 
                      ci=ci, ci.level=ci.level, 
                      title=title, xlab=xlab, ylab=ylab, 
                      bar.color=bar.color, cl=cl, ...))
}

