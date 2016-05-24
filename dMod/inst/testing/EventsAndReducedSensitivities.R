library(deSolve)
library(cOde)
library(dMod)

f <- NULL
f <- addReaction("2*A", "B", "k1*A*Egf", f)
f <- addReaction("A", "B", "k3*A", f)
f <- addReaction("B", "A", "k2*B*u", f)
f <- addReaction("", "Egf", "0", f)

forcings <- "u"
observables <- c(tot = "s1*(5*A + B) + off", rat = "A/B")
model <- generateModel(f, forcings, fixed = NULL, jacobian="inz.lsodes")
events <- data.frame(var = "Egf", time =2, value = 1, method = "add")
forcdata <- data.frame(name = "u", time = 0:5, value = 1)

x <- Xs(model$func, model$extended, forcings = forcdata, events = events, optionsSens = list(atol=1e-6, rtol=1e-6, verbose=TRUE, lrw = 4000))
g <- Y(observables, f)
y <- function(times, pars, ...) g(x(times, pars, ...), pars)

times <- seq(0, 5, len=101)
pars <- c(k1 = 1, k2 = 1, k3=.1, A = 0, B = 1, s1 = 1, off = 1)

innerpars <- getSymbols(c(f, observables), exclude=forcings); names(innerpars) <- innerpars
innerpars <- replaceSymbols(names(innerpars), paste("exp(", names(innerpars), ")"), innerpars)
innerpars["Egf"] <- "1"

p <- P(innerpars)

out <- wide2long(y(times, p(pars)))

timesD <- c(0, 1, 2, 3, 4, 5)
data <- subset(out, time %in% timesD)
data$value <- data$value + rnorm(length(data$value), 0, 0.01)
data <- cbind(data, sigma = 0.01)
data <- data[c(1, 3, 4, 5, 6, 7, 10),]

print(plotCombined(list(a = y(times, p(pars))), list(a = data)))

obj <- function(pouter, fixed=NULL, deriv=TRUE) {
  
  prediction <- y(0:5, p(pouter, fixed), deriv = deriv)
  out <- wrss(res(data, prediction))

  out
  
}

