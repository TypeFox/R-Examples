##load data
data(mesa.model)
data(est.mesa.model)

##find regression parameters using GLS
x.reg <- predict(mesa.model, est.mesa.model, only.pars = TRUE)
str(x.reg$pars)

\dontrun{
  ##compute predictions at all locations, including beta-fields
  pred.mesa.model <- predict(mesa.model, est.mesa.model,
                             pred.var=TRUE)
}
##Let's load precomputed results instead.
data(pred.mesa.model)

##study results
print(pred.mesa.model)
