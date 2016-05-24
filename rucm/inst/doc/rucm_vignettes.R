## ----install, echo = T, message = F, eval = T----------------------------
#install.packages("rucm")
library(rucm)

## ----modelNile, echo = TRUE, cache=FALSE---------------------------------
modelNile <- ucm(formula = Nile~0, data = Nile, level = TRUE)
modelNile #Printing method for class ucm
plot(Nile, ylab = "Flow of Nile")
lines(modelNile$s.level, col = "blue")
legend("topright", legend = c("Observed flow","S_level"), col = c("black","blue"), lty = 1)

## ----forecast, echo = TRUE, eval = TRUE----------------------------------
modelNile <- ucm(formula = Nile~0, data = Nile, level = TRUE, slope = TRUE)
predict(modelNile$model, n.ahead = 12) # Forecasting

