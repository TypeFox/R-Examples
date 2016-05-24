`sumres` <-
function(x) {
  sr <- summary(residuals(x))
  srm <- mean(residuals(x))
  if (abs(srm) < 1e-10){
    sr <- sr[c(1:3,5:6)]
  }
  sr
}

