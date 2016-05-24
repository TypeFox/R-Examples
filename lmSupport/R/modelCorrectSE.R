modelCorrectSE <- function(Model, Digits=3)
#2011-02-27:  added digits as a parameter
{
  newSEs = sqrt(diag(hccm(Model)))
  modelsum = summary(Model)
  thetable = modelsum$coefficients
  thetable[,2] = newSEs
  thetable[,3] = thetable[,1] / thetable[,2]
  thetable[,4] = 2*(pt(abs(thetable[,3]), df=modelsum$df[2], lower.tail=FALSE))

  cat('Uncorrected Tests of Coefficients\n\n')
  print(modelsum$coefficients, digits=Digits)
  cat('\nWhite (1980) Heteroscedasticity-corrected SEs and Tests\n\n')
  print(thetable, digits=Digits)
}