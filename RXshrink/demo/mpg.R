  require(RXshrink)
  # input Gasoline mileage data of Hocking(1976).
  data(mpg)
  # Specify linear regression model with only four predictors...
  form <- mpg~cylnds+cubins+hpower+weight
  # Fit of this model using 2-parameter Generalized Ridge Regression
  rxrobj <- RXridge(form, data=mpg)
  rxrobj
  plot(rxrobj)
  cat("\n Press ENTER to make True Risk Calculations...")
  scan()
  # Define true parameter values...
  trugam <- matrix(c(-.5,-.1,.1,-.6),4,1)
  trusig <- 0.4
  # Create true shrinkage MSE risk scenario.
  trumse <- RXtrisk(form, data=mpg, trugam, trusig, Q=-1, steps=4)
  trumse
  plot(trumse)
  cat("\n Press ENTER to simulate True Squared Error Loss...")
  scan()
  trusel <- RXtsimu(form, data=mpg, trugam, trusig, Q=-1, steps=4)
  trusel
  plot(trusel)