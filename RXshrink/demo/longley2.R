  require(RXshrink)
  # input revised Longley dataset of Hoerl(2000).
  data(longley2)
  # Specify form of regression model (linear here)...
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  # Fit of this model using 2-parameter Generalized Ridge Regression
  rxrobj <- RXridge(form, data=longley2)
  rxrobj
  names(rxrobj)
  plot(rxrobj)
  cat("\n Press ENTER for Least Angle Regression demo...")
  scan()
  # Fit of the above model with Least Angle Regression
  rxlobj <- RXlarlso(form, data=longley2)
  rxlobj
  names(rxlobj)
  plot(rxlobj)
  cat("\n Press ENTER for Least Angle fit to Uncorrelated Components...")
  scan()
  # Fit Least Angle Regression to Uncorrelated Components (closed form)...
  rxuobj <- RXuclars(form, data=longley2)
  rxuobj
  plot(rxuobj)