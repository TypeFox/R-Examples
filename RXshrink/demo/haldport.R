  require(RXshrink)
  # read in Portland cement hardening data of Hald(1952).
  data(haldport)
  form <- heat~p3ca+p3cs+p4caf+p2cs
  rxrobj <- RXridge(form, data=haldport)
  rxrobj
  plot(rxrobj)
  cat("\n Press ENTER for True MSE with known parameters...")
  scan()
  # define true parameter values.
  trugam <- matrix(c(.8,.0,.3,.5),4,1)
  trusig <- 0.2
  # create true shrinkage MSE risk scenario.
  trumse <- RXtrisk(form, data=haldport, trugam, trusig, Q=-5)
  trumse
  plot(trumse)
  cat("\n Press ENTER to simulate True Squared Error Loss...")
  scan()
  trusim <- RXtsimu(form, data=haldport, trugam, trusig, Q=-5)
  trusim
  plot(trusim)
  cat("\n Press ENTER for simulated ridge shrinkage as if parameters were unknown...")
  scan()
  haldpsim <- haldport
  haldpsim[,5] <- trusim$ydat[,1]
  rxsobj <- RXridge(form, data=haldpsim)
  rxsobj
  plot(rxsobj)