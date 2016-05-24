diag.panel.hist <-
function# define scatterplot diagonal formatting
### internal function for sisus
(x
### internal variable
, ...
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  usr = par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$counts; y = y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="steelblue3", ...)
  ### internal variable
}
