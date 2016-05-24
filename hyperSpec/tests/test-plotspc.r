context ("plotspc")

test_that("BARBITURATES", {
  spc <- do.call (collapse, barbiturates [1:3])

  plotspc (spc, col = matlab.dark.palette (3), stacked = TRUE, lines.args = list (type = "h"))

}) 


