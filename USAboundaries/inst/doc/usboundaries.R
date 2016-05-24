## ------------------------------------------------------------------------
library(USAboundaries)
library(sp) # for plotting

states <- us_states()
plot(states)

## ------------------------------------------------------------------------
states_1790 <- us_states("1790-07-04")
plot(states_1790)

## ------------------------------------------------------------------------
counties <- us_counties(states = c("South Carolina", "North Carolina"))
plot(counties)

