#' @demo Concentration of transfer spendings in the German
#'   Bundesliga; analysis from the manuscript "Measuring
#'   Concentration in Data with an Exogenous Order" by Abedieh,
#'   Eugster, and Augustin (2011).

library("SportsAnalytics")



### Data: ############################################################

data("BundesligaFinalStandings")
data("BundesligaTransferSums")


spendings <- function(season) {
  s <- subset(BundesligaTransferSums, Season == season)
  structure(s$Spendings, names = as.character(s$Team))
}

standings <- function(season) {
  s <- subset(BundesligaFinalStandings, Season == season)
  as.character(s$Team)[s$Position]
}



### Step-by-step for season 2009/2010: ###############################

s <- "2009-2010"

sp <- spendings(s)         # Transfer spendings
st <- standings(s)         # Final league standing

o <- match(st, names(sp))  # Order of spendings according to standings

sp[o]                      # Spendings sorted according to standings



### Concentration ratios:

concentration_ratio(sp, 3)     # Classical concentration ratio
concentration_ratio(sp, 3, o)  # Concentration ratio with exogene order



### Concentration curves:

cr_sp <- function(g, ...) {
  concentration_ratio(sp, g, ...)
}

CR <- sapply(seq(along = sp), cr_sp)
OR <- sapply(seq(along = sp), cr_sp, o)



## Plot (Figure~1):
par(mfrow = c(2, 1))
plot(CR, type = "b", xlab = "Number of teams (by spendings)", ylab = expression(CR[g]))
plot(OR, type = "b", xlab = "Number of teams (by standings)", ylab = expression(OR[g]))



### Concentration indices:

herfindahl(sp)
rosenbluth(sp)
exogeny(sp, ex = o)



### Analysis for all seasons 1992/1993--2009/2010: ###################

### Concentration curves:

cc <- function(season, type = c("CR", "OR")) {
  type <- match.arg(type)

  x <- spendings(season)[standings(season)]

  o <- switch(type,
              "CR" = order(x, decreasing = TRUE),
              "OR" = seq(along = x))

  c(0, sapply(seq(along = x),
              function(g) concentration_ratio(x, g, o)))

}


seasons <- levels(BundesligaFinalStandings$Season)[-c(1:2)]

CR_curves <- sapply(seasons, cc, type = "CR")
OR_curves <- sapply(seasons, cc, type = "OR")


## Figure~2 and Figure~3:
matplot(CR_curves, type = "l", lty = 1)
matplot(OR_curves, type = "l", lty = 1)



### Concentration indices (Table~1):

ci <- function(season) {
  x <- spendings(season)[standings(season)]
  c(H = herfindahl(x),
    RB = rosenbluth(x),
    OI = exogeny(x))
}

CI <- t(sapply(seasons, ci))
CI


### Equivalent number of ... (Table~2):

1/CI
