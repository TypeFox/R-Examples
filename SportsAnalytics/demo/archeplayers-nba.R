#' @demo Archetypal basketball players based on player statistics;
#'   analysis from the manuscript "Archetypal athletes" by Eugster
#'   (2011)

library("SportsAnalytics")
library("archetypes")
library("RColorBrewer")

col_pal <- brewer.pal(7, "Set1")
col_black <- rgb(0, 0, 0, 0.2)



### Data: ############################################################

data("NBAPlayerStatistics0910")

dat <- subset(NBAPlayerStatistics0910,
              select = -c(Ejections, FlagrantFouls))

mat <- as.matrix(subset(dat, select = -c(League, Name, Team, Position)))


pcplot(mat, col = col_black, las = 2)



### Archetypes: ######################################################

set.seed(4321)
as <- stepArchetypes(mat, k = 1:10)

rss(as)
screeplot(as)


a4 <- bestModel(as[[4]])


### Archetypal basketball players:

parameters(a4)
barplot(a4, mat, percentiles = TRUE)



### Player interpretation: ###########################################

players <- function(which) {
  players <- list()
  players$which <- which
  players$mat <- mat[which, ]
  players$coef <- coef(a4, "alphas")[which, ]
  players$dat <- dat[which, ]

  players
}


### Archetypal players:

which <- apply(coef(a4, "alphas"), 2, which.max)
atypes <- players(which)

cbind(subset(atypes$dat, select = c(Name, Team, Position)),
      atypes$coef)



### Good players:

good_players <- function(atype, threshold) {
  which <- which(coef(a4, "alphas")[, atype] > threshold)

  good_coef <- coef(a4, "alphas")[which, ]
  good_dat <- subset(dat[which, ], select = c(Name, Team, Position))
  good_dat <- cbind(good_dat, good_coef)
  good_dat <- good_dat[order(-good_coef[, atype]), ]

  good_dat
}


good_threshold <- 0.95

players <- lapply(2:4, good_players, good_threshold)
players
