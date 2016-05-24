#' @demo Two-dimensional archetypal basketball players; illustrative
#'   example from manuscript "Archetypal athletes" by Eugster (2011)

library("SportsAnalytics")
library("archetypes")
library("RColorBrewer")
library("vcd")

col_pal <- brewer.pal(7, "Set1")
col_black <- rgb(0, 0, 0, 0.2)



### Data: ############################################################

data("NBAPlayerStatistics0910")

dat <- subset(NBAPlayerStatistics0910,
              select = c(Team, Name, Position,
                         TotalMinutesPlayed, FieldGoalsMade))

mat <- as.matrix(subset(dat, select = c(TotalMinutesPlayed, FieldGoalsMade)))
rownames(mat) <- NULL


plot(mat, pch = 19, col = col_black)



### Archetypes: ######################################################

set.seed(4321)
as <- stepArchetypes(mat, k = 1:5)

screeplot(as)

a3 <- bestModel(as[[3]])


### Archetypal basketball players:

parameters(a3)
barplot(a3, mat, percentiles = TRUE)


### Visualization:

xyplot(a3, mat, data.col = col_black,
       chull = chull(mat), atypes.col = col_pal[1])
points(parameters(a3), pch = 4, col = col_pal[1])


### Alpha coefficients:

ternaryplot(coef(a3, "alphas"), dimnames = 1:3, pch = 19,
            border = "gray", col = col_black, main = NULL,
            labels = "outside", grid = TRUE)



### Player interpretation: ###########################################

players <- function(which) {
  players <- list()
  players$which <- which
  players$mat <- mat[which, ]
  players$coef <- coef(a3, "alphas")[which, ]
  players$dat <- dat[which, ]

  players
}

players_xyplot <- function(players, col) {
  xyplot(a3, mat, data.col = col_black,
         atypes.col = col_pal[1], atypes.pch = NA)

  points(parameters(a3), pch = 4, col = col_pal[1])
  text(parameters(a3), labels = sprintf("A%s", 1:3), pos = c(1, 3, 4))

  points(players$mat, pch = 19, col = col)
}

players_ternaryplot <- function(players, col) {
  cols <- rep(col_black, nrow(mat))
  cols[players$which] <- col

  ternaryplot(coef(a3, "alphas"), pch = 19, dimnames = 1:3,
              border = "gray", col = cols, main = NULL,
              labels = "outside", grid = TRUE)
}



### Archetypal players:

which <- apply(coef(a3, "alphas"), 2, which.max)
atypes <- players(which)

atypes$dat
atypes$coef

players_xyplot(atypes, col_pal[1])
players_ternaryplot(atypes, col_pal[1])



### Players near archetype 1:

which <- which(coef(a3, "alphas")[, 1] > 0.8)
atypes1 <- players(which)

atypes1$dat
atypes1$coef

players_xyplot(atypes1, col_pal[2])
players_ternaryplot(atypes1, col_pal[2])



### Players near archetype 2:

which <- which(coef(a3, "alphas")[, 2] > 0.8)

set.seed(1234)
which <- sample(which, 5)

atypes2 <- players(which)



atypes2$dat
atypes2$coef

players_xyplot(atypes2, col_pal[5])
players_ternaryplot(atypes2, col_pal[5])



### Players near archetype 3:

which <- which(coef(a3, "alphas")[, 3] > 0.8)
atypes3 <- players(which)

atypes3$dat
atypes3$coef

players_xyplot(atypes3, col_pal[3])
players_ternaryplot(atypes3, col_pal[3])



### Random selection of players in the data sets' center:

which <- which(apply(coef(a3, "alphas") < 0.5, 1, all))

set.seed(1234)
which <- sample(which, 5)

atypes0 <- players(which)

atypes0$dat
atypes0$coef

players_xyplot(atypes0, col_pal[7])
players_ternaryplot(atypes0, col_pal[7])

