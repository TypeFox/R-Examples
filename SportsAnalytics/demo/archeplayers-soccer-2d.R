#' @demo Two-dimensional archetypal soccer players; illustrative
#'   example from talk "On the power of modern statistical methodology
#'   in soccer analysis -- Archetypal soccer players" by Eugster,
#'   Abedieh, Schnell, and Augustin (2011)

library("SportsAnalytics")
library("archetypes")
library("RColorBrewer")
library("vcd")

col_pal <- brewer.pal(3, "Set1")
col_black <- rgb(0, 0, 0, 0.2)



### Data: ############################################################

data("EURO4PlayerSkillsSep11")


### Two-dimensional dataset:

dat <- subset(EURO4PlayerSkillsSep11,
              Position != "Goalkeeper" & TopSpeed > 0,
              select = c(Name, Team, Position, Technique, TopSpeed))

mat <- as.matrix(subset(dat, select = c(Technique, TopSpeed)))
rownames(mat) <- NULL


plot(mat, pch = 19, col = col_black,
     xlim = c(60, 100), ylim = c(60, 100))



### Archetypes: ######################################################

set.seed(1234)
as <- stepArchetypes(mat, k = 1:5)

screeplot(as)

a3 <- bestModel(as[[3]])


### Archetypal soccer players:

parameters(a3)
barplot(a3, mat, percentiles = TRUE)



### Visualization:

xyplot(a3, mat, data.col = col_black,
       # chull = chull(mat),
       atypes.col = col_pal[1], atypes.pch = NA,
       xlim = c(60, 100), ylim = c(60, 100))

points(parameters(a3), pch = 4, col = col_pal[1])
text(parameters(a3), labels = sprintf("A%s", 1:3), pos = c(1, 3, 4))



### Alpha coefficients:

ternaryplot(coef(a3, "alphas")[, c(1, 3, 2)],
            dimnames = sprintf("A%s", c(1, 3, 2)),
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
         atypes.col = col_pal[1], atypes.pch = NA,
         xlim = c(60, 100), ylim = c(60, 100))

  points(parameters(a3), pch = 4, col = col_pal[1])
  text(parameters(a3), labels = sprintf("A%s", 1:3), pos = c(1, 3, 4))

  points(players$mat, pch = 19, col = col)
}

players_ternaryplot <- function(players, col) {
  cols <- rep(col_black, nrow(mat))
  cols[players$which] <- col

  ternaryplot(coef(a3, "alphas")[, c(1, 3, 2)],
              dimnames = sprintf("A%s", c(1, 3, 2)), pch = 19,
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



### Best players:

which <- which(coef(a3, "alphas")[, 1] == 0)
best <- players(which)

best$dat
best$coef

players_xyplot(best, col_pal[3])
players_ternaryplot(best, col_pal[3])



### Random players:

#which <- sample(5, nrow(mat))
which <- c(170, 217, 429, 1512, 1610)
rand <- players(which)

rand$dat
rand$coef

players_xyplot(rand, col_pal[2])
players_ternaryplot(rand, col_pal[2])


