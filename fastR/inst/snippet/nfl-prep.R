nfl <- nfl2007                         # shorten name of data set
head(nfl,3)
nfl$dscore <- nfl$HomeScore - nfl$VisitorScore
w <- which(nfl$dscore > 0)
nfl$winner <- nfl$Visitor; nfl$winner[w] <- nfl$Home[w]
nfl$loser <- nfl$Home; nfl$loser[w] <- nfl$Visitor[w]

# did the home team win?
nfl$homeTeamWon <- nfl$dscore > 0
head(nfl,3)
