bta <- BTabilities(nfl.model)
nflRatings<- data.frame(
    team = rownames(bta),
    rating = bta[,"ability"],
    se = bta[,"s.e."],
    wins = as.vector(table(nfl$winner)),
    losses = as.vector(table(nfl$loser))
    )
rownames(nflRatings) = NULL

nfl$winnerRating <- nflRatings$rating[as.numeric(nfl$winner)]
nfl$loserRating <- nflRatings$rating[as.numeric(nfl$loser)]
nfl$upset <- nfl$loserRating > nfl$winnerRating
nflRatings[rev(order(nflRatings$rating)),]
