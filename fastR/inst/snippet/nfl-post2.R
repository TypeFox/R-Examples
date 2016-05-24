require(faraway)           # for ilogit(), the inverse logit
nfl$pwinner <- ilogit(nfl$winnerRating - nfl$loserRating)
# how big an upset was the Super Bowl?
nfl[nrow(nfl),]
