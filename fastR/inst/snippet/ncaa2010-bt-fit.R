# fit a Bradley-Terry model
require(BradleyTerry2)
ncaa.model <- BTm( cbind(homeTeamWon, 1-homeTeamWon), 
                   home, away, data=ncaa, refcat="Duke")
