# fit Bradley-Terry model
require(BradleyTerry2)
BTm(cbind(homeTeamWon,!homeTeamWon), Home, Visitor, 
          data=nfl, id='team') -> nfl.model
