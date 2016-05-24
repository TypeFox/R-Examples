library(BradleyTerry2)

##  This reproduces the analysis in Sec 10.6 of Agresti (2002).

data(baseball, package = "BradleyTerry2")

##  Simple Bradley-Terry model, ignoring home advantage:
baseballModel1 <- BTm(cbind(home.wins, away.wins), home.team, away.team,
                      data = baseball, id = "team")
summary(baseballModel1)

##  Now incorporate the "home advantage" effect
baseball$home.team <- data.frame(team = baseball$home.team, at.home = 1)
baseball$away.team <- data.frame(team = baseball$away.team, at.home = 0)
baseballModel2 <- update(baseballModel1, formula = ~ team + at.home)
summary(baseballModel2)

##  Compare the fit of these two models:
anova(baseballModel1, baseballModel2)
