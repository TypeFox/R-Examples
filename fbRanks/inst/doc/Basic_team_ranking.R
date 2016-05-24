### R code from vignette source 'Basic_team_ranking.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
options(prompt=" ", continue=" ", width=100)
library(fbRanks)


###################################################
### code chunk number 2: loadfile (eval = FALSE)
###################################################
## pkgpath=find.package("fbRanks")
## file.loc=paste(pkgpath,"\\doc\\scores-web.csv",sep="")
## file.copy(file.loc,".")


###################################################
### code chunk number 3: load-data
###################################################
temp=create.fbRanks.dataframes(scores.file="scores-web.csv")
scores=temp$scores


###################################################
### code chunk number 4: show-data
###################################################
head(scores[,1:5])


###################################################
### code chunk number 5: rank-teams
###################################################
#get the ranks without any explantory variables
ranks1=rank.teams(scores=scores)


###################################################
### code chunk number 6: print-ranks (eval = FALSE)
###################################################
## #print all ranks
## print(ranks1)


###################################################
### code chunk number 7: load-teams
###################################################
temp=create.fbRanks.dataframes(scores.file="scores-web.csv", teams.file="teams-web.csv")
scores=temp$scores
teams=temp$teams
head(teams[,c("name","age","region","fall.league")])


###################################################
### code chunk number 8: rank-teams2
###################################################
#get the ranks without any explantory variables
ranks2=rank.teams(scores=scores, teams=teams)


###################################################
### code chunk number 9: print-ranks2
###################################################
print(ranks2, fall.league="RCL D1 U12")


###################################################
### code chunk number 10: expvars
###################################################
names(scores)


###################################################
### code chunk number 11: rank-teams3 (eval = FALSE)
###################################################
## #get the ranks with surface effect; note surface must be a column 
## #in scores (or teams) dataframe for this to work
## ranks3=rank.teams(scores=scores,teams=teams,add="surface")


###################################################
### code chunk number 12: rank-teams4
###################################################
ranks4=rank.teams(scores=scores,teams=teams,add=c("surface","adv"))


###################################################
### code chunk number 13: coef.turf
###################################################
coef(ranks4$fit$cluster.1)["surface.fTurf"]


###################################################
### code chunk number 14: coef.home
###################################################
coef(ranks4$fit$cluster.1)["adv.fhome"]


###################################################
### code chunk number 15: rank-summer
###################################################
#Drop the home/away advantage since all the summer games are neutral
#Add max.date so tell the function to only use matches up to max.date
ranks.summer=rank.teams(scores=scores,teams=teams,add=c("surface"), max.date="2012-9-5")


###################################################
### code chunk number 16: simulate
###################################################
simulate(ranks.summer, venue="RCL D1")


###################################################
### code chunk number 17: predict
###################################################
predict(ranks.summer, venue="RCL D1", date=as.Date("2012-09-16"))


###################################################
### code chunk number 18: fantasy-tournament
###################################################
fantasy.teams=c("Seattle United Copa B00","Seattle United Tango B00",
                "Seattle United Samba B00","Seattle United S Black B00")
home.team=combn(fantasy.teams,2)[1,]
away.team=combn(fantasy.teams,2)[2,]
fantasy.games=data.frame(
  date="2013-1-1", 
  home.team=home.team, 
  home.score=NaN, 
  away.team=away.team,
  away.score=NaN, surface="Grass", 
  home.adv="neutral", away.adv="neutral")


###################################################
### code chunk number 19: simulate-fantasy
###################################################
simulate(ranks4, newdata=fantasy.games, points.rule="tournament10pt")


###################################################
### code chunk number 20: Reset
###################################################
options(prompt="> ", continue="+ ")


