ncaa <- ncaa2010   # make a local copy (with shorter name)
# at a neutral site?
ncaa$neutralSite <- grepl('notes',ncaa$n,ignore.case=TRUE)
# did home team win?
ncaa$homeTeamWon <- ncaa$hscore > ncaa$ascore 
# remove teams that didn't play >= 5 at home and >=5 away
# (typically div II teams that played a few div I teams)
h <- table(ncaa$home); a <- table(ncaa$away)
deleteTeams <- c(names(h[h<=5]), names(a[a<=5]))
ncaa <- ncaa[!( ncaa$home %in% deleteTeams | 
    ncaa$away %in% deleteTeams ), ]
# remove unused levels from home and away factors
teams <- union(ncaa$home, ncaa$away)
ncaa$home <- factor(ncaa$home, levels=teams)
ncaa$away <- factor(ncaa$away, levels=teams)
