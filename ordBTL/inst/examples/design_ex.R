# load german football-league (Bundesliga) data
library("wikibooks")
data(Bundesliga)

# add new variable Y3 reflecting the response which is coded as 
# 1 if the home team wins
# 2 if the game ends up with a tie
# 3 if the home team loses
diff <- Bundesliga$Tore.Heim - Bundesliga$Tore.Gast
Bundesliga$Y3 <- as.ordered(ifelse(diff >= 1, 1, 
                                   ifelse(diff <= -1, 3, 2)))
buli0506 <- subset(Bundesliga, Saison=="2005/2006")
str(buli0506)

# Design matrix without home advantage
des.nohome <- design(buli0506, var1="Heim", var2="Gast", 
                     home.advantage="no")
str(des.nohome)

# Design matrix with one home advantage parameter for all objects
des.onehome <- design(buli0506, var1="Heim", var2="Gast", 
                      home.advantage="yes")
str(des.onehome)

# Design matrix with home advantage parameters for each object
des.teamhome <- design(buli0506, var1="Heim", var2="Gast",
                       home.advantage="specific")
str(des.teamhome)

# Design matrix with additional covariable "Spieltag"
des.covs <- design(buli0506, var1="Heim", var2="Gast", 
                   use.vars=c("Spieltag"), home.advantage="no")
str(des.covs)
