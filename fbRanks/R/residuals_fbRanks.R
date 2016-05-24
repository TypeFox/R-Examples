residuals.fbRanks = function(object, ...){
x=object
team.data = x$teams
glm.fits=x$fit
scores=x$scores

#Set up the team names; used for printing
cc = tolower(names(team.data))=="name"
all.team.names=team.data[,cc]
nteams = length(all.team.names)

#determine which teams to include based on any filter arguments the user included
ts.filters=team.and.score.filters(list(scores=scores, teams=team.data),...)
include.teams=ts.filters$include.teams

team.residuals=list()
include.teams=sort(include.teams) #make alphabetical
#Return a list with an element for each team
for(i in include.teams){
  #need to match up the requested residuals against the residuals column
  # what matches should we show the residuals for
  # include.teams is what teams the user requested based on filters (league, name, etc.)
  match.in.request = (scores$home.team %in% i) | (scores$away.team %in% i)
  team.residuals[[i]]=scores[match.in.request,!(names(scores) %in% c("home.residuals","away.residuals"))]
  theresids = scores$home.residuals[match.in.request]
  theresids[team.residuals[[i]]$away.team %in% i]=scores$away.residuals[match.in.request & (scores$away.team %in% i)]
  team.residuals[[i]]$attack.residuals = theresids
  theresids = scores$away.residuals[match.in.request]
  theresids[team.residuals[[i]]$away.team %in% i]=scores$home.residuals[match.in.request & (scores$away.team %in% i)]
  team.residuals[[i]]$defense.residuals = theresids
  
}


return(team.residuals)
}