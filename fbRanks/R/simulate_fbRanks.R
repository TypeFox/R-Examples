################################################################
## simulate matches using a fbRanks model
# User can pass in newdata or specify dates and other filter info to construct
# newdata from the scores data.frame
# in the fbRanks object
# groups.column, if newdata has brackets/groups/leagues you can specify which column to use for that info
#    example, groups.column="venue"
# show.matches, show the predicted results of the matches
# points.rule, specify the rule for points.  This can be a text string if coded in or a list if not
# tie.rule, currently just GD with a GD max specified by tie.rule.gd.max
# newdata, same as for predict()
################################################################
simulate.fbRanks=function(object, nsim=100000, seed=NULL, ..., newdata=list(home.team="foo", away.team="bar"),
                         bracket.names=NULL,
                         max.date="2100-6-1", min.date="1900-5-1", silent=FALSE,
                         points.rule="tournament 10pt",tie.rule=list(tie.rule.gd.max=10),
                         non.equal.games.rule=list(n.games=3,rule="proportional"),
                         show.matches=FALSE, groups.column=NULL){
  x=object
  
  if(!missing(bracket.names)){
    if(!is.list(bracket.names)){
      stop("bracket.names must be a list with each list element a different vector of team names",call.=FALSE)
    }
    if(!missing(groups.column))
      cat("You passed in both bracket.names and groups.column.  Bracket.names will be ignored and groups.column used.\n")
  }
  
  if(!is.list(non.equal.games.rule)){
    stop("non.equal.games.rule must be a list with n.games (a number) and rule (a character name).\n",call.=FALSE)
  }
  if(!all(c("n.games","rule") %in% names(non.equal.games.rule))){
    stop("non.equal.games.rule must be a list with n.games (a number) and rule (a character name).\n",call.=FALSE)
  }
  if(!(non.equal.games.rule$rule %in% c("proportional","ignore")))
    stop("non.equal.games.rule$rule must be proportional or ignore; that's the only rules coded currently.\n",call.=FALSE)

if(is.list(points.rule)){
  if(!all(c("pwin","pdraw","pshutout","pgoal","pgoal.max") %in% names(points.rule)))
    stop("If points.rule is a list it must have elements pwin, pdraw, pshutout, pgoal, pgoal.max.\n",call.=FALSE)
  pwin=points.rule$pwin; pdraw=points.rule$pdraw; pshutout=points.rule$pshutout
  pgoal=points.rule$pgoal; pgoal.max=points.rule$pgoal.max
}else{
  if(length(points.rule)!=1)
    stop("points.rule must be the name of a points rule or a list with elements pwin, pdraw, pshutout, pgoal, pgoal.max.\n",call.=FALSE)
  ok.points.rules=c("tournament.10pt","tournament10pt","tournament 10pt","league.3pt","league3pt","league 3pt")
  if(!(points.rule %in% ok.points.rules))
    stop("Known points rules are tournament.10pt and league.3pt.\nYou can specify different rules with a list with elements pwin, pdraw, pshutout, pgoal, pgoal.max.\n",call.=FALSE)
  # #Tournament points rule #1
  # a. 6 points for a win
  # b. 3 points for a draw
  # c. 0 (zero) points for a loss
  # d. 1 point for each goal scored (up to a maximum of 3 per game for both teams)
  # e. 1 point for a shutout - holding an opponent scoreless(in the event of a 0-0 tie, both teams will be awarded 4 points)
  if(points.rule %in% c("tournament.10pt","tournament10pt","tournament 10pt")){
    pwin=6; pdraw=3; pshutout=1; pgoal=1; pgoal.max=3
  }
  # #League points rule #1
  # a. 3 points for a win
  # b. 1 points for a draw
  # c. 0 (zero) points for a loss
  if(points.rule %in% c("league.3pt","league3pt","league 3pt")){
    pwin=3; pdraw=1; pshutout=0; pgoal=0; pgoal.max=0
  }
  }
  if(!is.list(tie.rule))
    stop("tie.rule must be a list with element tie.rule.gd.max.\n",call.=FALSE)
  if(is.null(tie.rule$tie.rule.gd.max))
    stop("tie.rule must be a list with element tie.rule.gd.max.\n",call.=FALSE) 
  #I don't have much of a tie rule entered; just GD rule
  tie.rule.gd.max=tie.rule$tie.rule.gd.max
  if(!missing(newdata)){
    scores=create.newdata.dataframe(x, newdata, min.date, max.date, ...)
  }else{ #newdata not given so we get the newdata from the fit object; teams names will be ok    
    scores = x$scores
#     #filter the scores data.frame based on any filters the user specified
#     include.teams=team.name.filter(list(scores=scores, teams=x$teams),...)
#     scores= scores[scores$home.team %in% include.teams | scores$away.team %in% include.teams,]
#     extras=list(...)
#     scores.filter = names(extras)[names(extras) %in% names(scores)]
#     for(sfil in scores.filter){
#       scores=scores[scores[[sfil]] %in% extras[[sfil]],]
#     }

    #determine which teams to include based on any filter arguments the user included
    include.scores=team.and.score.filters(list(scores=scores, teams=x$teams),...)$include.scores
    scores=scores[include.scores,,drop=FALSE]
    
    #filter on the date range
    #scores date is already a date class
    scores = scores[scores$date>=min.date,,drop=FALSE]
    scores = scores[scores$date<=max.date,,drop=FALSE]
  }
  
  if(!missing(groups.column)){
    if(length(groups.column)!=1)
      stop("groups.column must be a single column name in newdata or scores.\n",call.=FALSE)
    if(!(groups.column %in% names(scores)))
      stop("groups.column must be a column name in newdata or scores.\n",call.=FALSE)
  }
  
#fbRanks object to use;  this is what is output by a call to rank.teams()
fbRanks=x

#set up the groupings names
if(is.null(groups.column)){ 
  glevels=0
  scores.teams=unique(c(as.character(scores$home.team),as.character(scores$away.team)))
  if(is.null(bracket.names)){
    bracket.names=list(Bracket=scores.teams)
}else{
  if(!all(unlist(bracket.names) %in% scores.teams))
    stop("Not all the bracket.names are in the newdata scores file.",call.=FALSE)  
}
  }else{
if(is.factor(scores[groups.column])){ glevels=levels(scores[groups.column])
}else{ glevels=unique(scores[[groups.column]]) }
bracket.names=list()
for(gval in glevels){
  bscores=scores
  bscores=scores[scores[groups.column]==gval,,drop=FALSE]
  bteams=unique(c(as.character(bscores$home.team),as.character(bscores$away.team)))
  bracket.names[[gval]]=bteams
}
}
if(is.null(names(bracket.names)))
  names(bracket.names)=paste("Bracket",LETTERS[1:length(bracket.names)])

all.standings=list()
glevels=0  #4-7 hack to use bracket.names list instead
for(gval in glevels){
  bscores=scores
  if(!identical(glevels,0)) bscores=scores[scores[groups.column]==gval,,drop=FALSE]
  #Predict using the Founders bracket info
  pout=predict(fbRanks,newdata=bscores,silent=TRUE,n=nsim)
  #team names for this bracket
  bteams=unique(c(as.character(pout$scores$home.team),as.character(pout$scores$away.team)))

  n.teams=length(bteams)  #in bracket
  n.games=dim(pout$scores)[1] #in bracket
  #Set up a holder for the home team for each of the bracket games
  points.home=matrix(0,n.games,nsim) 
  rownames(points.home)=rownames(pout$home.goals) #rownames will be the home team names for each game in bracket
  points.home[pout$home.goals>pout$away.goals]=pwin #points for a win
  points.home[pout$home.goals==pout$away.goals]=pdraw #points for a draw
  points.home[pout$away.goals==0]=points.home[pout$away.goals==0]+pshutout  #points for a shutout
  points.home=points.home+pgoal*pmin(matrix(pgoal.max,n.games,nsim),pout$home.goals)   #points for each goal
  gd.home=pout$home.goals-pout$away.goals
  gd.home[gd.home>tie.rule.gd.max]=tie.rule.gd.max
  
  #repeat for the away team
  points.away=matrix(0,n.games,nsim)
  rownames(points.away)=rownames(pout$away.goals)
  points.away[pout$away.goals>pout$home.goals]=pwin #points for a win
  points.away[pout$away.goals==pout$home.goals]=pdraw #points for a draw
  points.away[pout$home.goals==0]=points.away[pout$home.goals==0]+pshutout #points for a shutout
  points.away=points.away+pgoal*pmin(matrix(pgoal.max,n.games,nsim),pout$away.goals) #points for each goal
  gd.away=pout$away.goals-pout$home.goals
  gd.away[gd.away>tie.rule.gd.max]=tie.rule.gd.max
  
  #set up a holder for the points for each team in bracket/venue from each simulation
  points=gds=matrix(0,n.teams,nsim)
  rownames(points)=rownames(gds)=bteams  #rownames are now the team names
  #for each team add together the points it got for games when it was home and games when it was away
  for(team in bteams){
    #number of games per team because teams don't play the same number
    n.games.team = sum(rownames(points.home) %in% team) + sum(rownames(points.away) %in% team)
    #sum up the points when playing as 'home'
    if(team %in% rownames(points.home)){
      points[team,]=apply(points.home[rownames(points.home)==team,,drop=FALSE],2,sum)
      gds[team,]=apply(gd.home[rownames(gd.home)==team,,drop=FALSE],2,sum)
    }
    #sum up the points when playing as 'away'
    if(team %in% rownames(points.away)){
      points[team,]=points[team,]+apply(points.away[rownames(points.away)==team,,drop=FALSE],2,sum)
      gds[team,]=gds[team,]+apply(gd.away[rownames(gd.away)==team,,drop=FALSE],2,sum)
    }
    #Some teams play different numbers of games.  Need to adjust the points in that case
    if(non.equal.games.rule$rule=="proportional"){
      points[team,]=points[team,]*(non.equal.games.rule$n.games/n.games.team)
    }
  }
} #hack 4-7-2013; set glevels to 0 so everything is simed together;
  #afterwards apply the brackets
  for(gval in names(bracket.names)){
    bteams=bracket.names[[gval]]
  n.teams=length(bteams)
  #This is a hack to eliminate any ties (after GD tie breaker rule). Basically it applies a coin-flip to ties
    #this is a hack since this isn't how tournaments do this.  They use a head-to-head rule then GF rule usually
  points=points+apply(gds,2,function(x){x/(10*(sum(x-min(x))))+runif(length(x),-1,1)/1000})
  
  #get just the points for the group/bracket
  bpoints=points[bteams,]
  #Now get the standings; this is terse code; 
  #what it is doing is for each column, run the sort function and get the ordering
  standing.function=function(x){
    #this applies a gd tie breaker with random selection if the there are ties after applying the GD breaker
    places=1:n.teams
    places[sort(x,index.return=TRUE,decreasing=TRUE)$ix]=places
    places
  }
  standings=apply(bpoints,2,standing.function)
  rownames(standings)=bteams
  
  #Set up the matrix that will hold the standings
  pstandings=matrix(0,n.teams,n.teams)
  rownames(pstandings)=bteams
  #column names will the 1st-nth
  colnames(pstandings)=paste(1:n.teams,c("st","nd","rd",rep("th",length(bteams)-3)),sep="")

  for(team in bteams){
    #need to wrap standings in factor so I can set levels.  Otherwise won't get the 0s when say a team is never 4th
    pstandings[team,]=table(factor(standings[team,],levels=1:length(bteams)))/nsim
    }
  
  if(!silent){
    cat("\n");
    if(length(bracket.names)!=1) cat(gval)
    cat("\n");
  print(round(pstandings*100,digits=0))
    cat("\n");
  }
  #Predict games
  if(show.matches){ cat("\n"); predict(fbRanks,newdata=bscores) }
  
  if(length(bracket.names)!=1) all.standings[[gval]]=pstandings else all.standings=pstandings
}

  invisible(all.standings)
}