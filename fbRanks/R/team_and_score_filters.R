team.and.score.filters = function(x,...){
#include.teams is a vector of teams names based on the filter info in ...
#include.scores is a TRUE,FALSE vector for each row in x$scores
  
#x is a list with appropriate data.frames for teams and scores, could be fbRanks object
  if(!all(c("scores","teams") %in% names(x))) stop("team.and.score.filters function requires a list with scores and teams data frames.\n", call.=FALSE)
  team.data = x$teams
  scores=x$scores
  #names are not case sensitive
  names(team.data)=tolower(names(team.data))
  names(scores)=tolower(names(scores))

  if(!any(names(team.data)=="name")) stop("The teams data frame requires a column named \"name\" for the team (display) name.\n",call.=FALSE)
  if(!all(c("home.team","away.team") %in% names(scores))) stop("The scores data frame requires columns named \"home.team\" and \"away.team\".\n",call.=FALSE)
  
  #Set up the team names; used for printing
  cc = names(team.data)=="name"
  all.team.names=team.data[,cc]
  nteams = length(all.team.names)
  #Find any extra arguments passed in and check if it corresponds to columns in the team file or match file
  extra=list(...)
  names(extra)=tolower(names(extra))
  if(any(duplicated(names(extra)))){
  duplicated.extra=names(extra)[duplicated(names(extra))]
  stop(paste("There are duplicated column names passed in for filtering. Names are case insensitive.\nThe duplicated names are:",paste(duplicated.extra,collapse=", ")),call.=FALSE)
  }
  if(any(!(names(extra) %in% c(names(team.data),names(scores)))) ){
    cat("Extra arguments (besides team and venue) passed in should correspond to columns in the team or match files.\n")
    cat("The following extra arguments are not in the team or match files (names are not case sensitive): ")
    bad.vals = names(extra)[!(names(extra) %in% c(names(team.data),names(scores)))]
    cat(paste(bad.vals,collapse=", "))
    cat("\n")
    stop(call.=FALSE)
  } 
  if(any((names(extra) %in% names(team.data)) & (names(extra) %in% names(scores))) ){
    cat("You cannot filter on names that appear as column names in BOTH the team file and match file.\n")
    cat("The following appear as column names in both team or match files (names are case sensitive): ")
    bad.vals = names(extra)[(names(extra) %in% names(team.data)) & (names(extra) %in% names(scores))]
    cat(paste(bad.vals,collapse=", "))
    cat("\n")
    stop(call.=FALSE)
  } 
  extra.team = extra[names(extra) %in% names(team.data)]
  extra.scores = extra[names(extra) %in% names(scores)]
  
  #Set up the extra team filters the original values are in extra.team
  extra.team.filter = extra.team
  for(j in names(extra.team)){
    cc = names(team.data)==j
    if(identical(tolower(extra.team.filter[[j]]),"all")){
      extra.team.filter[[j]] = unique(as.vector(as.matrix(team.data[,cc])))
    }
    if(all(!(extra.team.filter[[j]] %in% unique(as.vector(as.matrix(team.data[,cc])))))){
      cat(paste("None of the values passed in for ", j, " are in the team file.\n",sep=""))
      cat("Exiting since nothing would be printed.\n")
      stop(call.=FALSE)
    }
    if(any(!(extra.team.filter[[j]] %in% unique(as.vector(as.matrix(team.data[,cc])))))){
      bad.vals = extra.team.filter[[j]][!(extra.team.filter[[j]] %in% unique(as.vector(as.matrix(team.data[,cc]))))]
      cat(paste("FYI: values ", paste(bad.vals, collapse=",")," passed in for ", j, " are not found in the ", j," column of the ", as.character(match.call()$x),"$teams dataframe.", sep=""))
    }
  }
  #Set up the extra score filters the original values are in extra.scores
  extra.scores.filter = extra.scores
  for(j in names(extra.scores)){
    cc = names(scores)==j
    if(identical(tolower(extra.scores.filter[[j]]),"all")){
      extra.scores.filter[[j]] = unique(as.vector(as.matrix(scores[,cc])))
    }
    if(all(!(extra.scores.filter[[j]] %in% unique(as.vector(as.matrix(scores[,cc])))))){ 
      cat(paste("None of the values passed in for ", j, " are in the ", as.character(match.call()$x),"$scores dataframe.\n",sep=""))
      cat("Exiting since nothing would be printed.\n")
      stop(call.=FALSE)
    }
    if(any(!(extra.scores.filter[[j]] %in% unique(as.vector(as.matrix(scores[,cc])))))){
      bad.vals = extra.scores.filter[[j]][!(extra.scores.filter[[j]] %in% unique(as.vector(as.matrix(scores[,cc]))))]
      cat(paste("FYI: values ", paste(bad.vals, collapse=",")," passed in for ", j, " are not found in the ", j," column of the ", as.character(match.call()$x),"$scores dataframe.", sep=""))
      cat("\nProceeding using the other values.\n\n")
    }
  }
  
  # Set up the team name filter
  team.filter = "all" 
  if(identical(tolower(team.filter),"all")){
    cc = tolower(names(team.data))=="name"
    team.filter = unique(as.vector(as.matrix(team.data[,cc])))
  } 
  
  #which subset of teams to include in ranks printing
  tmp.fun=function(x,y){ any(x %in% y) }
  include.extra.team = rep(TRUE,nteams)
  for(j in names(extra.team)){
    filt=extra.team.filter[[j]]
    include.extra.team = include.extra.team & (team.data[[j]] %in% filt) 
  }
  
  #Filtering on columns in the match file is different because we need to find
  #each a set of teams associated with the values in the filter and add those teams to our list of teams to print
  include.extra.scores = rep(TRUE,nteams)
  label.extra.scores = list()
  for(j in names(extra.scores)){
    label.extra.scores[[j]]=list()
    label.extra.scores[[j]][all.team.names]=""
    teams.in.this.filter = list() #the list of teams associated with each value given for the name in extra.scores
    for(i in extra.scores.filter[[j]]){  #the values for name j in the scores filter
      rows.to.include = scores[scores[[j]]==i,c("home.team", "away.team")]
      teams.in.this.filter[[i]] = unique(as.vector(as.matrix(rows.to.include)))
      label.extra.scores[[j]][ teams.in.this.filter[[i]] ]=lapply(label.extra.scores[[j]][teams.in.this.filter[[i]]],paste,i,sep=",")
    }
    label.extra.scores[[j]]=lapply(label.extra.scores[[j]],str_sub,2)
    teams.in.this.filter=unique(unlist(teams.in.this.filter)) #all teams that match that filter
    include.extra.scores = include.extra.scores & (all.team.names %in% teams.in.this.filter)
  }
  
  #filtering on team name is easy; we just see if the team name appears in the team filter 
  cc = tolower(names(team.data))=="name"
  include.team = apply(team.data[,cc,drop=FALSE], 1, tmp.fun, team.filter) 
  
  #this is a vector of team names
  include.teams=all.team.names[(include.team & include.extra.scores & include.extra.team)]
  
  #this gets a score filter for any extra names that appear in the scores file
  scores.filter = names(extra)[names(extra) %in% names(scores)]
  include.scores=scores$home.team %in% include.teams | scores$away.team %in% include.teams #for each row
  for(sfil in scores.filter){
    include.scores=include.scores & (scores[[sfil]] %in% extra[[sfil]])
  }
  
  return(list(include.teams=include.teams, include.scores=include.scores))
}