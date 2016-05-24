#This will take a fbRanks object along with a newdata object
#and create a scores dataframe with the display names in x$teams
create.newdata.dataframe=function(x, newdata, min.date, max.date, ...){  
    if(is.null(newdata$home.team)) stop("You need to specify the home teams.\n")
    if(is.null(newdata$away.team)) stop("You need to specify the away teams.\n")
    if(is.null(newdata$home.score)) newdata$home.score=NaN
    if(is.null(newdata$away.score)) newdata$away.score=NaN 
    if(is.null(newdata$date)) newdata$date=NA
    for(i in names(newdata))
      if(is.factor(newdata[[i]])) newdata[[i]]=as.character(newdata[[i]])
    scores=data.frame(newdata, stringsAsFactors=FALSE)
    #because the user passed in team names, we need to match names against the team file
    #and report any problems
    tmp.resolver=data.frame(name=x$teams$name, alt.name=x$teams$name, stringsAsFactors=FALSE)
    tmp=resolve.team.names(scores, tmp.resolver, team.data=NULL, use.team.select=FALSE)
    if(!tmp$ok) stop("Bad team names in newdata argument.  They do not match names in the teams data frame.\n", call.=FALSE)
    scores=tmp$scores
           
    #determine which teams to include based on any filter arguments the user included
    include.scores=team.and.score.filters(list(scores=scores, teams=x$teams),...)$include.scores
    scores=scores[include.scores,,drop=FALSE]
    
    #Fix the date
    if(!all(is.na(scores$date))){
      scores$date=as.Date(scores$date, x$date.format) 
      if(any(is.na(scores$date))) stop(paste("scores dates must be entered in the following format:",format(Sys.Date(),x$date.format),"\n"))
      scores = scores[scores$date>=min.date,,drop=FALSE]
      scores = scores[scores$date<=max.date,,drop=FALSE]
    }
    return(scores)
  }
    