#This is a utility function to replace the names in a scores data frame with the display names in teams data frame
resolve.team.names=function(scores, team.resolver, team.data=NULL, use.team.select=TRUE){

  if(use.team.select) have.tcltk=require(tcltk, quietly=TRUE) else have.tcltk=FALSE
  if(missing(team.data) | is.null(team.data)){ 
    missing.team.data = TRUE
    team.data=data.frame(name=team.resolver$name, stringsAsFactors=FALSE)
  }else{ missing.team.data=FALSE }
  #Set up the team names
  names(team.resolver)=tolower(names(team.resolver))
  if(!is.null(team.data)) names(team.data)=tolower(names(team.data))
  team.resolver$name=str_trim(team.resolver$name)
  display.names=unique(team.resolver$name)
  nteams = length(unique(display.names))
  
  #Check if there are problem with the team resolver
  updated=FALSE
  bad.names = team.resolver$name[!(team.resolver$name %in% team.data$name)]
  if(length(bad.names)!=0){
    if(use.team.select & missing.team.data){
      cat("In order to use the team select GUI, you need to pass in the team info dataframe.\n")
    }
    if(use.team.select & !have.tcltk){
      cat("\nTo use the GUI to select missing team names, you need to install tcltk package.\nOr set use.team.select=FALSE to turn off this warning.\n")
    }
    if(use.team.select & have.tcltk & !missing.team.data){
      tmp=team.name.select(unique(bad.names), team.resolver, team.data, scores, type="disp.name")
      team.resolver=tmp$team.resolver
      team.data=tmp$team.data
      updated=tmp$updated
    }
  } 
  
  bad.names = team.resolver$name[!(team.resolver$name %in% team.data$name)]
  if(length(bad.names)!=0){
    cat("The following team names in team resolver don't match names in team file:\n")
    cat(as.character(bad.names),sep=", ")
    cat("\nTeam resolver and team data dataframes being returned.\n")
    return(list(scores=scores, ok=FALSE, team.resolver=team.resolver, team.data=team.data, updated=updated))
  }

  #Set up alternate names list
alt.names = list()
for(i in display.names){
  #row where name is display.names
  alt.names[[ i ]] = as.character(team.resolver$alt.name[team.resolver$name==i])
  empty.alt.name = (str_trim(alt.names[[ i ]]) == "") | is.na(str_trim(alt.names[[ i ]]))
  alt.names[[ i ]] = alt.names[[ i ]][!empty.alt.name]
  alt.names[[ i ]] = str_strip.white(alt.names[[ i ]])
}
alt.names = lapply(alt.names,unique)
if(any(duplicated(unlist(alt.names),incomparables=character(0)))){
  cat("Something is wrong with the team resolver file.  There are duplicated alt names for different teams.\n")
  cat("Duplicated alt names are ")
  cat(paste(unlist(alt.names)[which(duplicated(unlist(alt.names),incomparables=character(0)))], collapse=", ")); cat("\n")
  stop()
}
  
#Replace any alternate names with the display name
tmp.fun=function(x,y){ unlist(lapply(lapply(y,is.element,x),any)) }
for(col in c("home.team","away.team")){
  bad.name=c(); bad.line=c()
  not.display.name=which(!(scores[[col]] %in% display.names))
  clean.alt.names = lapply(alt.names,tolower)
  for(i in not.display.name){
    alt.matches = tmp.fun(tolower(str_trim(scores[[col]][i])), clean.alt.names)
    if(sum(alt.matches) == 1) scores[[col]][i] = names(alt.names[alt.matches])
  }
}

  
#Check for bad names
bad.names=c(scores$home.team[!(scores$home.team %in% display.names)], scores$away.team[!(scores$away.team %in% display.names)])
  
  if(length(bad.names)!=0){
  if(use.team.select & is.null(team.data)){
    cat("In order to use the team select GUI, you need to pass in the team info dataframe.\n")
  }
  if(use.team.select & !have.tcltk){
    cat("\nTo use the GUI to select missing team names, you need to install tcltk package.\nOr set use.team.select=FALSE to turn off this warning.\n")
  }
  if(use.team.select & have.tcltk & !is.null(team.data)){
    tmp=team.name.select(unique(bad.names), team.resolver, team.data, scores)
    team.resolver=tmp$team.resolver
    team.data=tmp$team.data
    bad.names=tmp$skipped.teams
    updated=updated | tmp$updated  #might have set updated earlier
  }
}

  if(updated){ 
    #then we need to redo the replacing alt names because some of the away names are now in the team resolver
  # Now do the away since the team resolver has been updated teams.  Will repeat this for the away teams
  #Set up the team names
  names(team.resolver)=tolower(names(team.resolver))
  team.resolver$name=str_trim(team.resolver$name)
  display.names=unique(team.resolver$name)
  nteams = length(unique(display.names))
  
  #Set up alternate names list
  alt.names = list()
  for(i in display.names){
    #row where name is display.names
    alt.names[[ i ]] = as.character(team.resolver$alt.name[team.resolver$name==i])
    empty.alt.name = (str_trim(alt.names[[ i ]]) == "") | is.na(str_trim(alt.names[[ i ]]))
    alt.names[[ i ]] = alt.names[[ i ]][!empty.alt.name]
    alt.names[[ i ]] = str_strip.white(alt.names[[ i ]])
  }
  alt.names = lapply(alt.names,unique)
  if(any(duplicated(unlist(alt.names),incomparables=character(0)))){
    cat("Something is wrong with the team resolver file.  There are duplicated alt names for different teams.\n")
    cat("Duplicated alt names are ")
    cat(paste(unlist(alt.names)[which(duplicated(unlist(alt.names),incomparables=character(0)))], collapse=", ")); cat("\n")
    stop()
  }
  
  #Replace any alternate names with the display name
  tmp.fun=function(x,y){ unlist(lapply(lapply(y,is.element,x),any)) }
  for(col in c("home.team","away.team")){
    bad.name=c(); bad.line=c()
    not.display.name=which(!(scores[[col]] %in% display.names))
    clean.alt.names = lapply(alt.names,tolower)
    for(i in not.display.name){
      alt.matches = tmp.fun(tolower(str_trim(scores[[col]][i])), clean.alt.names)
      if(sum(alt.matches) == 1) scores[[col]][i] = names(alt.names[alt.matches])
    }
  }
  } #end of updating the alt.names again

  #Now check to see if there are any bad names left
  ok=TRUE
  bad.names=scores$home.team[!(scores$home.team %in% display.names)]
  if(length(bad.names)!=0){
      cat("The following home team names don't match names in team file:\n")
      cat(as.character(bad.names),sep=", ")
      cat("\nerrors are on lines ")
      cat((1:dim(scores)[1])[scores$home.team %in% bad.names]); cat("\n")
      ok=FALSE
  }
  bad.names=scores$away.team[!(scores$away.team %in% display.names)]
  if(length(bad.names)!=0){
    cat("The following away team names don't match names in team file:\n")
    cat(as.character(bad.names),sep=", ")
    cat("\nerrors are on lines ")
    cat((1:dim(scores)[1])[scores$away.team %in% bad.names]); cat("\n")
    ok=FALSE
  }  
  if(!ok){ 
  cat("Bad team names in scores file.  Score file being returned.\n")
}
return(list(scores=scores, ok=ok, team.resolver=team.resolver, team.data=team.data, updated=updated))
}
