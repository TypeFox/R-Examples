create.fbRanks.dataframes=function(scores.file, team.resolver=NULL, teams.file=NULL, date.format="%Y-%m-%d", na.remove=FALSE){
  #This creates scores and teams dataframes from the scores files and team files
  #team.resolver is 2 columns: name=team name (display name), alt.name=name in match file
  #team info is many columns.  name=team name (display name)
  updated=FALSE #if the team.resolver or teams.file changed with user input
  
  if(missing(scores.file)) stop("The name of a scores file (csv) must be passed in.\n",call.=FALSE)
  
  #read in the scores files (match files)
  #all the column names will be made lower and striped of extra white space
  scores=data.frame()
  scores.warning=list()
  for(ffile in scores.file){
    filename=ffile
    if(length(scores)==0){ 
      scores=read.csv(file=filename, colClasses=c("character"),strip.white=TRUE, stringsAsFactors=FALSE)
      colnames(scores)=tolower(str_strip.white(names(scores)))
      if(!all(c("date","home.team","home.score","away.team","away.score") %in% colnames(scores)))
         stop(paste(filename,"is missing date, home.team, home.score, away.team, or away.score column.\n"), call.=FALSE)
      if(na.remove) scores = scores[!(scores$home.score=="NaN" & scores$away.score=="NaN"),]
      }else{
      next.f=read.csv(file=filename, colClasses=c("character"),strip.white=TRUE, stringsAsFactors=FALSE)
      colnames(next.f)=tolower(str_strip.white(names(next.f)))
      if(!all(c("date","home.team","home.score","away.team","away.score") %in% colnames(next.f)))
         stop(paste(filename,"is missing date, home.team, home.score, away.team, or away.score column.\n"), call.=FALSE)
      if(any(!(names(next.f) %in% names(scores)))){
        new.f.names=names(next.f)[!(names(next.f) %in% names(scores))]
        scores[new.f.names]=NA
        scores.warning[[ffile]]=list(f.names.not.in.others=new.f.names)
      }
      if(any(!(names(scores) %in% names(next.f)))){
        new.s.names=names(scores)[!(names(scores) %in% names(next.f))]
        next.f[new.s.names]=NA
        if(length(scores.warning[[ffile]])==0){
          scores.warning[[ffile]]=list(other.names.not.in.f.file=new.s.names)
        }else{ scores.warning[[ffile]][["other.names.not.in.f.file"]]=new.s.names }
      }
      if(any(next.f$home.team=="NA" | next.f$away.team=="NA")){
        next.f=next.f[!(next.f$home.team=="NA" | next.f$away.team=="NA"),,drop=FALSE]
        if(length(scores.warning[[ffile]])==0){
          scores.warning[[ffile]]=list(NA.in.f.file=TRUE)
        }else{ scores.warning[[ffile]][["NA.in.f.file"]]=TRUE }
      }
      bad.team.names=str_detect(next.f$home.team, "Place Flight") | str_detect(next.f$away.team, "Place Flight") |
                    str_detect(tolower(next.f$away.team), "quarterfinal") | str_detect(tolower(next.f$away.team),"semi-final")
      if(any(bad.team.names)){
        next.f=next.f[!bad.team.names,,drop=FALSE]
        if(length(scores.warning[[ffile]])==0){
          scores.warning[[ffile]]=list(non.team.name.in.f.file=TRUE)
        }else{ scores.warning[[ffile]][["non.team.name.in.f.file"]]=TRUE }
      }
      if(na.remove) next.f = next.f[!(next.f$home.score=="NaN" & next.f$away.score=="NaN"),]
      scores=rbind(scores,next.f) 
    } 
  }
  if(length(scores.warning)!=0){
    for(i in names(scores.warning)){
      if(!is.null(scores.warning[[i]][["f.names.not.in.others"]]))
        cat(paste("In file",i,"the following column names are not in some of the other scores files:",paste(scores.warning[[i]][["f.names.not.in.others"]],collapse=", "),"\nThese columns added to other files with value NA.\n"))
      if(!is.null(scores.warning[[i]][["other.names.not.in.f.file"]]))
        cat(paste("File",i," missing the following column names that are in some of the other scores files:",paste(scores.warning[[i]][["other.names.not.in.f.file"]],collapse=", "),"\nThese added with value NA.\n"))
      if(!is.null(scores.warning[[i]][["NA.in.f.file"]]))
        cat(paste("File",i," has NA for some team names.\nThese were removed.\n"))
      if(!is.null(scores.warning[[i]][["non.team.name.in.f.file"]]))
        cat(paste("File",i," has some team names that are not teams (e.g. First in bracket A).\nThese were removed.\n"))
    }
  }
  
  #check that there are no duplicated names
  if(any(duplicated(names(scores))))
    stop("Duplicated column names are not allowed in the match files.  Names are case insensitive and extra space is striped.\n", call.=FALSE)
  
  #Check that the required columns: date, home.team, home.score, away.team, away.score are present
  req.columns=c("date","home.team","home.score","away.team","away.score")
  if(!all(req.columns %in% names(scores)))
    stop("The scores file must include the columns: date, home.team, home.score, away.team, away.score\n", call.=FALSE)
  
  #strip any extra white space from the score data
  for(i in names(scores)) scores[[i]]=str_strip.white(scores[[i]])
  
  #convert dates to R date format
  if(any(is.na(as.Date(scores$date,date.format)))){
    bad.line = which(is.na(as.Date(scores$date,date.format)))
    cat(paste("something wrong with date in line(s)\n",paste(bad.line,collapse=", "),"\nin the match file.\n Dates must be in the following format",format(Sys.Date(),date.format),"unless you specify a different date format.\n match file being returned"))
    return(list(scores=scores, updated=FALSE, ok=FALSE))
  }  
  scores$date=as.Date(scores$date, date.format)
  
  #reorder in increasing date
  scores=scores[order(scores$date),,drop=FALSE]
  
  #convert score to numbers
  scores$home.score=as.numeric(scores$home.score)
  scores$away.score=as.numeric(scores$away.score)
  
  #Check for bad scores
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5){
      if(is.nan(x)) return(TRUE)
      if(is.na(x)) return(FALSE)
      abs(x - round(x)) < tol
    }
  bad.scores=(!sapply(scores$home.score,is.wholenumber) | !sapply(scores$away.score,is.wholenumber))
  if(any(bad.scores)){
    cat("There are non-numeric or non-whole number scores in the score file:\n")
    cat("errors are on lines ")
    cat((1:dim(scores)[1])[bad.scores]); cat("\n")
    cat("Non-numeric scores in scores file. Scores file is being returned\n.")
    return(list(scores=scores, updated=FALSE, ok=FALSE))
  }
  
  if(missing(teams.file)){
    cat("Alert: teams info file was not passed in.\nWill construct one from the scores data frame but teams in the scores file must use a unique name.\n")
    team.data=data.frame(name=sort(unique(c(scores$home.team,scores$away.team))), stringsAsFactors=FALSE)
    #      team.data$alt.name.0=team.resolver$name      
  }else{
    #read in the team info file(s)
    #all the column names will be made lower and striped of extra white space
    team.data=data.frame()
    for(ffile in teams.file){
      filename=ffile
      if(length(team.data)==0){ 
        team.data=read.csv(file=filename, strip.white=TRUE, colClasses = "character", stringsAsFactors=FALSE)
        colnames(team.data)=tolower(str_strip.white(names(team.data)))
      }else{
        next.f=read.csv(file=filename, colClasses=c("character"),strip.white=TRUE, stringsAsFactors=FALSE)
        colnames(next.f)=tolower(str_strip.white(names(next.f)))
        if(any(!(names(next.f) %in% names(team.data))))
          team.data[names(next.f)[!(names(next.f) %in% names(team.data))]]=NA
        if(any(!(names(team.data) %in% names(next.f))))
          next.f[names(team.data)[!(names(team.data) %in% names(next.f))]]=NA
        team.data=rbind(team.data,next.f)
      } 
    }
    #Check that there are no duplicated column names
    if(any(duplicated(names(team.data))))
      stop("Duplicated column names are not allowed in the team files.  Names are case insensitive and extra space is striped.\n", call.=FALSE)
    
    #Check that the required column: name is present
    req.columns=c("name")
    if(!all(req.columns %in% names(team.data)))
      stop("The team file must have a column called \"name\" for the team (display) names.\n", call.=FALSE)
    
    #strip any extra white space from the team data
    for(i in names(team.data)) team.data[[i]]=str_strip.white(team.data[[i]])
    
    #check that scores and teams data frame have unique column names
    if(any(duplicated(c(names(scores),names(team.data)))))
      stop("Column names that appear in both the scores and teams files are not allowed.  Names are case insensitive and extra space is striped.\n", call.=FALSE)
  }
  
  if( missing(team.resolver) | is.null(team.resolver) ){
    cat("Alert: teams resolver was not passed in.\nWill construct one from the team info data frame.\n")
    team.resolver=data.frame()
    if(any("alt.name" %in% str_sub(colnames(team.data),1,8))){
      cc=str_sub(colnames(team.data),1,8)=="alt.name"
      for(i in 1:dim(team.data)[1]){
        alt.names=team.data[i,cc][team.data[i,cc]!=""]
        if(!(team.data$name[i] %in% alt.names)) alt.names=c(team.data$name[i], alt.names)
        for(j in alt.names)
          team.resolver=rbind(team.resolver, data.frame(name=team.data$name[i],alt.name=j,stringsAsFactors=FALSE))
      }
    }else{
      team.resolver=data.frame(name=team.data$name,alt.name=team.data$name,stringsAsFactors=FALSE)
    }
  }else{
    #read in the team resolver file(s)
    #all the column names will be made lower and striped of extra white space
    team.resolver=read.csv(file=team.resolver, strip.white=TRUE, colClasses = "character", stringsAsFactors=FALSE)
    if(dim(team.resolver)[2] != 2) 
      stop("The team resolver file should have only 2 columns: name and alt.name",.call=FALSE)
    colnames(team.resolver)=tolower(str_strip.white(names(team.resolver)))
    if(!all(c("name","alt.name") %in% colnames(team.resolver))) 
      stop("The team resolver file should have only 2 columns called name and alt.name",.call=FALSE)
    
  }
  
  raw.scores=scores
  tmp=resolve.team.names(scores, team.resolver, team.data)
  scores=tmp$scores
  updated=tmp$updated
  team.resolver=tmp$team.resolver
  team.data=tmp$team.data
  ok=tmp$ok
  return(list(scores=scores, raw.scores=raw.scores, teams=team.data, team.resolver=team.resolver, updated=updated, ok=ok))
}