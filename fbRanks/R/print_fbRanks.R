print.fbRanks = function(x, ..., scaling="median total", base=2, log.transform.attack=FALSE,header=TRUE,silent=FALSE, type="", file=""){
  if(type=="html") require(xtable)
  
  team.data = x$teams
  names(team.data)=tolower(names(team.data))
  glm.fits=x$fit
  scores=x$scores
  names(scores)=tolower(names(scores))
  
  #Set up the team names; used for printing
  cc = tolower(names(team.data))=="name"
  all.team.names=team.data[,cc]
  nteams = length(all.team.names)
  
  #team.and.scores.filters will also print warnings about improper filter specs
  include.teams=try(team.and.score.filters(list(scores=scores, teams=team.data),...))
  if(class(include.teams)=="try-error"){ return(invisible(list(ranks="nothing to print", raw="nothing to print", original="nothing to print", ok=FALSE)))
  }else{ include.teams=include.teams$include.teams }
  extra=list(...)
  names(extra)=tolower(names(extra))
  extra.team = extra[names(extra) %in% tolower(names(team.data))]
  extra.scores = extra[names(extra) %in% tolower(names(scores))]
  
  #Create a labels for any score filter to add to the team list
  #For example, if venue is filtered on, we add a column to the team with the venue shown
  extra.scores.filter = extra.scores
  for(j in names(extra.scores)){
    cc = names(scores)==j
    if(identical(tolower(extra.scores.filter[[j]]),"all")){
      extra.scores.filter[[j]] = unique(as.vector(as.matrix(scores[,cc])))
    }
  }
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
  }
  
  #Need to get rid of scores data past the fit so that only teams in the appropriate data range appear in include.teams
  scores = scores[scores$date>=x$min.date,,drop=FALSE]
  scores = scores[scores$date<=x$max.date,,drop=FALSE]
  
  #set flag that something is printed
  ok=TRUE
  
  #Print header
  if(!silent) cat("\n")
  if(header & !silent){
    cat(paste("Team Rankings based on matches ", format(x$min.date, x$date.format)," to ",format(x$max.date, x$date.format),"\n",sep=""))
    # if(!identical(orig.team.filter,"all")){ cat("Team filter: "); cat(paste(orig.team.filter,collapse=",")); cat("\n")  }
    for(j in names(extra)){
      cat(j); cat(": "); cat(paste(extra[[j]],collapse=",")); cat("\n")
    }
  }
  
  nclus = length(glm.fits)
  disp.ranks=raw.ranks=original.ranks=list() 
  nothing.to.print=TRUE
  #list of coefficients
  coef.list=coef(x)$coef.list
  
  for(clus in 1:nclus){
    glm.fit=glm.fits[[clus]]
    
    attack.scores = coef.list[[clus]]$attack      
    defense.scores = coef.list[[clus]]$defend
    team.names = glm.fit$xlevels$attack
         
      #compute n games for each team; table returns names sorted alphabetically
      n=table(c(scores$home.team[!is.na(scores$home.score)],
                scores$away.team[!is.na(scores$away.score)]))
      n=n[names(n) %in% team.names]  

    #original ranks coming out of the glm
    orig.ranks=data.frame(
      team=team.names,total=attack.scores-defense.scores,
      attack=attack.scores,
      defense=defense.scores,
      n.games=n, 
      row.names=NULL)
    
    #scale the scores based on user's choice
    if(scaling == "median total") interc=median(attack.scores-1*defense.scores,na.rm=TRUE)/2
    if(scaling == "median attack") interc=median(attack.scores,na.rm=TRUE)
    if(scaling == "median defend") interc=median(-1*defense.scores,na.rm=TRUE)
    if(is.numeric(scaling)) interc=scaling
    attack.scores = attack.scores-interc
    defense.scores = -1*defense.scores-interc
    #Stick bad scores at the bottom but list with the biggest attack or defense scores higher
    attack.scores[is.na(attack.scores) | abs(attack.scores)>10]=-1000*abs(min(attack.scores,na.rm=TRUE))
    defense.scores[is.na(defense.scores) | abs(defense.scores)>10]=-1000*abs(min(defense.scores,na.rm=TRUE))
    total=(attack.scores + defense.scores)/log(base)
    total= sort(total, index.return=TRUE, decreasing=TRUE)
    
    #set up the ranks data.frame which will be printed
    team.ranks=data.frame(
      team=team.names[total$ix],total=total$x,
      attack=attack.scores[total$ix],
      defense=defense.scores[total$ix],
      n.games=n[total$ix], 
      row.names=NULL)
    
    #add the extra columns requested in the print call
    for(j in names(extra.scores)){
      team.ranks[j]=unlist(label.extra.scores[[j]][team.names])[total$ix]
    }
    for(j in names(extra.team)){
      team.ranks[j]=team.data[[j]][match(team.names,all.team.names)][total$ix]
    }
    
    #only return the teams in the include.teams vector
    team.ranks=team.ranks[team.ranks$team %in% include.teams,,drop=FALSE]
    
    if(dim(team.ranks)[1]>0){ #something to print
      nothing.to.print=FALSE
      row.names(team.ranks)=1:dim(team.ranks)[1]
      
      #display.team.ranks is prettied
      bad.est = abs(team.ranks$attack)>10
      display.team.ranks=team.ranks
      display.team.ranks$attack[bad.est]=NA
      bad.est = abs(team.ranks$defense)>10
      display.team.ranks$defense[bad.est]=NA
      bad.est = abs(team.ranks$attack)>10 | abs(team.ranks$defense)>10
      display.team.ranks$total=round(display.team.ranks$total,digits=2)
      display.team.ranks$total[bad.est]=NA
      if(log.transform.attack){
        display.team.ranks$attack=round(display.team.ranks$attack,digits=2)
        display.team.ranks$defense=round(display.team.ranks$defense,digits=2)
      }else{
        display.team.ranks$attack=round(exp(display.team.ranks$attack),digits=2)
        display.team.ranks$defense=round(exp(display.team.ranks$defense),digits=2)
      }
      
      if(dim(display.team.ranks)[1]!=0){
        row.names(display.team.ranks)=1:dim(display.team.ranks)[1]
        if(!silent){
          if(nclus>1 & !(all(include.teams %in% team.ranks$team))) cat(paste("\nRankings for cluster ",clus,"\n",sep=""))
          print(format(display.team.ranks,scientific=FALSE),right=FALSE)
        }
        if(type=="html"){
          clus.head=""
          if(nclus>1 & !(all(include.teams %in% team.ranks$team))) clus.head=paste("\nRankings for cluster ",clus,"\n",sep="")
          rank.html=print(xtable(display.team.ranks), type = "html", print.results=FALSE)
          rank.html=c(clus.head,rank.html,"<br><br>")
        }
      }
      disp.ranks[[clus]] = display.team.ranks
      raw.ranks[[clus]] = team.ranks
      original.ranks[[clus]]=orig.ranks
    }else{
      disp.ranks[[clus]] = "no teams in this cluster fit the filters"
      raw.ranks[[clus]] = "no teams in this cluster fit the filters"
      original.ranks[[clus]] = "no teams in this cluster fit the filters"
    }
  } #end over list of glm fits
  if(nothing.to.print){
    cat("\nNo teams fit the filters so nothing is printed.\n")
    ok=FALSE
  }
  
  if(type=="html"){
    head1.txt=c("<html>","<body>",paste("<h2>Strength Ratings</h2>",sep=""))
    head2.txt=paste("<h3>","Ratings based on matches from ",x$min.date," to ",x$max.date,"</h3>",sep="")
    footer1.txt=c("\n</body>\n</html>")
    cat(head1.txt,head2.txt,rank.html,"total=log(attack)+log(defense)",footer1.txt,file=file)  
  }
  
  if(nclus==1){ disp.ranks=disp.ranks[[1]]; raw.ranks=raw.ranks[[1]];original.ranks=original.ranks[[1]] }
  invisible(list(ranks=disp.ranks, raw=raw.ranks, original=original.ranks, ok=ok))
}