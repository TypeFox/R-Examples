################################################################
## Predict matches using a fbRanks model
# User can pass in newdata or specify dates and other filter info to construct newdata from the scores data.frame
# in the fbRanks object
################################################################
predict.fbRanks=function(object, ..., newdata=list(home.team="foo", away.team="bar"),
                         max.date="2100-6-1", min.date="1900-5-1", rnd=TRUE, silent=FALSE, show.matches=TRUE, verbose=FALSE,
                         remove.outliers=TRUE, n=100000){
  x=object
  the.fits=x$fit
  clusters=x$graph
  team.data=x$teams
  cc = tolower(names(team.data))=="name"
  all.team.names=team.data[,cc]
  nteams = length(all.team.names)
  
  
  #Check for mis-entered years
  if(missing(max.date)) max.date=as.Date(max.date) #if missing, then the default is in default data format
  else max.date=as.Date(max.date, x$date.format)
  if(is.na(max.date)) stop(paste("max.date must be entered in the following format:",format(Sys.Date(),x$date.format),"\n"))
  if(missing(min.date)) min.date=as.Date(min.date) #if missing, then the default is in default data format
  else min.date=as.Date(min.date, x$date.format)
  if(is.na(min.date)) stop(paste("min.date must be entered in the following format:",format(Sys.Date(),x$date.format),"\n"))
  
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
  
  if(dim(scores)[1]==0) stop("No matches to predict. Either you didn't specify home and away teams or the scores database \n  in your fbRanks object has no matches for your specified dates and filter values.")
  el=cbind(as.character(scores$home.team),as.character(scores$away.team))
  
  tmp.fun = function(x,y){ any(x %in% y) }
  if(!silent){
    cat("Predicted Match Results for ")
    cat(format(min.date,x$date.format)); cat(" to ");  cat(format(max.date,x$date.format));  cat("\n")
    cat("Model based on data from ")
    cat(format(x$min.date,x$date.format)); cat(" to ");  cat(format(x$max.date,x$date.format));  cat("\n---------------------------------------------\n")
  }
  
  #the list of coefficients for each cluster
  coef.list = coef(object)$coef.list
  for(clus in 1:length(the.fits)){
    fit=the.fits[[clus]]
    names.in.clus = clusters$names[clusters$membership == clus]
    rows.clus = apply(el,1,tmp.fun,names.in.clus)
    if(!any(rows.clus)) next
    scores.clus = scores[rows.clus,,drop=FALSE]
    
    #predictor variables other than required
    #How names work; in fitted object predictors (except attack and defend) have the name xyz.f
    #In the scores file, they have the name xyz if applied to both home and away
    #and home.xyz if only for home and away.xyz if for away.  Must be a pair home. and away. in this case
    pred.name=attr(fit$terms, "term.labels")[!(attr(fit$terms, "term.labels") %in% c("attack","defend"))]
    s.pred.name=home.away.predictors=both.predictors=character(0)
    if(length(pred.name)!=0){ #there are other predictors
      s.pred.name=pred.name
      #the name in the scores file; I put .f on the predictor in rank.teams() to make sure name is distinct
      s.pred.name[str_sub(pred.name,-2)==".f"]=str_sub(pred.name[str_sub(pred.name,-2)==".f"],end=-3)
      scores.predictor.names=names(scores.clus)[!(names(scores.clus) %in% c("date","home.team","home.score","away.team","away.score"))]
      home.away.predictors = scores.predictor.names[str_sub(scores.predictor.names,1,5)=="home." | str_sub(scores.predictor.names,1,5)=="away."]
      both.predictors = scores.predictor.names[!(scores.predictor.names %in% home.away.predictors)]
      home.away.predictors = unique(str_sub(home.away.predictors,6))
      bad.pred.name = s.pred.name[!(s.pred.name %in% both.predictors | s.pred.name %in% home.away.predictors)]
      if(length(bad.pred.name)!=0){
          stop(paste("The predictor",paste(bad.pred.name,collapse=", "),"is missing from newdata argument.\n It needs to be specified because it is used in the model.\n"))
      }
      for(ha.pred in s.pred.name[s.pred.name %in% home.away.predictors]){
        tmp.name = paste("home.",ha.pred,sep="")
        if(!(tmp.name %in% names(scores.clus)))
          stop(paste(tmp.name,"if missing from the newdata argument and it is required since it is in the model"),call.=FALSE)
        tmp.name = paste("away.",ha.pred,sep="")
        if(!(tmp.name %in% names(scores.clus)))
          stop(paste(tmp.name,"if missing from the newdata argument and it is required since it is in the model"),call.=FALSE)    
      }
    }
    
    #newdata1 is home as attack team
    newdata1= data.frame(attack=factor(scores.clus$home.team,levels=fit$xlevels$attack),
                         defend=factor(scores.clus$away.team,levels=fit$xlevels$defend))
    #How names work; in fitted object predictors (except attack and defend) have the name xyz.f
    #In the scores file, they have the name xyz if applied to both home and away
    #and home.xyz if only for home and away.xyz if for away.  Must be a pair home. and away. in this case
    for(ii in s.pred.name){ #this is the name as it appears in scores file
      i = which(ii == s.pred.name)
      #predictor name is different if it is a home/away predictor
      if(s.pred.name[i] %in% home.away.predictors) s.name = paste("home.",s.pred.name[i],sep="")
      else s.name = s.pred.name[i]
      if(pred.name[i] %in% names(fit$xlevels)){ #if it is a factor
        newdata1=cbind(newdata1,factor(scores.clus[[s.name]],levels=fit$xlevels[[pred.name[i]]]))
      }else{ newdata1=cbind(newdata1,scores.clus[[s.name]]) }#not a factor
      colnames(newdata1)[dim(newdata1)[2]]=pred.name[i]
    }
    
    #newdata2 is away as attack team
    newdata2= data.frame(defend=factor(scores.clus$home.team,levels=fit$xlevels$defend),
                         attack=factor(scores.clus$away.team,levels=fit$xlevels$attack))
    for(ii in s.pred.name){ #this is the name as it appears in scores file
      i = which(ii == s.pred.name)
      #predictor name is different if it is a home/away predictor
      if(s.pred.name[i] %in% home.away.predictors) s.name = paste("away.",s.pred.name[i],sep="")
      else s.name = s.pred.name[i]
      if(pred.name[i] %in% names(fit$xlevels)){ #if it is a factor
        newdata2=cbind(newdata2,factor(scores.clus[[s.name]],levels=fit$xlevels[[pred.name[i]]]))
      }else{ newdata2=cbind(newdata2,scores.clus[[s.name]]) }#not a factor
      colnames(newdata2)[dim(newdata2)[2]]=pred.name[i]
    }
    
#     #Set up the factor coef so that I can reference them later
#     for(coef.name in names(fit$xlevels)){
#       coef.scores = coef.list[[clus]][[coef.name]]
#       assign(paste(coef.name,".scores",sep=""),coef.scores)
#     }
#     #next lines are so the package passes the tests during build and doesn't complain that attack.scores and defend.scores do not exist
    attack.scores=coef.list[[clus]]$attack
    defend.scores=coef.list[[clus]]$defend
    
    
    #I want to exclude predictions when the attack or defense score does not fit the normality assumption
    bad.attack=fit$xlevels$attack[!detect.normality.outliers(attack.scores)]
    bad.defend=fit$xlevels$defend[!detect.normality.outliers(defend.scores)]
    
    if(remove.outliers){
      newdata1$attack[newdata1$attack %in% bad.attack]=NA
      newdata1$defend[newdata1$defend %in% bad.defend]=NA
      newdata2$attack[newdata2$attack %in% bad.attack]=NA
      newdata2$defend[newdata2$defend %in% bad.defend]=NA
    }
    
    #Create predictions
    prate=0
    for(i in attr(terms(fit),"term.labels")){
      prate=prate+coef(object)$coef.list[[clus]][[i]][newdata1[[i]]]
    }
    home.score=exp(prate)
     #simulate n home.goals
    home.goals = matrix(rpois(n*dim(newdata1)[1],exp(prate)),dim(newdata1)[1],n,byrow=FALSE)
    rownames(home.goals)=newdata1$attack
    
    prate=0
    for(i in attr(terms(fit),"term.labels")){
      prate=prate+coef(object)$coef.list[[clus]][[i]][newdata2[[i]]]
    }
    away.score=exp(prate)
    #simulate n away goals
    away.goals = matrix(rpois(n*dim(newdata2)[1],exp(prate)),dim(newdata2)[1],n,byrow=FALSE)
    rownames(away.goals)=newdata2$attack
    
    #Store attack and defend strengths in the scores that are returned for information
    home.attack=attack.scores[match(newdata1$attack,fit$xlevels$attack)]
    home.defend=defend.scores[match(newdata2$defend,fit$xlevels$attack)]
    away.attack=attack.scores[match(newdata2$attack,fit$xlevels$attack)]
    away.defend=defend.scores[match(newdata1$defend,fit$xlevels$attack)]
        
    #which rows in the scores (to predict) are in the glm cluster
    scores$pred.home.score[rows.clus]=home.score
    scores$pred.away.score[rows.clus]=away.score
    scores$home.residuals[rows.clus]=scores$home.score[rows.clus]-home.score
    scores$away.residuals[rows.clus]=scores$away.score[rows.clus]-away.score
    scores$home.attack[rows.clus]=home.attack
    scores$home.defend[rows.clus]=home.defend
    scores$away.attack[rows.clus]=away.attack
    scores$away.defend[rows.clus]=away.defend
       
    home.win=100*apply(home.goals>away.goals,1,sum)/n
    away.win=100*apply(away.goals>home.goals,1,sum)/n
    tie=100-home.win-away.win
    home.shutout=100*apply(home.goals==0,1,sum)/n
    away.shutout=100*apply(away.goals==0,1,sum)/n
    scores$home.win[rows.clus]=home.win
    scores$away.win[rows.clus]=away.win
    scores$tie[rows.clus]=tie
    scores$home.shutout[rows.clus]=home.shutout
    scores$away.shutout[rows.clus]=away.shutout
    
  } #end clus
  home.goals.sum = c(sum(scores$pred.home.score[!is.na(scores$home.score)],na.rm=TRUE),sum(scores$home.score[!is.na(scores$pred.home.score)],na.rm=TRUE))
  away.goals.sum = c(sum(scores$pred.away.score[!is.na(scores$away.score)],na.rm=TRUE),sum(scores$away.score[!is.na(scores$pred.away.score)],na.rm=TRUE))
  if(!silent){
    exclude.game=is.na(scores$pred.home.score)|is.na(scores$pred.away.score)|is.na(scores$home.score)|is.na(scores$away.score)
    n.games=sum(!exclude.game)
    pred.home=scores$pred.home.score[!exclude.game]
    pred.away=scores$pred.away.score[!exclude.game]
    pred.home.wins=pred.home>pred.away
    home.wins=scores$home.score[!exclude.game]>scores$away.score[!exclude.game]
    pred.away.wins=pred.home<pred.away
    away.wins=scores$home.score[!exclude.game]<scores$away.score[!exclude.game]
    pred.ties=abs(pred.home-pred.away)<.5 #roughly where prob of win > 50%
    pred.home.wins[pred.ties]=FALSE
    pred.away.wins[pred.ties]=FALSE
    ties=scores$home.score[!exclude.game]==scores$away.score[!exclude.game]
    shutout=as.numeric(scores$home.score[!exclude.game]==0) + as.numeric(scores$away.score[!exclude.game]==0)

    if(verbose){
      cat("home goals predicted vs actual: ");cat(round(home.goals.sum,digits=1));cat("\n")
      cat("away goals predicted vs actual: ");cat(round(away.goals.sum,digits=1));cat("\n")
      cat("home wins (predicted vs actual): ");cat(c(round(sum(scores$home.win[!exclude.game])/100,digits=1),sum(home.wins)));cat("\n")
      cat("pred ties (predicted vs actual): ");cat(c(round(sum(scores$tie[!exclude.game])/100,digits=1),sum(ties)));cat("\n")
      cat("away wins (predicted vs actual): ");cat(c(round(sum(scores$away.win[!exclude.game])/100,digits=1),sum(away.wins)));cat("\n")
      cat("pred shutouts (predicted vs actual): ");cat(c(round(sum(scores$home.shutout[!exclude.game])/100,digits=1)+round(sum(scores$away.shutout[!exclude.game])/100,digits=1),sum(shutout)));cat("\n")
      cat("correct predictions: "); cat(sum(pred.home.wins+home.wins==2 | pred.ties+ties==2 | pred.away.wins+away.wins==2)/
        sum(!exclude.game));cat("\n")
      cat("correct home wins: "); cat("fraction: ");cat(sum(pred.home.wins+home.wins==2)/
        sum(pred.home.wins));cat(" "); cat("number: "); cat(sum(pred.home.wins+home.wins==2)/
        sum(home.wins)); cat("\n")
      cat("correct away wins: "); cat(sum(pred.away.wins+away.wins==2)/
        sum(pred.away.wins));cat(" "); cat(sum(pred.away.wins+away.wins==2)/
        sum(away.wins));cat("\n")
      cat("correct ties: "); cat(sum(pred.ties+ties==2)/
        sum(pred.ties));cat(" "); cat(sum(pred.ties+ties==2)/
        sum(ties));cat("\n") 
    }
  }
  
  if(!silent & show.matches){
    for(i in 1:dim(scores)[1]){
      if(rnd){
        if(!is.na(scores$date[i])) cat(format(scores$date[i],x$date.format));cat(" ")
        cat(as.character(scores$home.team[i])); cat(" vs "); cat(as.character(scores$away.team[i])); 
        cat(", HW ");cat(round(scores$home.win[i]));cat("%, AW ")
        cat(round(scores$away.win[i]));cat("%, T "); cat(round(scores$tie[i])); cat("%, pred score ")
        cat(round(scores$pred.home.score[i],digits=1)); cat("-");cat(round(scores$pred.away.score[i],digits=1)); cat("");
        if(!is.nan(scores$home.score[i]) & !is.nan(scores$away.score[i])){
          cat("  actual: ")
          if(scores$home.score[i]>scores$away.score[i]) cat("HW")
          if(scores$home.score[i]<scores$away.score[i]) cat("AW")
          if(scores$home.score[i]==scores$away.score[i]) cat("T")
          cat(" ("); cat(scores$home.score[i])
          cat("-");cat(scores$away.score[i]); cat(")")
        }
        cat("\n")
      }else{
        cat(as.character(scores.clus$home.team[i]));cat(" ");cat(format(home.score[i],digits=2))
        cat(" - ");cat(format(away.score[i],digits=2)); cat(" "); cat(as.character(scores.clus$away.team[i])); cat("\n")
      }
    }
  }
  
  invisible(list(scores=scores, home.score=home.score, away.score=away.score, home.goals.sum=home.goals.sum, home.goals=home.goals, away.goals.sum=away.goals.sum, away.goals=away.goals))
}