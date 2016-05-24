rank.teams=function(scores=NULL, teams=NULL, 
                    family="poisson", fun="glm",
                    max.date="2100-6-1", min.date="1900-5-1", date.format="%Y-%m-%d",
                    time.weight.eta=0,
                    add=NULL, silent=FALSE, ...){
#g(mu)=log(mu)
#score(i,j) ~ Pois(alpha_i * beta_j); i is self, j is opposing team
#x_i, the self team, a 1 x # teams matrix with 1 in i-th column, 0 otherwise
#x_j, the opposing team, a 1 x # teams matrix with 1 in j-th column, 0 otherwise
#glm(formula = score.i.vs.j ~ factor(attack.i) + factor(defense.j), data=scores, family = family)

#do some error checking on the scores data frame
if(missing(scores)) stop("A scores data frame must be passed in.\n",call.=FALSE)
if(!is.data.frame(scores)) stop("scores must be a data frame.\n",call.=FALSE)
colnames(scores)=tolower(str_strip.white(names(scores)))
if(any(duplicated(names(scores))))
  stop("Duplicated column names are not allowed in the score data frame.  Names are case insensitive and extra space is striped.\n", call.=FALSE)
req.columns=c("date","home.team","home.score","away.team","away.score")
if(!all(req.columns %in% names(scores)))
  stop("The scores data frame must include the columns: date, home.team, home.score, away.team, away.score\n", call.=FALSE)
#no factors allowed
for(i in names(scores)){
  if(is.factor(scores[[i]])) scores[[i]]=as.character(scores[[i]])
}
if(any(is.na(as.Date(scores$date,date.format))))
  stop("The scores dates are not in the same format as date.format.\n",call.=FALSE)
if(!all(fun %in% c("glm","speedglm","glmnet")))
  stop('The fun argument must be either "glm". "speedglm" or "glmnet".\n',call.=FALSE)

#do some error checking on the teams data frame
if(missing(teams)){
  if(!silent) cat("Alert: teams data frame was not passed in. Will attempt to construct one from the scores data frame.You should ensure that teams use only one name in scores data frame.\n")
  teams=data.frame(name=sort(unique(c(scores$home.team,scores$away.team))),stringsAsFactors=FALSE)
}
if(!is.data.frame(teams)) stop("teams must be a data frame.\n",call.=FALSE)
colnames(teams)=tolower(str_strip.white(names(teams)))
if(any(duplicated(names(teams))))
  stop("Duplicated column names are not allowed in the teams data frame.  Names are case insensitive and extra space is striped.\n", call.=FALSE)
req.columns=c("name")
if(!all(req.columns %in% names(teams)))
  stop("The teams data frame must include the column \"name\" which is the team name.\n", call.=FALSE)
#no factors allowed
for(i in names(teams)){
  if(is.factor(teams[[i]])) teams[[i]]=as.character(teams[[i]])
}

#all column names in scores and teams must be unique
if(any(duplicated(c(names(scores),names(teams)))))
  stop("Column names that appear in both the scores and teams data frames are not allowed.  Names are case insensitive and extra space is striped.\n", call.=FALSE)

#Check for mis-entered years
if(is.na(as.Date(max.date, date.format))) stop(paste("max.date must be entered in the following format:",format(Sys.Date(),date.format),"\n or pass in the argument date.format to specify a different date format.\n"))
if(is.na(as.Date(min.date, date.format))) stop(paste("min.date must be entered in the following format:",format(Sys.Date(),date.format),"\n or pass in the argument date.format to specify a different date format.\n"))
#store dates as Date objects
max.date = as.Date(max.date, date.format)
min.date = as.Date(min.date, date.format)

#Error check on the team names in scores and in the team file
display.names=unique(teams$name)
nteams = length(unique(display.names))

ok=TRUE
#Check for bad names
bad.names=scores$home.team[!(scores$home.team %in% display.names)]
if(length(bad.names)!=0){
  cat("The following home team names don't match names in team file:\n")
  cat(as.character(bad.names),sep=", ")
  cat("\nerrors are on lines ")
  cat((1:dim(scores)[1])[!(scores$home.team %in% display.names)]); cat("\n")
  ok=FALSE
}
bad.names=scores$away.team[!(scores$away.team %in% display.names)]
if(length(bad.names)!=0){
  cat("The following away team names don't match names in team file:\n")
  cat(as.character(bad.names),sep=", ")
  cat("\nerrors are on lines ")
  cat((1:dim(scores)[1])[!(scores$away.team %in% display.names)]); cat("\n")
  ok=FALSE
}
if(!ok){ 
  cat("Bad team names in scores file.  Score file being returned.\n")
  return(scores=scores)
}

#determine which teams to include based on any filter arguments the user included
tmp=team.and.score.filters(list(scores=scores, teams=teams),...)
scores=scores[tmp$include.scores,,drop=FALSE]
teams=teams[teams$name %in% tmp$include.teams,,drop=FALSE]

#Check that any names in add appears in scores
#create add.to.formula which is a list with the pair of predictor names for home and away
#if the predictors are different, pair with be home.foo,away.foo (where foo is a predictor name)
#if the predictors are the same (e.g. surface), pair will be foo,foo
#this list says in which column to look for the additional predictor
if(!is.null(add) & !is.character(add))
  stop("add must be a vector of names to add to the model.",call.=FALSE)
scores.names= names(scores)[!(names(scores) %in% c("date","home.team","home.score","away.team","away.score"))]
add.to.formula=list()
for(tmp in add){
  if(paste("home.",tmp,sep="") %in% scores.names){
    if(!(paste("away.",tmp,sep="") %in% scores.names)) stop(paste("home.",tmp," appears in scores file, therefore, away.",tmp," must also appear.\n", sep=""), call.=FALSE)
    if(tmp %in% scores.names) stop(paste("home.",tmp, ", home.",tmp, " and ", tmp, "appear in scores file, therefore the predictor to add is ambiguous.\n", sep=""), call.=FALSE)
    }
  if(paste("away.",tmp,sep="") %in% scores.names){
    if(!(paste("home.",tmp,sep="") %in% scores.names)) stop(paste("away.",tmp," appears in scores file, therefore, home.",tmp," must also appear.\n", sep=""), call.=FALSE)
    }
  if(str_sub(tmp,1,5)=="home.")
    stop(paste(tmp, "should not have the home. prefix.  Predictors that are different for home and away are\nentered as home.pred and away.pred (a pair) in scores file, and only pred is passed into add.to.formula.\n", sep=""), call.=FALSE)
  if(str_sub(tmp,1,5)=="away.")
    stop(paste(tmp, "should not have the away. prefix.  Predictors that are different for home and away are\nentered as home.pred and away.pred (a pair) in scores file, and only pred is passed into add.to.formula.\n", sep=""), call.=FALSE)
#if we got here then the pair home.tmp and away.tmp are in the scores so add as pair
  if(paste("home.",tmp,sep="") %in% scores.names){
  add.to.formula[[tmp]]=c(paste("home.",tmp,sep=""),paste("away.",tmp,sep=""))
}else{
  if(!(tmp %in% scores.names))
    stop(paste(tmp, " is in add but no column of that name in the scores file.\n", sep=""), call.=FALSE)
  add.to.formula[[tmp]]=c(tmp,tmp)  
}
}

#Replace teams with a data.frame with only the teams in scores

#Set up the data to be used for the ranking
#rows with both NA in home.score and away.score are not useful. Remove
scores = scores[!(is.na(scores$home.score) & is.na(scores$away.score)),]
scores.glm = scores

#Set up the level names for the home and away teams
level.names = sort(unique(c(as.character(scores$home.team),as.character(scores$away.team))))
scores.glm$home.team=factor(scores$home.team, levels=level.names)   
scores.glm$away.team=factor(scores$away.team, levels=level.names)   

#subset scores to date range
scores.glm = scores.glm[scores.glm$date>=min.date,,drop=FALSE]
scores.glm = scores.glm[scores.glm$date<=max.date,,drop=FALSE]

# #Set up the factor names for glm
#commented out since cluster code below overrides this
# attack.team=c(as.character(scores.glm$home.team),as.character(scores.glm$away.team))
# defense.team=c(as.character(scores.glm$away.team),as.character(scores.glm$home.team))
# scores.team=c(scores.glm$home.score,scores.glm$away.score)
# attack=factor(attack.team, levels=level.names)
# defend=factor(defense.team, levels=level.names)

# Detect any unconnected clusters and rank them as separate groups
el=cbind(as.character(scores.glm$home.team),as.character(scores.glm$away.team))
g1 = graph.edgelist(el,directed=FALSE)
clus=clusters(g1) #by factor order
clus.names = get.vertex.attribute(g1, "name")

glm.fit=list()
glmer.fit=list()
#set up the residuals column for scores
scores$home.residuals=NA
scores$away.residuals=NA
el.full=cbind(as.character(scores$home.team),as.character(scores$away.team))
date.filter = scores$date>=min.date & scores$date<=max.date

 tmp.fun = function(x,y){ any(x %in% y) }
 for(cluster in 1:clus$no){
   names.in.clus = clus.names[clus$membership == cluster]
   #el is a 2 column matrix with home.team in col 1 and away.team in col 2 for scores.glm
   #rows.clus says which rows in scores.glm are in the cluster
   rows.clus = apply(el,1,tmp.fun,names.in.clus)
   #scores.clus is a scores data.frame just for the cluster
   scores.clus = scores.glm[rows.clus,,drop=FALSE]
   time.diff=as.numeric(max.date-scores.clus$date)
   time.diff=c(time.diff,time.diff) #since away and home
   attack.team=c(as.character(scores.clus$home.team),as.character(scores.clus$away.team))
   defense.team=c(as.character(scores.clus$away.team),as.character(scores.clus$home.team))
   scores.team=c(scores.clus$home.score,scores.clus$away.score)
   level.names = sort(unique(c(as.character(scores.clus$home.team),as.character(scores.clus$away.team))))  
   attack=factor(attack.team, levels=level.names)
   defend=factor(defense.team, levels=level.names)
   
   #set up xlevels and terms.dataClasses if fun != glm since speedglm and glmnet don't return this and coef.fbRanks needs it
   if(fun %in% c("speedglm","glmnet")){
     xlevels = list()
     xlevels$attack = level.names
     xlevels$defend = level.names
     terms.dataClasses = c(scores.team="numeric",attack="factor",defend="factor")
   }

   #Now set up any additional factors for the model
   #tack ".f" to name of additional factor
   add.to.f=""
   for(add.f in names(add.to.formula)){
     ok=TRUE
     #add.to.formula is a list for each name of extra predictor.
     #Each list is 2 values of the column name for home and away predictor
     new.col=c(scores.clus[[add.to.formula[[add.f]][1]]], scores.clus[[add.to.formula[[add.f]][2]]])
     if(is.character(new.col)){ #then it is a factor
       new.col=factor(new.col) #make it a factor
       if(fun %in% c("speedglm","glmnet"))
         xlevels[[paste(add.f,".f",sep="")]]=levels(new.col) #speedglm and glmnet don't return this
       if(length(levels(new.col))==1){
         ok=FALSE #don't add that factor to the fit
       }
     }
     assign(paste(add.f,".f",sep=""), new.col) #add ".f" to name
     #add.to.f is what to add to the model formula
     if(ok){
       add.to.f=paste(add.to.f,"+",paste(add.f,".f",sep=""),sep="")
       if(fun %in% c("speedglm","glmnet"))
         terms.dataClasses[paste(add.f,".f",sep="")]=class(new.col)
     }
   }

   if("glm" %in% fun){
     #set up the formula
     #no intercept in the Dixon and Coles model
     my.formula = paste("scores.team~-1+attack+defend",add.to.f,sep="")
     glm.fit[[paste("cluster.",cluster,sep="")]]=
       glm(formula=as.formula(my.formula), family=family, na.action="na.exclude",weights=exp(-1*time.weight.eta*time.diff))
     #which rows in the full scores, with all dates, are in the glm cluster; rows.clus.full is needed to subset scores
     rows.clus.full = apply(el.full,1,tmp.fun,names.in.clus) & date.filter
     #get residuals
     glm.fit[[paste("cluster.",cluster,sep="")]]
     theresids=residuals(glm.fit[[paste("cluster.",cluster,sep="")]],type="response")
     scores$home.residuals[rows.clus.full]=theresids[1:(length(theresids)/2)]
     scores$away.residuals[rows.clus.full]=theresids[(length(theresids)/2+1):length(theresids)]
   }
   if("speedglm" %in% fun){
     require(speedglm)
     weights=exp(-1*time.weight.eta*time.diff)
     bad.nas = is.na(scores.team)
     if(sum(bad.nas)==length(scores.team)){
       glm.fit[[paste("cluster.",cluster,sep="")]]="all NAs"
       next
     }
     #set up the data frame
     da=data.frame(scores.team=scores.team, attack=attack, defend=defend, weights=weights)
     for(add.f in names(add.to.formula)){       
       da[,paste(add.f,".f",sep="")]=get(paste(add.f,".f",sep=""))
     }  
     da=da[!bad.nas,]
     
     if(family=="poisson") linkfun = poisson(log)
     #set up the formula
     my.formula = paste("scores.team~-1+attack+defend",add.to.f,sep="")
     glm.fit[[paste("cluster.",cluster,sep="")]]=
       speedglm(formula=as.formula(my.formula), data=da, family=linkfun, weights=da$weights)
     #need to add xlevels to speedglm fit since needed for returning coef later
     glm.fit[[paste("cluster.",cluster,sep="")]]$xlevels=xlevels
     attr(glm.fit[[paste("cluster.",cluster,sep="")]]$terms,"dataClasses")=terms.dataClasses
     
     #which rows in the full scores, with all dates, are in the glm cluster; rows.clus.full is needed to subset scores
     rows.clus.full = apply(el.full,1,tmp.fun,names.in.clus) & date.filter
     #residuals() won't work with speedglm object
     attack.coef=defend.coef=rep(0,length(names.in.clus))
     loc=match(paste("attack",names.in.clus,sep=""), names(glm.fit[[paste("cluster.",cluster,sep="")]]$coef))
     attack.coef=glm.fit[[paste("cluster.",cluster,sep="")]]$coef[loc]
     names(attack.coef)=names.in.clus
     attack.coef[is.na(attack.coef)]=0
     loc=match(paste("defend",names.in.clus,sep=""), names(glm.fit[[paste("cluster.",cluster,sep="")]]$coef))
     defend.coef=glm.fit[[paste("cluster.",cluster,sep="")]]$coef[loc]
     names(defend.coef)=names.in.clus
     defend.coef[is.na(defend.coef)]=0
     scores$home.residuals[rows.clus.full]=scores$home.score[rows.clus.full]-exp(attack.coef[match(scores$home.team[rows.clus.full],names.in.clus)]+defend.coef[match(scores$away.team[rows.clus.full],names.in.clus)])
     scores$away.residuals[rows.clus.full]=scores$away.score[rows.clus.full]-exp(attack.coef[match(scores$away.team[rows.clus.full],names.in.clus)]+defend.coef[match(scores$home.team[rows.clus.full],names.in.clus)])     
#      theresids=residuals(glm.fit[[paste("cluster.",cluster,sep="")]],type="response")
#      scores$home.residuals[rows.clus.full]=theresids[1:(length(theresids)/2)]
#      scores$away.residuals[rows.clus.full]=theresids[(length(theresids)/2+1):length(theresids)]
   }
   if("glmnet" %in% fun){
     require(glmnet)
     weights=exp(-1*time.weight.eta*time.diff)
     bad.nas = is.na(scores.team)
     if(sum(bad.nas)==length(scores.team)){
       glm.fit[[paste("cluster.",cluster,sep="")]]="all NAs"
       next
     }
#this code gets rid of teams with 0 goals fielded or 0 goals allowed
     #isn't quite the way to do this because it sets the strength to 0 because those columns are still there
     #need to get rid of those coef altogether
#      bad.attack.teams=bad.defend.teams=c()
#      for(i in level.names){
#        if(sum(scores.team[attack.team==i], na.rm=TRUE)==0) bad.attack.teams=c(bad.attack.teams,i)
#        if(sum(scores.team[defense.team==i], na.rm=TRUE)==0) bad.defend.teams=c(bad.defend.teams, i)       
#      }
    
     #clean out any NAs as glmnet doesn't like those
#     bad.scores = (attack.team %in% bad.attack.teams) | (defense.team %in% bad.defend.teams)
     bad.scores=FALSE
     attack=attack[!bad.nas & !bad.scores]
     defend=defend[!bad.nas & !bad.scores]
     weights = weights[!bad.nas & !bad.scores]
     scores.team=scores.team[!bad.nas & !bad.scores]
     
     #set sparse matrix of explanatory variable. factors are refered to by number
     #this code is implicitly assuming that teams in attack col are same as those in defend col
     n = length(attack) #n games in clus
     sx.i=1:n; sx.j = as.numeric(attack); sx.x = rep(1,n)
     sx.dn = paste("attack",level.names,sep="")
     add.col = length(level.names)
     sx.i=c(sx.i,1:n); sx.j = c(sx.j,as.numeric(defend) + add.col); sx.x = c(sx.x,rep(1,n))
     sx.dn = c(sx.dn,paste("defend",level.names,sep=""))
     #first.factor is keeping track of which columns correspond to first factor and thus need
     #to be set to 0 to ensure identifiability
     first.factor = add.col+1
     #add.col keeps track of what column in x we should start on
     add.col = add.col + length(level.names)
     for(add.f in names(add.to.formula)){
       add.f=paste(add.f,".f",sep="")
       if(is.factor(get(add.f))){
         sx.i=c(sx.i,1:n); sx.j = c(sx.j,as.numeric(get(add.f)[!bad.nas & !bad.scores]) + add.col); sx.x = c(sx.x,rep(1,n))
         sx.dn = c(sx.dn,paste(levels(get(add.f)),add.f,sep=""))
         #first factors must be set to 0; keep track of which cols these appear in
         first.factor = c(first.factor,add.col+1)
         add.col = add.col + length(levels(get(add.f)))
       }else{
         sx.i=c(sx.i,1:n); sx.j = c(sx.j, 1 + add.col); sx.x = c(sx.x,get(add.f)[!bad.nas & !bad.scores])
         sx.dn = c(sx.dn,add.f)
         add.col = add.col + 1
       }
     }
     #set all the first factors (after attack) to 0; for identifiability
     #don't need to set first attack factor to 0 since we exclude the intercept
     #do this by removing any rows in x that have j == the column where the defend coef for team 1 appears
     sx.i = sx.i[!(sx.j %in% first.factor)]
     sx.x = sx.x[!(sx.j %in% first.factor)]
     sx.j = sx.j[!(sx.j %in% first.factor)]
     #add.col is the number of columns in x
     sx=sparseMatrix(i=sx.i, j=sx.j,x=sx.x,dims=c(n,add.col),dimnames=list(1:n,sx.dn))

     call.list = list(x=sx, y=scores.team, intercept=FALSE, family=family, alpha=1, lambda=0, weights=weights)
     glm.fit[[paste("cluster.",cluster,sep="")]]=do.call(glmnet,call.list)
     #add xlevels and terms to fit since coef.fbRanks will use that
     glm.fit[[paste("cluster.",cluster,sep="")]]$xlevels=xlevels
     my.formula = as.formula(paste("scores.team~-1+attack+defend",add.to.f,sep=""))     
     glm.fit[[paste("cluster.",cluster,sep="")]]$terms=terms(my.formula)
     attr(glm.fit[[paste("cluster.",cluster,sep="")]]$terms,"dataClasses")=terms.dataClasses
     
     #which rows in the full scores, with all dates, are in the glm cluster
     #rows.clus.full is needed to subset scores; scores is different than scores.glm since it is
     #not restricted by max.date and min.date
     rows.clus.full = apply(el.full,1,tmp.fun,names.in.clus) & date.filter
     #residuals() won't work with glmnet object
     attack.coef=defend.coef=rep(NA,length(names.in.clus))
     loc=match(paste("attack",names.in.clus,sep=""), rownames(coef(glm.fit[[paste("cluster.",cluster,sep="")]])))
     attack.coef=coef(glm.fit[[paste("cluster.",cluster,sep="")]],s=0)[loc]
     names(attack.coef)=names.in.clus
     loc=match(paste("defend",names.in.clus,sep=""), rownames(coef(glm.fit[[paste("cluster.",cluster,sep="")]])))
     defend.coef=coef(glm.fit[[paste("cluster.",cluster,sep="")]],s=0)[loc]
     names(defend.coef)=names.in.clus
     scores$home.residuals[rows.clus.full]=scores$home.score[rows.clus.full]-exp(attack.coef[match(scores$home.team[rows.clus.full],names.in.clus)]+defend.coef[match(scores$away.team[rows.clus.full],names.in.clus)])
     scores$away.residuals[rows.clus.full]=scores$away.score[rows.clus.full]-exp(attack.coef[match(scores$away.team[rows.clus.full],names.in.clus)]+defend.coef[match(scores$home.team[rows.clus.full],names.in.clus)])     
   }
#    if("glmer" %in% fun){
#      require(lme4)
#      #set up the formula
#      #my.formula = paste("scores.team~(1|attack)+(1|defend)",add.to.f,sep="")
#      my.formula = paste("scores.team~attack+(1|defend)",add.to.f,sep="")
#      if(!intercept) my.formula=paste(my.formula,"-1",sep="")
#      glmer.fit[[paste("cluster.",cluster,sep="")]]=
#        glmer(formula=as.formula(my.formula), family=family, na.action="na.exclude", weights=exp(-1*time.weight.eta*time.diff))
#     #to get effects use ranef(glmer.fit)$attack and ranef(glmer.fit)$defend
#      #which rows in the full scores, with all dates, are in the glm cluster; rows.clus.full is needed to subset scores
#      rows.clus.full = apply(el.full,1,tmp.fun,names.in.clus) & date.filter
#      theresids=residuals(glmer.fit[[paste("cluster.",cluster,sep="")]],type="response")
#      scores$home.residuals[rows.clus.full]=theresids[1:(length(theresids)/2)]
#      scores$away.residuals[rows.clus.full]=theresids[(length(theresids)/2+1):length(theresids)]
#    }
   
 }

# Set up the fbRanks object to return
rtn.list = list(fit=glm.fit, graph=list(graph=g1, membership=clus$membership, csize=clus$csize, no=clus$no, names=get.vertex.attribute(g1, "name")), 
                scores=scores, teams=teams,
                max.date=max.date, min.date=min.date,time.weight.eta=0, date.format=date.format)
class(rtn.list) = "fbRanks"

if(!silent) print(rtn.list)

invisible(rtn.list)
}

