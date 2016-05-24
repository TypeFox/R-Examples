# check input

RRcheck.p <- function(model,p){
  if (model %in% c("Warner","Crosswise") &&
        (length(p)!=1 || p<0 || p>1 || p==0.5)){
    warning(paste0("For model: '",model,"' , the randomization probability 'p' 
                   has to be in [0,1] and different to 0.5"))
  }else if (model %in% c("UQTknown","Kuk","SLD","CDM","UQTunknown")  && 
              (length(p)!=2 || sum(p<0) > 0 || sum(p>1) > 0 || p[1]==p[2]) ){
    warning(paste0("For model: '",model,"' , the randomization probability vector 'p' 
                   must contain two distinct values in [0,1]"))
  }else if (model=="FR" && (sum(p<0) > 0 || sum(p>1) > 0 || sum(p)>=1) ){
    warning("For the Forced Response (FR) model, 'p' must be a vector with a length 
            equal to the number of response categories, containing probabilities 
            in [0,1] which sum up to a value smaller than 1")
  }else if (model =="CDMsym"){
    if(length(p)!=4 || sum(p<0) > 0 || sum(p>1) > 0 ) 
      warning(paste0("For model: '",model,"' , the randomization probability vector 'p' 
                     must contain four distinct values in [0,1]. p[1]/p[2]: forced Yes/No 
                     response for group 1; p[3]/p[4]: forced Yes/No response for group 2."))
    if (p[1]==p[3] || p[2] == p[4])
      warning(paste0("For model: '",model,"' , the randomization probabilities must meet 
                     these conditions: p[1]!=p[3]; p[2]!=p[4] (i.e., the probability of 
                     forced Yes/No responses have to differ for the two groups)"))
    if (p[1]+p[2]>=1 || p[3]+p[4]>=1)
      warning(paste0("For model: '",model,"' , the randomization probabilities must meet 
                     these conditions: p[1]+p[2]<1 and p[3]+p[4]<1 (i.e., the sum of 
                     probabilities for a forced response has to be smaller than 1)"))
  }else if (model== "Mangat" &&  (length(p)!=1 || p<0 || p>1)){
    warning(paste0("For model: '",model,"' , the randomization probability 'p' has to be in [0,1]"))
  }else if (model %in% c("mix.norm") && (length(p)!=3 || p[1]<0 || p[1]>1)){
    warning(paste0("For model: '",model,"' , the randomization probability 'p' has to 
                   be a vector with 3 values: p[1] = p_truth in [0,1] ; 
                   p[2],p[3] = Mean, Variance of masking distribution  "))
  }else if (model %in% c("mix.exp") && (length(p)!=2 || p[1]<0 || p[1]>1)){
    warning(paste0("For model: '",model,"' , the randomization probability 'p' has to 
                   be a vector with 2 values:  p[1] = p_truth in [0,1] ; 
                   p[2] = Mean of masking distribution"))
  }else if (model %in% c("mix.unknown") && (length(p)!=2 || p[1]==p[2] || any(p<0) || any(p>1))){
    warning(paste0("For model: '",model,"' , the randomization probability 'p' has to be a vector 
                   with 2 probabilities between 0 and 1: 
                   p[1] = p_truth[group==1] != p[2] = p_truth[group==2] "))
  }else if(model %in% "custom" && (any(colSums(p) != 1 | nrow(p) != ncol(p)))){
    warning("For RR model 'custom', 'p' must be a quadratic misclassification matrix 
            with columns adding up to one!")
    
  }
}


RRcheck.xp <- function(model,y,p,vectorName){
  RRcheck.p (model,p)
  if (model %in% c("Warner","Crosswise","UQTknown","UQTknown","Mangat","SLD","CDM") &&
        (length(table(y))!=2 | min(y)!=0 | max(y)!= 1) ){
    warning(paste0("For the model: '",model,"' , '",vectorName,"' must be a numerical vector containing values 0 and 1"))
  }
  else if (model=="Kuk" && (max(y!=floor(y)) ||   min(y)<0   || 
                              max(y)<=0 || length(table(y))-1>max(y))){
    warning(paste0("For a single interview, '",vectorName,"' must be a numerical vector containing values 0 and 1 (1 is equivalent to naming a red card, if 'p[1]' and 'p[2]' give the proportions of red cards for participants with or without the sensitive attribute respectively. If the procedure was repeated 'm' times, 'response' contains values from 0 to 'm' (giving the number, how often a subject reported a red card)"))
  } 
  else if (model== "FR" && (  max(y!=floor(y)) || length(table(y))>length(p) ||  
                                min(y)<0 || max(y)> length(p)-1)){
    warning(paste0("For the model: 'FR', '",vectorName,"'  must be a numerical vector containing integers 0,1,...,m-1 ('m'=number of response categories, defined by length of vector 'p')."))        
  }else if(model == "custom" && (max(y!=floor(y)) || length(table(y))>ncol(p) ||  
                                 min(y)<0 || max(y)> ncol(p)-1)){
    warning(paste0("For the model: 'custom', '",vectorName,"'  must be a numerical vector containing integers 0,1,...,m-1 ('m'=number of response categories, defined by dimension of the matrix 'p')."))   
  }
  
}

RRcheck.xpgroup <- function(model,y,p,group,vectorName){
  RRcheck.xp (model,y,p,vectorName)
  if( is2group(model) ){
    group <- as.numeric(group)
    if ( missing(group) || is.null(group) || 
           length(group) != length(y) || !(length(table(group)) %in% 2:3) ||
           !(min(group) %in% 0:1) || max(group) != 2){
      warning(paste0("The model '",model,"' requires a vector 'group' specifying the group membership 
                     of each participant by 1 and 2 (numeric or factor)."))
    }
  }else if ( !  (missing(group) || is.null(group)) ){
    if ( length(group) != length(y) || !(length(table(group)) %in% 1:2) ||
           !(min(group) %in% 0:1) || max(group) != 1){
      warning(paste0("The model '",model,"' requires a vector 'group' specifying the group membership 
                     of each participant by 0 for DQ and 1 for the RR format."))
    }
  }
}

RRcheck.log <- function(model, y, p, group, n.response, vectorName){
#   if (! model  %in% c("FR", "custom")){
#     RRcheck.xp (model,y,p,vectorName)
#   } else 
    if ( any(y!=floor(y)) || length(table(y))>(n.response+1) ||  
                                  min(y)<0 || any(y>n.response)){
    warning(paste0("For a logistic regression with the model 'FR'/'custom', 
                   responses  must be between 0,...,n.response (default: n.response=1)."))        
  }
  RRcheck.p(model,p)
  if ( is2group(model)) {
    group <- as.numeric(group)
    if ( missing(group) || is.null(group) || 
           length(group) != length(y) || !(length(table(group)) %in% 2:3 ) || !(min(group) %in% 0:1) || max(group) != 2){
      warning(paste0("The model '",model,"' requires a vector 'group' specifying the group membership of each participant by 1 and 2 (numeric or factor). Additionally, the vector 'group' may contain the number 0 for participants who were given a direct question."))
    }else if (length(table(group)) == 3 )
      warning("The group vector contains 3 categories, which are interpreted as follows: 0=direct answer ; 1=group 1 with randomization probability p[1] ; 2=group 2 with randomization probability p[2]",call. = F)
  }else if (!missing(group) || !is.null(group)){
    if (length(group) != length(y) || !(length(table(group)) %in% 1:2 ) || !(min(group) %in% 0:1) || max(group) != 1)
      warning("The 'group' vector may contain 1 or 2 categories, which are interpreted as follows: 0=direct question format ; 1=randomized response format")
    if (length(table(group)) == 2 )
      warning("The group vector contains 2 categories, which are interpreted as follows: 0=direct question format ; 1=randomized response format",call. = F)
  }
#   if (length(table(y))!=2){
#     warning("For logistic regression, only dichotomous responses are allowed")
#   }
}

RRcheck.param<- function(pi){
  pi[pi <0] <- 0
  pi[pi >1] <- 1
  return(pi)
}

RRcheck.pi <- function(model,pi,n){
  if (floor(n)!=n || n<=0){
    warning("Sample size 'n' has to be an integer");
  }
  if (model %in% c("custom", "FR")){
    if (length(pi) == 1 ){
      if (pi<=0 | pi>=1)
        warning("The proportion 'pi' of the prevalence of the sensitive attribute 
                (response = 1) must be within the interval (0,1)")
    }else if (sum(pi<0) > 0 || sum(pi>1) > 0 || sum(pi)!=1)
      warning("For the 'FR'/'custom' model, 'pi' must be a vector with a length equal to the number of response categories, containing probabilities in [0,1] which sum up to 1")
  }else if (model %in% c("mix.norm") && (length(pi)!=2 || pi[2]<=0)){
    warning("For the model 'mix.norm', the true state 'pi' of the sensitive attribute is defined as pi[1] = Mean / pi[2] = Variance of the 'true' normal distrubtion.")
  }else if (model %in% c("mix.exp") && (length(pi)!=1 || pi<=0)){
    warning("For the model 'mix.exp', the true state 'pi' of the sensitive attribute is defined as pi[1] = Mean of the 'true' exponential distrubtion.")
  }else if (!(model %in% c("mix.norm","mix.exp","mix.unknown"))){
    if (pi<0 || pi>1){
      warning("True proportion 'pi' has to be in the interval [0,1]");
    }
  }
}

RRcheck.rate <- function(rate){
  if (rate<0 || rate >1 || length(rate) != 1){
    warning("'carriersComplianceRate' and 'nonCarriersCompRate' must be in [0,1]")
  }
}

RRcheck.groupRatio <- function(groupRatio){
  if (groupRatio<0 || groupRatio >1 || length(groupRatio) != 1){
    warning("'groupRatio' must be in (0,1)")
  }
}

RRcheck.cor <- function(X,m,models,p.list,nameVariables,groups){
  if ( any(is.na(models))){
    models.char <- paste(models,collapse=", ")
    warning(paste("\nVector 'models' does not define valid models. Available models: 
                  Warner, Mangat ,Kuk, FR, UQTknown, Crosswise, mix.norm, mix.exp, 
                  mix.unknown, and direct (i.e., no randomization). 
                  Currently defined: ",models.char))
  }
  if (m!=length(models) )
    warning("Number of variables does not match to length of 'models'")
  if ( m!=length(p.list))
    warning("Number of variables does not match to length of 'p.list'")
  for (i in 1:m){
    RRcheck.xp(models[i],X[,i],p.list[[i]],nameVariables[i]) 
  }
  
  numMultiplGr <- sum(is2group(models))
  if ( numMultiplGr > 0){
    if (missing(groups) || is.null(groups)){
      warning("Missing argument: 'groups' needed for multiple group models such as CDM, CDMsym, SLD or UQTunknown")
    }
    if (ncol(as.matrix(groups)) != numMultiplGr ){
      warning(paste0("'groups' not correctly specified: one column for each 
                     of the RR models CDM, CDMsym, UQTunknown and SLD needed"))
    }
    if (  length(table(groups)) != 2 || min(groups) != 1 || max(groups) != 2){
      warning(paste0("'groups' not correctly specified, only 1 and 2 should be used"))
    }
  }
#   else if (!is.null(groups)){
#     warning("'groups' is ignored, no multiple group RR model specified")
#   }
}


RRcheck.start <- function(model,x,start){
  if (is2group(model)){
    if ( (ncol(x)+1) != length(start) || sum(is.na(start))>0){
      start <- NULL
      warning(paste0("The vector 'start' must have a length of ",ncol(x)+1,
                     " for the model ",model,". Starting value are set to 
                     default (the estimation by RRuni)."))
    }
  }else{
    if (ncol(x) != length(start) || sum(is.na(start))>0){
      start <- NULL
      warning(paste0("The vector 'start' must have a length of ",ncol(x),
                     " for the model ",model,". Starting value are set 
                     to default (the estimation by RRuni)."))
    }
  }
  return(start)
}

RRcheck.lin <- function(y, w, u, models, p.list,group){
  if (! is.numeric(y)){
    warning("The dependent variable 'y' has to be numeric.")
  }
  n <- nrow(y)
  if ( nrow(w)!=n){
        warning("The RR predictors in 'w' must have the same length as the dependent variable 'y'")
      }
  if (! (missing(u) || is.null(u)) && nrow(u) != n){
      warning("The nonRR predictors in 'u' must have the same length as the dependent variable 'y'")
    }
if (  length(models)!=ncol(w)) {
    warning("Number of RR predictors in 'w' does not match the number specified by the vector 'models'")
  }
  if (length(models) != length(p.list)){
    warning("Number of RR predictors specified by the vector 'models' does not match the randomization probabilities in 'p.list'")
  }
  ## GROUP TEST
  if ( ncol(group)!=length(models) || min(group)<0 || max(group)>2 || max(group!=floor(group))){
      warning(paste0("'group' must have ",length(models)," columns containing 0 for DQ and 1 or 2 depending on RR group membership"))
    } 
  if ( sum( models %in% c("CDM","CDMsym"))>0){
    warning("The CDM and CDMsym models can not be used in the RR linear regression framework because no appropriate missclassification matrix PW can be constructed")
  }
  for (i in 1:length(models)){
    RRcheck.xpgroup(models[i],w[,i],p.list[[i]],group[,i],colnames(w)[i])
  }
}