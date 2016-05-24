###############################################################################
## abn-internal.R --- 
## Author          : Fraser Lewis
## Last modified   : 03/08/2012
###############################################################################

########################################################################################################################
## a set of simple commonsense validity checks on the data.df and data.dists arguments
check.valid.data <- function(data.df=NULL,data.dists=NULL,group.var=NULL){

    ## check data is in a data.frame
    if(!is.data.frame(data.df)){stop("The data must be in a data.frame");}

    ## check data for missing values
    if(sum(complete.cases(data.df))!=dim(data.df)[1]){stop("The data contains missing values! These must be removed.");}
 
    ## check that distributions are in a list  
    if(!is.list(data.dists)){stop("data.dist must be a list");}
   
    if(!is.null(group.var)){## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
     data.df<-data.df[,-which(names(data.df)==group.var)];}

     if(length(names(data.dists))!=length(names(data.df))){stop("data.dists must have named entries");} 
    
     ## check names in list are in correct order and exact match to data.frame
     for(i in 1:dim(data.df)[2]){if(names(data.dists)[i]!=names(data.df)[i]){stop("names in list must match names in data.frame");}}
     
     ## check names of likelihood function are valid
     allowed.dists<-c("gaussian","binomial","poisson");## n.b. these must match inla() family=""
     for(i in 1:length(data.dists)){
                        if(length(which(allowed.dists%in%data.dists[[i]]))!=1){## each variable must have one of the allowed.dists
                        message<-paste("Each variable must have a distribution of either",paste(allowed.dists,collapse=", "));
                        cat(message,"\n");
                        stop("");}
     }

    binomial.vars.indexes<-NULL;
    poisson.vars.indexes<-NULL;
    gaussian.vars.indexes<-NULL;
    
    ## check that data is consistent with distribution given
    for(i in 1:dim(data.df)[2]){## for each variable
                        cur.var<-data.df[,i];
                        if(data.dists[[i]]=="gaussian"){
                            if(is.factor(cur.var)){ cat((names(data.df)[i]),"is invalid - it must not be a factor.\n");stop("");}
                            if(length(unique(cur.var))<=2){ cat((names(data.df)[i]),"is invalid as it has two or less unique values!\n");stop("");}
                        gaussian.vars.indexes<-c(gaussian.vars.indexes,i); 
                        }
                        if(data.dists[[i]]=="binomial"){
                            if(!is.factor(cur.var)){ cat((names(data.df)[i]),"is invalid - it must be a factor\n");stop("");}
                            if(length(unique(cur.var))!=2){ cat((names(data.df)[i]),"is invalid as it must be binary. Multi-category variables should be split into separate binary variables.\n");stop("");}
                        binomial.vars.indexes<-c(binomial.vars.indexes,i); 
                        }
                      
                        if(data.dists[[i]]=="poisson"){
                            if(is.factor(cur.var)){ cat((names(data.df)[i]),"is invalid - it must not be a factor\n");stop("");} 
                            if(length(unique(cur.var))<=2){ cat((names(data.df)[i]),"is invalid as it has two or less unique values!");stop("");}
                        poisson.vars.indexes<-c(poisson.vars.indexes,i);
                        }
           }            
      
return(list(gaus=gaussian.vars.indexes,bin=binomial.vars.indexes,pois=poisson.vars.indexes)); ## return the indexes of any binary variables

} #end of check.valid.data()
#########################################   
## a set of simple commonsense validity checks on the directed acyclic graph definition matrix
check.valid.dag <- function(dag.m=NULL,data.df=NULL, is.ban.matrix=NULL,group.var=NULL){
   
    if(!is.null(group.var)){## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    data.df<-data.df[,-which(names(data.df)==group.var)];}


    ## if dag.m null then create unlimited - empty - network
    if(is.null(dag.m)){## want ban matrix 
                       dag.m<-matrix(rep(0,dim(data.df)[2]^2),ncol=dim(data.df)[2]);
                       colnames(dag.m)<-rownames(dag.m)<-names(data.df); ##names must be set
    return(dag.m); ## finished
    }

    ## check dag is in a matrix
    if(!is.matrix(dag.m)){stop("The DAG definition dag.m must be in a matrix");}

    ## check data for missing names
    if(is.null(colnames(dag.m)) || is.null(rownames(dag.m))){stop("dag.m must have both row and column names set");}

    ## check dimension
    if(dim(dag.m)[1]!=dim(data.df)[2] || dim(dag.m)[2]!=dim(data.df)[2] ){stop("dag.m as dimension inconsistent with data.df");}

    ## check binary
    for(i in 1:dim(dag.m)[1]){for(j in 1:dim(dag.m)[2]){if(dag.m[i,j]!=0 && dag.m[i,j]!=1){stop("dag.m must comprise only 1's or 0's");}}}
 
    ## check diagnonal and cycles - but ignore these checks for a ban matrices
    if(!is.ban.matrix){
    
    for(i in 1:dim(dag.m)[1]){if(dag.m[i,i]==1){stop("dag.m is not a valid DAG - a child cannot be its own parent!");}}
    
    ## coerce to int for sending to C
    dim<-dim(dag.m)[1];## number of cols (or rows)
    dag.m<-as.integer(dag.m);## this creates one long vector - filled by cols from dag.m = same as internal C reprentation so fine.

    res<-.Call("checkforcycles",dag.m,dim
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              ) 
    }
  

}


#########################################   
## a set of simple checks on the list given as parent limits
check.valid.parents<-function(data.df=NULL,max.parents=NULL,group.var)
{ 
  if(!is.null(group.var)){## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    data.df<-data.df[,-which(names(data.df)==group.var)];}
 #print(data.df);print(max.parents);
 ## if a constant then make integer vector
 if(is.numeric(max.parents) && length(max.parents)==1){return(as.integer(rep(max.parents,dim(data.df)[2])));}

 ## if a list must be named list with names as in original data
 if(is.list(max.parents) && length(max.parents)==dim(data.df)[2]){
 for(i in 1:dim(data.df)[2]){if(names(max.parents)[i]!=names(data.df)[i]){stop("names in max.parents list must match names in data.frame data.df");}}
 if(!is.numeric(unlist(max.parents))){stop("max.parents is not valid - must be numeric");}
 max.parents.int<-unlist(max.parents);
 if(length(max.parents.int)!=dim(data.df)[2]){stop("max.parents list is wrong length");}
 max.parents.int<-as.integer(max.parents.int); 
 return(max.parents.int);
 }
 
 stop("max.parents is not valid");

}

#########################################   
## a set of simple checks on the list given as parent limits
check.which.valid.nodes<-function(data.df=NULL,which.nodes=NULL,group.var)
{
  if(!is.null(group.var)){## have a grouping variable so temporarily drop this from data.df - LOCAL TO THIS FUNCTION ONLY
    data.df<-data.df[,-which(names(data.df)==group.var)];}

 if(is.null(which.nodes)){which.nodes<-1:dim(data.df)[2];return(as.integer(which.nodes));} ## if null then assume ALL nodes

 if(is.numeric(which.nodes) && max(which.nodes)<=dim(data.df)[2] && min(which.nodes)>=1 && length(which.nodes)<=dim(data.df)[2]){
    return(as.integer(which.nodes));} else{stop("which.nodes is invalid");}          

}

#########################################   
## a simple check on the grouping variable
check.valid.groups<-function(group.var,data.df,cor.vars){
 
  if(is.null(group.var)){return(list(data.df=data.df,grouped.vars=as.integer(c(-1)),group.ids=as.integer(rep(-1,dim(data.df)[1]))));}## have no groups so just return dummy values 
  if(!(is.character(group.var) && (length(group.var)==1))){stop("name of group variable is not a character?!");}
  if(!length(which(group.var%in%names(data.df)==TRUE))){stop("name of group variable does not match any of those in data.df");}
  group.var.vals<-data.df[,group.var];## get group id data
  data.df<-data.df[,-which(names(data.df)==group.var)];## drop the group variable from original data.frame and overwrite
  
                 
  ## have groups so some checks 
  
       if(    is.factor(group.var.vals) 
           && length(group.var.vals)==dim(data.df)[1]
           && length(unique(group.var.vals))>1 
         ){ ## is factor and of correct length and at least two groups
        } else {stop("grouping variable must be: i) a factor; ii) same length as data.df; and iii) contain more than one group");}
  
  group.var<-as.integer(group.var.vals);## get group memberships in terms of ints

  ## now find out which variables the grouping is to be applied to
  var.noms<-names(data.df);
  if(length(which(cor.vars%in%var.noms==TRUE))!=length(cor.vars)){stop("variables in cor.vars do not match those in data.df");}

  if(max(table(cor.vars))>1){stop("have repeated variables in cor.vars!");}
 
 ## to get to here group.var must be ok and also variable names so return integer code for the variables
 cor.var.indexes<-as.integer(sort(match(cor.vars,var.noms)));## get the index in names(data.df) for each variable and then sort into order
                                                               

return(list(data.df=data.df,grouped.vars=cor.var.indexes,group.ids=group.var));

}

##########################################
## create ordered vector with integers denoting the distribution
##########################################
get.var.types<-function(data.dists=NULL)
{
 store<-rep(NA,length(data.dists));

 for(i in 1:length(data.dists)){
  if(data.dists[[i]]=="binomial"){store[i]<-1;}
  if(data.dists[[i]]=="gaussian"){store[i]<-2;}
  if(data.dists[[i]]=="poisson"){ store[i]<-3;}
 }
 
 return(store);

}
##########################################
## tidy up cache
##########################################
tidy.cache<-function(thecache){
 if(!is.null(thecache[["error.indexes"]])){
             error.combs<-thecache[["error.indexes"]];
             corrected<-list();
             corrected[["children"]]<-thecache[["children"]][-error.combs];
             corrected[["node.defn"]]<-thecache[["node.defn"]][-error.combs,];
             corrected[["mlik"]]<-thecache[["mlik"]][-error.combs];
             return(corrected);## return new cache with appropriate nodes removed 
 }

}
##########################################
## set up graphic windows 
##########################################
graphicsetup<-function(){

## 1. check if an existing Cairo window is available
if(length(which(names(dev.list())=="Cairo"))==0){## no Cairo so open new windows 
           if(.Platform$OS.type=="unix"){mytype<-"x11";} else {mytype<-"win";}
           Cairo(type=mytype);## open up new Cairo window
           #X11();
} else {## have an existing Cairo window so use that 
     dev.set(dev.list()[which(names(dev.list())=="Cairo")[1]]);##set active device to the first available Cairo
     #dev.set(dev.list()[which(names(dev.list())=="X11")[1]])
}

}
# #########################################
# # adjust for split variables
# # this is used after creating all the parent combinations to remove all those combinations which are not allowed
# # e.g. a three cat variable split into X1, X2, X3 can only have one of X1,X2,X3 =1 in any parent combination
# #########################################
#adjust.for.split.variables<-function(split.vars.vec,defn.res.list)
#{
#}
