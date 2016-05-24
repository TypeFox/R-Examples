smd.c <- function(Group.T=NULL, Group.C=NULL, Mean.T=NULL, Mean.C=NULL, s.C=NULL, n.C=NULL, Unbiased=FALSE)
{
# Function to calculate the standardized mean difference.
if(length(Group.T)>=1 & length(Group.C)>=1 & (!is.null(Mean.T)|!is.null(Mean.C)|!is.null(s.C))) stop("Since you've specified raw data, you do not need to specify summary values.", call.=FALSE)
if(length(Group.T)==1 | length(Group.C)==1) stop("You only have 1 individual in one or both groups, you cannot calculate a variance in such a situation.", call.=FALSE)
if(length(Mean.T)==1 & length(Mean.C)==1 & !is.null(s.C) & is.null(n.C) & Unbiased==TRUE) stop("You need to specify the control group sample size in order to obtain the unbiased estiamte of the standardized mean difference.", call.=FALSE)

# Here determine 'smd.c' on summary data.
if(length(Mean.T)==1 & length(Mean.C)==1)
{
if(!is.null(Group.T) | !is.null(Group.C)) stop("Group.C and Group.T should be NULL, since you've specified the group means directly. Alternatively, you can specify the groups directly and make Mean.C and Mean.T NULL.", call.=FALSE)
d <- (Mean.T - Mean.C)/s.C
if(Unbiased==TRUE) d <- d*gamma((n.C-1)/2)/(sqrt((n.C-1)/2)*gamma(((n.C-1)-1)/2))
return(d)   
}

# Here determine 'smd' for raw data.
if(length(Group.T)>1 & length(Group.C)>1)
{
if(!is.null(s.C)) stop("Since you've specified raw group data, you should not specify any standard deviations.", call.=FALSE)
x.bar.T <- mean(Group.T)
x.bar.C <- mean(Group.C)
s.Con <- sqrt(var(Group.C))
d <- (x.bar.T - x.bar.C)/s.Con
if(Unbiased==TRUE) d <- d*gamma((n.C-1)/2)/(sqrt((n.C-1)/2)*gamma(((n.C-1)-1)/2))
return(d)
}
}
