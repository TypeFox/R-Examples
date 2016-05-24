smd <- function(Group.1=NULL, Group.2=NULL, Mean.1=NULL, Mean.2=NULL, s.1=NULL, s.2=NULL, s=NULL, n.1=NULL, n.2=NULL, Unbiased=FALSE)
{
# Function to calculate the standardized mean difference.
if(length(Group.1)>=1 & length(Group.2)>=1 & (!is.null(Mean.1)|!is.null(Mean.2)|!is.null(s.1)|!is.null(s.2)|!is.null(s))) stop("Since you've specified raw data, you do not need to specify summary values.", call.=FALSE)
if(length(Group.1)==1 | length(Group.2)==1) stop("You only have 1 individual in one or both groups, you cannot calculate a variance in such a situation.", call.=FALSE)
if(length(Mean.1)==1 & length(Mean.2)==1 & !is.null(s) & (is.null(n.1) | is.null(n.2)) & Unbiased==TRUE) stop("You need to specify the group sample sizes in order to obtain the unbiased estiamte of the standardized mean difference.", call.=FALSE)

# Here determine 'smd' on summary data.
if(length(Mean.1)==1 & length(Mean.2)==1)
{
 if(!is.null(Group.1) | !is.null(Group.2)) stop("Group.1 and Group.2 should be NULL, since you've specified the group means directly. Alternatively, you can specify the groups directly and make Mean.1 and Mean.2 NULL.", call.=FALSE)
     if(is.null(s.1) & is.null(s.2) & length(s)==1)
        {
        d <- (Mean.1 - Mean.2)/s
        if(Unbiased==TRUE) d <- d*gamma((n.1+n.2-2)/2)/(sqrt((n.1+n.2-2)/2)*gamma(((n.1+n.2-2)-1)/2))
        return(d)   
        }   
        
if(length(s.1)==1 & length(s.2)==1 & is.null(s))
{
if(is.null(n.1) | is.null(n.2) | length(n.1)==0 | length(n.2)==0) stop("You did not specify the per group sample sizes (i.e., 'n.1' and/or 'n.2')", .call=FALSE)
# Calculate the pooled variance.
sd <- sqrt((s.1^2*(n.1-1) + s.2^2*(n.2-1))/(n.1 + n.2 - 2))
d <- (Mean.1 - Mean.2)/sd
if(Unbiased==TRUE) d <- d*gamma((n.1+n.2-2)/2)/(sqrt((n.1+n.2-2)/2)*gamma(((n.1+n.2-2)-1)/2))
return(d)
}
}

# Here determine 'smd' for raw data.
if(length(Group.1)>1 & length(Group.2)>1)
{
if(!is.null(s.1) | !is.null(s.2) | !is.null(s)) stop("Since you've specified raw group data, you should not specify any standard deviations.", call.=FALSE)
x.bar.1 <- mean(Group.1)
x.bar.2 <- mean(Group.2)
SS1 <- sum((Group.1-x.bar.1)^2)
SS2 <- sum((Group.2-x.bar.2)^2)
n.1 <- length(Group.1)
n.2 <- length(Group.2)
s <- sqrt((SS1 + SS2)/(n.1+n.2-2))
d <- (x.bar.1 - x.bar.2)/s
if(Unbiased==TRUE) d <- d*gamma((n.1+n.2-2)/2)/(sqrt((n.1+n.2-2)/2)*gamma(((n.1+n.2-2)-1)/2))
return(d)
}
}
