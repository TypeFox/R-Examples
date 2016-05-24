"cleanstrata" <-
function(strata)
{
a=sort(unique(strata)) 
b=1:length(a) 
names(b)=a 
as.vector(b[as.character(strata)]) 
}

