"disjunctive" <-
function(strata)
{ss=cleanstrata(strata)
m=matrix(0,length(strata),length(unique(strata)))
for(i in 1:length(ss)) m[i,ss[i]]=1
m
}

