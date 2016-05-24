pdarroch=function(n1,n2,m2)
{

# POSTERIOR PROBABILITIES OF N FOR THE DARROCH MODEL
# UNDER UNIFORM PRIOR (N_Max=400)

prob=c(rep(0,max(n1+n2-m2,1)-1),choose(n1,m2)*choose(max((n1+n2-m2),1):400-n1,n2-m2)/choose(max((n1+n2-m2),1):400,n2))
prob/sum(prob)
}
