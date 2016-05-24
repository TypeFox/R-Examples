pbino=function(nplus)
{

# POSTERIOR PROBABILITIES OF N FOR THE BINOMIAL CAPTURE MODEL
# UNDER UNIFORM PRIOR (N_MAX=400)

prob=c(rep(0,max(nplus,1)-1),1/(max(nplus,1):400+1))
prob/sum(prob)
}
