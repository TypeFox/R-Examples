pcapture=function(T,nplus,nc)
{

# POSTERIOR PROBABILITIES OF N FOR THE T-STAGE CAPTURE-RECAPTURE MODEL
# UNDER UNIFORM PRIOR (N_Max=400)

lprob=lchoose(max(nplus,1):400,nplus)+lgamma(T*max(nplus,1):400-nc+1)-lgamma(T*max(nplus,1):400+2)
prob=c(rep(0,max(nplus,1)-1),exp(lprob-max(lprob)))
prob/sum(prob)
}
