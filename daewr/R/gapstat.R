gapstat<-function(beta,pse) {
# computes the standardized gap score
p<-length(beta)
psehe<-pse
# gets positive coefficients
sel<-beta >= 0
betap<-beta[sel]
# sorts positive elements
betap<-sort(betap)
# gets Beta_s
betas<-betap[1]
#gets negative coefficients
sel<-beta < 0
betan<-beta[sel]
nn<-length(betan)
# sorts negative coefficients
betan<-sort(betan)
#gets Beta_L
betal<-betan[nn]
# gets Z_L and Z_S
zl<-qnorm((nn-.375)/(p+.25))
zs<-qnorm((nn+1-.375)/(p+.25))
# calculates gap statistic
gap<-((betas-betal)/psehe)/(zs-zl)
return(gap)
                      }
