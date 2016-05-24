# Guo & Krishnamoorthy function
# it provides the solution of the equation
# $P(F_t(qu;df,x)=p)$
# where $F_t$ is the cdf (calculated in $qu$) of a non-central Student r.v. with $df$ degrees of freedom
# and non-centrality parameter=$x$.
# In R code, gkf provides the solution of $pt(qu,df,x)=p$
gkf<-function(p,q,df,eps=0.00001)
{
# eps=maximum error allowed
qq<-q-qnorm(p)
pp<-pt(q,df,qq)
while(abs(pp-p)>eps)
{
qq<-qq-sign(p-pp)*abs(qnorm(p)-qnorm(pp))
pp<-pt(q,df,qq)
}
qq
}
