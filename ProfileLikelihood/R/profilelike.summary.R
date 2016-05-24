profilelike.summary <-
function(k, theta=theta, profile.lik.norm=profile.lik.norm, round=2){
mle <- round(max(theta[profile.lik.norm==max(profile.lik.norm)]),round)
theta1.x.p <- theta[profile.lik.norm >=(1/k)]
li.k.p <- rep(k, length(theta1.x.p))
theta.x.p.norm <- theta[profile.lik.norm >=0.146] #6.849
li.x.p.norm <- rep(0.146, length(theta.x.p.norm))

LI.norm <- c(round(min(theta.x.p.norm),round), round(max(theta.x.p.norm),round) )
LI.k <- c(round(min(theta1.x.p),round), round(max(theta1.x.p),round) )

return(list(k=k, mle=mle, LI.k = LI.k,  LI.norm = LI.norm))
}

