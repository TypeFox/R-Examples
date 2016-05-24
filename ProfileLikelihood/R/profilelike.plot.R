profilelike.plot <-
function(theta=theta, profile.lik.norm=profile.lik.norm, round=2){

mle <- round(max(theta[profile.lik.norm==max(profile.lik.norm)]),round)
theta1.x.p <- theta[profile.lik.norm >=(1/8)]
theta2.x.p <- theta[profile.lik.norm >=(1/20)]
theta3.x.p <- theta[profile.lik.norm >=(1/32)]
theta.x.p.norm <- theta[profile.lik.norm >=0.146] #6.849
li.x.p.norm <- rep(0.146, length(theta.x.p.norm))
li8.x.p <- rep(1/8, length(theta1.x.p))
li20.x.p <- rep(1/20, length(theta2.x.p))
li32.x.p <- rep(1/32, length(theta3.x.p))

plot(theta, profile.lik.norm, type="l", lty=1, lwd=1, ylim=c(0,1), xlim=c(min(theta), max(theta)), ylab="", xlab=expression(theta))
lines(theta.x.p.norm, li.x.p.norm, lty=1, lwd=1, col="violet")
lines(theta1.x.p, li8.x.p, lty=1, lwd=1, col=4)
lines(theta2.x.p, li20.x.p, lty=1, lwd=1, col=2)
lines(theta3.x.p, li32.x.p, lty=1, lwd=1, col=3)
abline(v=0, lty=2, lwd=1.2, col="gray")
if(mle > 0){pos=0.15
} else{pos=0.8
}
text(theta[pos*length(theta)], 1, paste("Max at  ", round(max(theta[profile.lik.norm==max(profile.lik.norm)]),round)), cex=0.9, col=1)
text(theta[pos*length(theta)], 0.95, paste("1/6.8 LI (", round(min(theta.x.p.norm),round), ",", round(max(theta.x.p.norm),round), ")" ), cex=0.9, col="violet")
text(theta[pos*length(theta)], 0.90, paste("1/8 LI (", round(min(theta1.x.p),round), ",", round(max(theta1.x.p),round), ")" ), cex=0.9, col=4)
text(theta[pos*length(theta)], 0.85, paste("1/20 LI (", round(min(theta2.x.p),round), ",", round(max(theta2.x.p),round), ")" ), cex=0.9, col=2)
text(theta[pos*length(theta)], 0.80, paste("1/32 LI (", round(min(theta3.x.p),round), ",", round(max(theta3.x.p),round), ")" ), cex=0.9, col=3)

}

