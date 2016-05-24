##load a model object
data(mesa.model)

##create vector of initial values
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(c( rep(2, dim$nparam.cov-1), 0),
                c( rep(c(1,-3), dim$m+1), -3, 0))
rownames(x.init) <- loglikeSTnames(mesa.model, all=FALSE)

\dontrun{
  ##estimate parameters
  est.mesa.model <- estimate(mesa.model, x.init, hessian.all=TRUE)
}

##time consuming estimation, load pre-computed results instead
data(est.mesa.model)

#estimation results
print(est.mesa.model)

##compare the estimated parameters for the two starting points
est.mesa.model$summary$par.all
##and values of the likelihood (and convergence info)
est.mesa.model$summary$status

##extract the estimated parameters and approximate uncertainties
x <- coef(est.mesa.model)

##compare estimated parameters
##plot the estimated parameters with uncertainties
par(mfrow=c(1,1),mar=c(13.5,2.5,.5,.5))
with(x, plot(par, ylim=range(c(par-1.96*sd, par+1.96*sd)),
                xlab="", xaxt="n"))
with(x, points(par - 1.96*sd, pch=3))
with(x, points(par + 1.96*sd, pch=3))

abline(h=0, col="grey")
##add axis labels
axis(1, 1:length(x$par), rownames(x), las=2)

\dontrun{
  ##example using a few fixed parameters
  x.fixed <- coef(est.mesa.model)$par
  x.fixed[c(1,2,5:9)] <- NA
  est.fix <- estimate(mesa.model, x.init, x.fixed, type="p")
}
