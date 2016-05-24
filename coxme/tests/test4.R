library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# A further test of mixed sparse/ non-sparse computation
#  Make sure all the terms are added in and taken out properly,
#  when people enter and leave the risk set. 
#
# The choice of "sparse" below forces group 1 to be non-sparse and
#   group 2-4 to be sparse.  It will also force the coefficients to be
#   in 2,3,4,1 order

tdata1 <- data.frame(time  =c(5,4,1,1,2,2,2,2,3, 1:8), 
                     status=c(0,1,1,0,1,1,1,0,0, rep(1:0,4)),
                     x1    =c(0,1,2,0,1,1,0,1,0, rep(4:1,2)),
                     wt    =c(1,2,1,2,3,4,3,2,1, rep(1:2,4)),
                     x2    =c(1,3,5,2,3,6,4,3,1, rep(1:2,4)),
                     grp   =c(1,1,2,2,1,1,2,2,1, rep(3,4), rep(4,4)))
temp1 <-  5*tdata1$time -3  # each subject spends 4 days not at risk
tdata2 <- data.frame(t1 = c(rep(0, length(temp1)), temp1+4),
		     t2 = c(temp1, 10*tdata1$time),
                     status=c(tdata1$status, rev(tdata1$status)),
		     x1 = rep(tdata1$x1, 2),
		     wt = rep(tdata1$wt, 2),
		     x2 = rep(tdata1$x2, 2),
		     grp= rep(tdata1$grp, 2),
		     id = rep(1:17, 2))

# What's it look like?
#plot(c(0,80), c(0,17), type='n')
#segments(tdata2$t1, tdata2$id, tdata2$t2, tdata2$id)
#points(tdata2$t2[tdata2$status==1], tdata2$id[tdata2$status==1], pch=1)

theta <- .83
fit1 <- coxme(Surv(t1, t2, status) ~ x1 + x2 +(1|grp), data=tdata2,
	      vfixed=theta, ties='breslow', 
              sparse=c(2, .25))

# This fit will complain when it tries to invert the information matrix--
#  ignore the complaint
fit2 <- coxph(Surv(t1, t2, status) ~ I(grp==2) + I(grp==3) + I(grp==4) +
	      I(grp==1) + x1 + x2, data=tdata2, method='breslow',
	      init=c(unlist(ranef(fit1)), fixef(fit1)),
	      iter=0)

lp <- c(cbind(tdata2$x1, tdata2$x2) %*% fixef(fit1) + 
	unlist(ranef(fit1))[c(4,1,2,3)[tdata2$grp]])
aeq(lp, fit1$linear)
fit3 <- coxph(Surv(t1, t2, status) ~ offset(lp), data=tdata2, method='breslow')

aeq(fit1$loglik[3] + fit1$penalty, fit3$loglik)
aeq(fit3$loglik, fit2$loglik[2])


#
# And now, do it mostly by hand as the definitive check
#
lp2 <- fit1$linear - sum(fit1$mean * fixef(fit1))
temp.risk <- exp(lp2)
dtimes <- sort(unique(tdata2$t2[tdata2$status==1]))
nd <- length(dtimes)
denom <- double(nd)
xbar  <- matrix(0., nrow=nd, ncol=6)
xtemp <- cbind(1*(tdata2$grp==2), 1*(tdata2$grp==3), 1*(tdata2$grp==4),
	       1*(tdata2$grp==1), tdata2$x1, tdata2$x2)

for (i in 1:nd) {
    who <- (tdata2$t1 <dtimes[i] & tdata2$t2 >= dtimes[i])
    denom[i] <- sum(temp.risk[who])
    xbar[i,] <- colSums(temp.risk[who] * xtemp[who,])/denom[i]
    }
temp.log <- sum(lp2[tdata2$status==1] -
		log(denom[match(tdata2$t2[tdata2$status==1], dtimes)]))
aeq(fit1$loglik[3] + fit1$pen, temp.log)
aeq(fit1$penalty, sum(unlist(ranef(fit1))^2)/(2*theta))
# The linear predictors will differ by a constant, since fit2 will have
#   subtracted means from the factor terms.
aeq(diff(range(fit1$linear-fit2$linear)),0)

# Score statistic
temp.score <- xtemp[tdata2$status==1,] - 
	      xbar[match(tdata2$t2[tdata2$status==1], dtimes),]
aeq(colSums(temp.score) - c(unlist(ranef(fit1)),0,0)/theta, fit1$u)


dt <- coxph.detail(fit2, riskmat=T)
aeq(dt$x, xtemp[as.numeric(dimnames(dt$x)[[1]]),])
nevent <- table(tdata2$t2, tdata2$status)[,2]
nevent <- nevent[nevent>0]
aeq(dt$nevent, nevent)
aeq(dt$means, xbar)
aeq(colSums(temp.score), colSums(dt$score))

# Check out the information matrix.
#  The coxme model has extra added to the diagonal, and due to
#  the sparseness the off-diagonal elements for vars 1-3 are zero.
temp1 <- as.matrix(fit1$hmat) %*% diag(sqrt(diag(fit1$hmat)))
temp.imat <- temp1 %*% t(temp1)
dt.imat <- apply(dt$imat,1:2, sum)
aeq(diag(temp.imat), diag(dt.imat) + c(rep(1/theta,4),0,0))
aeq(temp.imat[4:6, 3], dt.imat[4:6,3])
aeq(temp.imat[5:6, 4:6], dt.imat[5:6,4:6])

