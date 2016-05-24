## This is just from ?anova.rq (extended)
library(quantreg)
data(barro)
fit0 <- rq(y.net ~  lgdp2 + fse2 + gedy2 , data = barro)
fit1 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro)

a01 <- anova(fit1,fit0)
a01

fit2 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,
	   tau = 0.75)
fit3 <- rq(y.net ~  lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro,
	   tau = 0.25)

a123  <- anova(fit1,fit2, fit3)
a.123 <- anova(fit1,fit2, fit3, joint=FALSE)
a.123

AE <- function(x,y) all.equal(x, y, tol = 1e-5)
##                                  ---------- {giving a bit more digits below}
stopifnot(
	  AE(100 * unname(coef(fit0)),
	     c(-0.74679759, 0.46539963, 0.15902838, -36.619915))
	  ,
	  AE(unlist(a01$table),
	     c(ndf=2, ddf=155, Tn = 18.878717, pvalue= 4.6e-08))
	  ,
	  AE(100* unname(coef( fit2 )),
	     c(13.103018, -1.4885239, -0.026452369,
	       0.3999839, 14.526663, -13.504643))
	  ,
	  AE(100* unname(coef( fit3 )),
	     c(6.0860719, -0.88350554, 0.24596781,
	       -14.962498, 15.592489, -15.861804))
	  ,
	  AE(unlist(a123$table),
	     c(ndf = 10, ddf = 473, Tn = 1.80385526, pvalue=0.0575117558))
	  ,
	  AE(a.123$table[,"Tn"],
	     c(1.0655561, 2.6398508, 0.78623238, 0.04467014, 0.065344348))
	  )
