F.cjs.gof <- function( cjsobj, resid.type="pearson", rule.of.thumb = 2, HL.breaks = "deciles" ){


if( !inherits( cjsobj, "cjs" ) ){
	stop("Wrong class: An object of class 'cjs' must be supplied. I.e., output from F.cjs.estim.")  
}


e <- cjsobj$fitted
o <- as.numeric( cjsobj$histories >= 1 )
p <- length(cjsobj$parameters)
active <- !is.na(e)

if( cjsobj$resid.type == "pearson" ){
	r <- cjsobj$residuals
} else {
	r <- residuals( cjsobj, type="pearson" )
}

#	McDonald and Regher overall test
gof.chi <- sum( r^2, na.rm = TRUE )
gof.df <- sum( !is.na(r) )   # conservative.
gof.p <- ifelse( gof.df > 0, 1-pchisq(gof.chi,gof.df), NA )


#	Osius and Rojek normal approximation to overall Chi-square.  
#	See p. 153 of Hosmer and Lemeshow (2000)
#	The real Osius-Rojek procedure
#	regresses the pseudo-response c on x, where x are the covariate patterns. Here, we do not have 
#	x's, so we fit an intercept only model which should be conservative in the sense that it 
#	will yeild a RSS that is larger than the RSS that would have been obtained had we found some
#	x's to use. 
e.factor <- c(e)[ !is.na(c(e)) ]
e.ord    <- order( e.factor )
e.factor <- e.factor[ e.ord ]
e.factor <- as.factor( round(e.factor, digits=6 ) )
e.unique <- as.numeric( as.character( unique( e.factor )))
o.unique <- c(o)[ !is.na(c(e)) ]
o.unique <- o.unique[ e.ord ]
o.unique <- tapply( o.unique, e.factor, FUN=sum )
m.j <- table(e.factor)
v.j <- m.j * e.unique * (1-e.unique)
c.j <- (o.unique - m.j*e.unique)^2 / v.j
u.j <- rep(1, length(c.j))
or.table <- rbind( o.unique, m.j*e.unique, m.j, v.j, c.j, u.j )
dimnames( or.table )[[1]] <- c("obs", "exp", "n", "var", "contrib", "use?")

# Now compute the correction/approximation to normal distribution.
or.chi <- sum( c.j )   # if all unique predicted values, this should equal gof from above.
J <- length(e.unique)
or.df <- J - p - 1
c.j <- (1-2*e.unique) / v.j
#	wfit <- lm( c.j ~ 1, weights=v.j ); rss <- sum( wfit$residuals^2 )   # we could do this,
or.rss <- sum( v.j*(c.j - mean(c.j))^2 )    # or this since we are fitting only the intercept.
or.correction <- 2*(J - sum(1/m.j))
or.z <- (or.chi - or.df) / sqrt( or.correction + or.rss ) 
or.p <- 2 * pnorm(abs(or.z), lower.tail=FALSE)



#	T4: sum over individuals
o.j <- apply( o*active, 2, sum, na.rm = TRUE )
e.j <- apply( e, 2, sum, na.rm = TRUE )
u.j <- as.numeric( e.j >= rule.of.thumb )
c.j <- ((o.j-e.j)*(o.j-e.j))/e.j
c.j[ u.j == 0 ] <- NA
t4.chi <- sum( c.j, na.rm=TRUE )
t4.df <- sum( u.j, na.rm=TRUE ) - 1
t4.table <- rbind( o.j, e.j, c.j, u.j )
dimnames( t4.table )[[1]] <- c("obs", "exp", "contrib", "use?")
t4.p <- ifelse( t4.df > 0, 1-pchisq(t4.chi,t4.df), NA )


#	T5: sum over occasions
o.j <- apply( o*active, 1, sum, na.rm = TRUE )
e.j <- apply( e, 1, sum, na.rm = TRUE )
u.j <- as.numeric( e.j >= rule.of.thumb )
c.j <- ((o.j-e.j)*(o.j-e.j))/e.j
c.j[ u.j == 0 ] <- NA
t5.chi <- sum( c.j, na.rm=TRUE )
t5.df <- sum( u.j, na.rm=TRUE ) - 1
t5.table <- rbind( o.j, e.j, c.j, u.j )
dimnames( t5.table )[[1]] <- c("obs", "exp", "contrib", "use?")
t5.p <- ifelse( t5.df > 0, 1-pchisq(t5.chi,t5.df), NA )


#	Hosmer-Lemeshow statistic
if( is.character( HL.breaks ) ){
	if( HL.breaks[1] == "deciles" ){
		HL.breaks <- quantile(e, probs=seq(.1, .9, by=.1), na.rm=TRUE)
	} else {
		stop("HL.breaks must be numeric or 'deciles'.")
	}
} else {
	if( min(HL.breaks) < 0 | max(HL.breaks) > 1 ) stop("HL.breaks must be between 0 and 1.")
}

HL.bins <- sort(unique( c(0,HL.breaks,1) ))

group <- cut( e, br=HL.bins )
o.j <- tapply( c(o), group, FUN=sum )
e.j <- tapply( c(e), group, FUN=sum )
n.j <- tapply( c(e), group, FUN=length ) # should be ~n/10 for all 
c.j <- (o.j - e.j)^2 / (e.j*(1-e.j/n.j))
u.j <- as.numeric( e.j >= rule.of.thumb )
c.j[ u.j == 0 ] <- NA
HL.chi <- sum( c.j, na.rm=TRUE )
HL.df <- sum( u.j, na.rm=TRUE ) - 2
HL.table <- rbind( o.j, e.j, c.j, u.j )
dimnames( HL.table )[[1]] <- c("obs", "exp", "contrib", "use?")
HL.p <- ifelse( HL.df > 0, 1-pchisq(HL.chi,HL.df), NA )

#	Area under ROC = Mann-Whitney U of predicted p's where samples are defined by 0's and 1's
e.ones <- e[ o == 1 ]
e.ones <- e.ones[ !is.na(e.ones) ]
e.zero <- e[ o == 0 ]
e.zero <- e.zero[ !is.na(e.zero) ]
mw.mat <- outer( e.zero, e.ones, FUN=function(p.0, p.1){ ifelse( abs(p.0 - p.1) <= 1e-6 , 0.5, as.numeric(p.0 < p.1))} )
roc <- sum( mw.mat ) / (length(e.ones) * length(e.zero))

#	Perhaps some version of the Stukel test (Hosmer-Lemeshow, 2000, p. 155) could be worked 
#	out for this setting, but without explicit x's, it would be hard.  I don't see it. But, 
#	perhaps you could work out the geometric form that cell predicted values should take on 
#	if the logistic link function is correct. If so, perhaps you could do a Stukel-like test 
#	by regressing the predicted values on a surrogate for this form. 


ans <- c( cjsobj, list(
	gof.chi = gof.chi, gof.df=gof.df, gof.pvalue=gof.p,
	or.table = or.table, or.chi = or.chi, or.z = or.z, or.df = or.df, or.pvalue = or.p, or.correction=or.correction, or.rss=or.rss,
	t4.table=t4.table, t4.chi=t4.chi, t4.df=t4.df, t4.pvalue=t4.p, 
	t5.table=t5.table, t5.chi=t5.chi, t5.df=t5.df, t5.pvalue=t5.p,
	HL.table=HL.table, HL.chi=HL.chi, HL.df=HL.df, HL.pvalue=HL.p,
	roc=roc))

class(ans) <- c("cjsgof","cjs","cr")

ans

}

