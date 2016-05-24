library(car)	# for OBrienKaiser data

## OBK.R	HE plots for repeated measures designs, OBrienKaiser data

# simplified analysis of OBrienKaiser data, collapsing over hour
OBK <- OBrienKaiser
OBK$pre <- rowMeans(OBK[,3:7])
OBK$post <- rowMeans(OBK[,8:12])
OBK$fup <- rowMeans(OBK[,13:17])
# remove separate hour scores
OBK <- OBK[,-(3:17)]
#contrasts(OBK$treatment) <- matrix(c(-2, 1, 1,  0, -1, 1), ncol=2)


# MANOVA model
mod.OBK <- lm(cbind(pre, post, fup) ~  treatment*gender,  data=OBK)


# for linear and quadratic effects of 'Time'
session <- ordered(c("pretest", "posttest", "followup"),
    levels=c("pretest", "posttest", "followup"))
# for profile contrasts
contrasts(session) <- matrix(c(-1,  1, 0,
		                      0, -1, 1), ncol=2)
#colnames(contrasts(session)) <- c("Post-Pre", "Fol-Post")

idata <- data.frame(session)

# Multivariate tests for repeated measures
aov.OBK <- Manova(mod.OBK, idata=idata, idesign=~session, type="III")
aov.OBK

# Univariate tests for repeated measures
summary(aov.OBK, multivariate=FALSE)


# HE plots for Between-S effects
heplot(mod.OBK, hypotheses=c("treatment1", "treatment2"),
	col=c("red", "black", "blue", "brown", "gray", "gray"),
	hyp.labels=c("(A,B)-Control", "A-B"),
	main="Between-S effects: Treat*Gender"
	)
pairs(mod.OBK, col=c("red", "black", "blue", "brown"))

heplot3d(mod.OBK, hypotheses=c("treatment1", "treatment2"),
  col=c("pink", "black", "blue", "green3", "gray40", "gray40"),
  hyp.labels=c("(A,B)-Control", "A-B"))

	
# HE plots for Within-S effects
	heplot(mod.OBK, idata=idata, idesign=~session, iterm="session",
		main="Within-S effects: Session: (Treat*Gender)")


### Plotting the Within-S effects by 'manually' transforming Y -> Y M

# Transform to profile contrasts for within-S effects
OBK$session.1 <- OBK$post - OBK$pre
OBK$session.2 <- OBK$fup - OBK$post

mod1.OBK <- lm(cbind(session.1, session.2) ~ treatment*gender,  data=OBK)
# HE plots for Within-S effects
heplot(mod1.OBK,
	main="Within-S effects: Session * (Treat*Gender)",
	remove.intercept=FALSE, type="III",
	xlab="Post-Pre", ylab="Fup-Post",
	term.labels=c("session", "treatment:session", "gender:session", 
				"treatment:gender:session"),
	col=c("red", "black", "blue", "brown"),
	xlim=c(-2,4), ylim=c(-2,3)
)
mark.H0()

# Main effect of session tests H0: Intercept=0
mod2.OBK <- lm(cbind(session.1, session.2) ~ 1,  data=OBK)
heplot(mod2.OBK,
	terms="(Intercept)", col=c("red", "blue"), type="3",
#	remove.intercept=FALSE,
	main="Within-S effects (profile contrasts): session",
	xlab="Post-Pre", ylab="Fup-Post",
	term.labels="session",
	xlim=c(-2,4), ylim=c(-2,3)
)
mark.H0()
