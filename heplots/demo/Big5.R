# Grice & Iwaski (2007) consider data on the sub-scales of the NEO-PIr,
# measuring the 'Big 5' personality traits (Neuroticism, Extraversion,
# Openness, Agreeableness and Conscientiousness) for college students
# classified into three groups (European Americans, Asian Americans
# and Asian Internationals).
# 
# References: 
# Grice & Iwasaki (2007). A truly multivariate approach to MANOVA.
# 	Applied Multivariate Research, 12(3), 199-226.
 
 
library(foreign)

Big5file <- "http://psychology.okstate.edu/faculty/jgrice/personalitylab/Iwasaki_Personality_Data.sav"
Big5 <-read.spss(Big5file)

# shorter group labels
Big5$GRP <- factor(Big5$GRP, labels=c("European", "AsianAmer","AsianIntl" ))
# use Helmert contrasts for group
contrasts(Big5$GRP) <- 
	matrix(c(2, -1, -1,
	         0, -1,  1), ncol=2)

str(Big5)

(Big5.mod <-lm(cbind(N,E,O,A,C) ~ GRP, data=Big5))
library(heplots)
Anova(Big5.mod)

colors <- c("red", "black", "blue", "brown")
hyp.labels <-c("Asian-European", "AsAmer-AsIntl")
heplot(Big5.mod, 
	xlab="Neuroticism", ylab="Extraversion", 
	hypotheses=c("GRP1", "GRP2"),
	col=colors)
	
# all pairwise plots
pairs(Big5.mod,
	hypotheses=c("GRP1", "GRP2"),
	col=colors)


# examine Extraversion, Openness, Anxiety in 3D 
colors <- c("pink", "gray", "blue", "brown")
heplot3d(Big5.mod, var=c(2:4),
	col=colors)

############## canonical discriminant HE plots #####################

library(candisc)
Big5.can <- candisc(Big5.mod)
Big5.can
heplot(Big5.can)

# 1-dim plot
Big5.can1 <- candisc(Big5.mod, data=Big5, ndim=1)
plot(Big5.can1)
