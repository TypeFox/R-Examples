## housing.R    Visualize models from example(housing, package="MASS")

# These examples fit a variety of models to the data(housing), giving a 4-way
# frequency table of 1681 individuals from the Copenhagen Housing Conditions
# Survey, classified by their Type of rental dwelling, perceived Influence on
# management of the property, and degree of Contact with other residents. The
# response variable here is Satisfaction of householders with their present
# housing circumstances.


library(vcdExtra)
data(housing, package="MASS")


oldop <-options(contrasts = c("contr.treatment", "contr.poly"))

##########################
#  Poisson models for Freq, equivalent to loglinear models
##########################
# Baseline model, with Satisfaction as a response
house.glm0 <- glm(Freq ~ Infl*Type*Cont + Sat, family = poisson,
                  data = housing)
modFit(house.glm0)

# labeling_args for mosaic()
largs <- list(set_varnames = c(Infl="Influence on management", 
			Cont="Contact among residents", Type="Type of dwelling", Sat="Satisfaction"),
	abbreviate=c(Type=3))

mosaic(house.glm0, labeling_args=largs, main='Baseline model: [ITC][Sat]')
# reorder variables in the mosaic, putting Sat last
mosaic(house.glm0, ~ Type+Infl+Cont+Sat, labeling_args=largs, main='Baseline model: [ITC][Sat]')

# what terms need to be added?
addterm(house.glm0, ~. + Sat:(Infl+Type+Cont), test = "Chisq")

# add all two way terms with Satisfaction
house.glm1 <- update(house.glm0, . ~ . + Sat*(Infl+Type+Cont))

# did it get better?
anova(house.glm0, house.glm1, test="Chisq")
mosaic(house.glm1, labeling_args=largs, main='Model [IS][TS][CS]', gp=shading_Friendly)

# Same model, fit by terative proportional scaling
(house.loglm <- loglm(Freq ~ Infl*Type*Cont + Sat*(Infl+Type+Cont), data = housing))

# Can we drop any terms?
dropterm(house.glm1, test = "Chisq")

# Need to add any terms?
addterm(house.glm1, ~. + Sat:(Infl+Type+Cont)^2, test  =  "Chisq")


##########################
# Effect plots, for glm1 model
##########################
library(effects)
house.eff <-allEffects(house.glm1)

# show the interactions of Infl, Cont and Type with Sat
plot(house.eff, 'Infl:Sat', x.var='Sat', xlab="Satisfaction")
# same plot in one panel, no std errors shown
plot(house.eff, 'Infl:Sat', x.var='Sat', xlab="Satisfaction", multiline=TRUE)

plot(house.eff, 'Cont:Sat', x.var='Sat', xlab="Satisfaction")

plot(house.eff, 'Type:Sat', x.var='Sat', xlab="Satisfaction")


##########################
# multinomial model
##########################
library(nnet)
# multinomial model, similar in spirit to house.glm1
(house.mult<- multinom(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing))
# Do we need a more complex model?
house.mult2 <- multinom(Sat ~ Infl*Type*Cont, weights = Freq,
                        data = housing)
anova(house.mult, house.mult2)


# effect plots for multinomial model
house.effm <- allEffects(house.mult)

plot(house.effm, 'Infl', xlab='Influence on management', style="stacked",
	main="Multinomial: Infl effect plot")

plot(house.effm, 'Cont', xlab='Contact among residents', style="stacked",
	main="Multinomial: Cont effect plot")

plot(house.effm, 'Type', xlab='Type of dwelling', style="stacked",
	main="Multinomial: Type effect plot")


##########################
# proportional odds model
##########################

(house.plr <- polr(Sat ~ Infl + Type + Cont,
                   data = housing, weights = Freq))

# Test proportional odds assumption by likelihood ratio test
# NB: multinom() objects do not have a df.residual component, so we have
#     to use the difference in edf to get df for the test
pchisq(deviance(house.plr) - deviance(house.mult), 
		df = house.mult$edf -house.plr$edf, lower.tail = FALSE)

# try more complex models
house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova

house.effp <- allEffects(house.plr)

plot(house.effp, 'Infl', xlab='Influence on management', style="stacked",
	main="Proportional odds: Infl effect plot")
plot(house.effp, 'Cont', xlab='Contact among residents', style="stacked",
	main="Proportional odds: Cont effect plot")
plot(house.effp, 'Type', xlab='Type of dwelling', style="stacked",
	main="Proportional odds: Type effect plot")

options(oldop)
