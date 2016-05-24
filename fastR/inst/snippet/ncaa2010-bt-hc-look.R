# the "order effect" is the coefficient on "at.home"
coef(ncaa.model2)["at.home"] -> oe; oe
# expressed a multiplicative odds factor
exp(oe)
# prob home team wins if teams are "equal"
require(faraway); ilogit(oe)   
