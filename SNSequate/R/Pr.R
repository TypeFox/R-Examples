### "Pr" provides the 3PL item response curve, or the probability that subject
### will give a correct answer to the item using the logistic model.

### Requires the value "theta" of ability, and the three item parameter values:
### "b" for difficulty, "a" for discrimination (by-default value of 1) and "c"
### for pseudo-guesing (by-default value of 0).


Pr<-function(theta,b,a=1,c=0) c+(1-c)/(1+exp(-a*(theta-b)))

