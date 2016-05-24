#
# Gerber, Alan S. and Donald P. Green. 2000. "The Effects of Canvassing, Telephone Calls, and
# Direct Mail on Voter Turnout: A Field Experiment." American Political Science Review 94: 653-663.
#
# Imai, Kosuke. 2005. "Do Get-Out-The-Vote Calls Reduce Turnout? The Importance of 
# Statistical Methods for Field Experiments".  American Political Science Review 99: 283-300.
#

set.seed(10391)
data(GerberGreenImai)

#replication of Imai's propensity score matching model
pscore.glm<-glm(PHN.C1 ~ PERSONS + VOTE96.1 + NEW + MAJORPTY + AGE +
                WARD + PERSONS:VOTE96.1 + PERSONS:NEW + AGE2, 
                family=binomial(logit), data=GerberGreenImai)

D<-GerberGreenImai$PHN.C1 #treatment phone calls
Y<-GerberGreenImai$VOTED98 #outcome, turnout

cat("\nTHIS MODEL FAILS TO BALANCE AGE\n")
X  <- fitted(pscore.glm)

#propensity score matching estimator
r1  <- Match(Y=Y, Tr=D, X=X, M=3)
summary(r1)

#check for balance before and after matching
mb1  <- MatchBalance(PHN.C1 ~ AGE + AGE2 + PERSONS + VOTE96.1 + NEW + MAJORPTY +
                     WARD + I(PERSONS*VOTE96.1) + I(PERSONS*NEW), match.out=r1,
                     data=GerberGreenImai)

