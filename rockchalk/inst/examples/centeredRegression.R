## Paul Johnson
## pauljohn@ku.edu 2012-03-09
##
## This is an R program that uses some functions in the newly revised
## rockchalk package to demonstrate my point that "centered variables"
## don't really make a difference in regressions with interaction.
##
## Centering does not help with multicollinearity, but I mean to say
## more than that. It does not help with regression interpretation,
## if one correctly understands what the parameter estimates and the
## predicted values mean in a regression with interaction.

## Here the idea is the following. The centered b's and se's "seem" different,
## but they are actually calculated on the EXACT SAME fitted plane. The
## difference is that the notcentered model has the y axis positioned at
## x1=0,x2=0, while in the centered model it is instead at
## x1=meanx1, x2=meanx2. The predicted values for any x1, x2 combination
## are EXACTLY the same with either model, as are the estimates of
## uncertainty (in the form of confidence intervals or standard errors).

## Thus it should be possible to take the estimates from the
## notcentered regression and calculate the slopes at x1=meanx1,
## x2=meanx2, and re-produce them.  AND I CAN!! This demonstrates
## that claim at the end.


library(rockchalk)
set.seed(222233)
dat3 <- genCorrelatedData(rho = .31, stde = 80, beta=c(0.1, 0.2, 0.3, -0.2))

## We now explore a regression model like this:

## y = b0 + b1 * x1 + b2 *x2 + b3 *x1*x2 + e

## The mean-centered model replaces x1 and x2 by the mean-centered
## versions  of those variables, (x1 - mean(x1)) and (x2 - mean(x2)).

## The usual problem that causes researchers to turn to "mean centering"
## is the following.  The linear model seems "good", but the interactive
## model seems "bad" without centering.  Here's an example.

## First, fit the model without the interaction term.
## y = b0 + b1 * x1 + b2 *x2 + e

m1 <- lm(y ~ x1 + x2, data= dat3)
summary(m1)
dev.new()
plotPlane(m1, plotx1="x1", plotx2="x2")


## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 424.1485    45.8004   9.261 5.29e-15 ***
## x1           -8.9042     0.8499 -10.477  < 2e-16 ***
## x2           -9.2701     0.8039 -11.531  < 2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 77.77 on 97 degrees of freedom
## Multiple R-squared: 0.8011,	Adjusted R-squared: 0.797
## F-statistic: 195.4 on 2 and 97 DF,  p-value: < 2.2e-16
##
## Yeah, its a "good model". Everything is "significant".
##
##
## Now the problem.
## Add an interaction. Watch, it ruins everything!
##

m2 <- lm(y ~ x1 * x2, data = dat3)
summary(m2)
dev.new()
plotPlane(m2, plotx1 = "x1", plotx2 = "x2", theta = -10)


## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)
## (Intercept) -74.45719  179.94096  -0.414  0.67995
## x1            1.57257    3.75575   0.419  0.67636
## x2            0.64618    3.55471   0.182  0.85614
## x1:x2        -0.20526    0.07181  -2.859  0.00522 **
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 75.05 on 96 degrees of freedom
## Multiple R-squared: 0.8167,	Adjusted R-squared: 0.811
## F-statistic: 142.6 on 3 and 96 DF,  p-value: < 2.2e-16


## Booh, the model's "no good". x1 and x2 "don't matter"
## any more.

## We'd better mean-center x1 and x2.

m2mc <- meanCenter(m2)
summary(m2mc)

## You can verify that manually, if you don't trust my meanCenter
## function. Results are the same, of course. First, mean center the
## variables with R's scale function:

dat3$x1c <- as.numeric(scale(dat3$x1, center = TRUE, scale = FALSE))
dat3$x2c <-  as.numeric(scale(dat3$x2, center = TRUE, scale = FALSE))
## The as.numeric is required because scale returns a matrix
## with one column, not a vector. How annoying!

m2manualmc <- lm(y ~ x1c*x2c, dat3)
summary(m2manualmc)
m2mcmanualpred <- fitted(m2manualmc)

## Confirmed that matches the output from meanCenter.  I suppose
## you are sorry you did not trust me now.


## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)
## (Intercept) -456.69849    8.01647 -56.970  < 2e-16 ***
## x1            -8.38383    0.84005  -9.980  < 2e-16 ***
## x2            -9.47969    0.77920 -12.166  < 2e-16 ***
## x1:x2         -0.20526    0.07181  -2.859  0.00522 **
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 75.05 on 96 degrees of freedom
## Multiple R-squared: 0.8167,	Adjusted R-squared: 0.811
## F-statistic: 142.6 on 3 and 96 DF,  p-value: < 2.2e-16


## Yeah! Mean centering saved the day! Woo Hoo! All your
## t values are big and your p values are small.

## Unfortunately, the difference in the output is just a mirage. The
## two fits are describing slope estimates at a DIFFERENT POINT on the
## same curving plane.

## To see that, note that the fitted models (centered or not centered) offer
## EXACTLY THE SAME predicted values! I don't mean "similar"
## I mean exactly the same. Look like this:

par(mfcol = c(1,2))

plotPlane(m2, "x1", "x2", plotPoints = FALSE, theta = -25,
          main = "Not Mean Centered", ticktype = "detailed")
plotPlane(m2mc, "x1", "x2", plotPoints = FALSE, theta = -25,
          main = "Mean Centered", ticktype = "detailed")

par(mfcol = c(1,1))

## Maybe you are not into 3d illustrations.  Too bad. I can show the
## same in two dimensions.  Let's create a scatterplot displaying the
## predicted values from the ordinary and the mean-centered models

plot(fitted(m2), fitted(m2mc), xlab = "predicted from uncentered x",
     ylab = "predicted from centered x",
     main = "(Not)Centered Predictions Identical")

##
##
## So, how can it be the two models are EXACTLY the same and yet
## the mean centered one "looks better"?
##
## Here's the simple answer. This is a nonlinear model, the slope
## changes from point to point.  If we pick a point where it is
## steep and make a "snapshot" of the slope, we conclude "the
## slope is huge!"  If we pick a point where the slope is small
## we conclude "there's no effect of x".  Its really just that
## simple. Mean centering amounts to deciding where to measure
## the slope, and sometimes we are lucky enough to accidentally
## pick a steep spot by mean centering.
##
## Now the proponents of the mean-centering approach respond to
## me as follows. "You ignore the fact that the standard errors
## are smaller in the mean centered model. So there is a difference."

## In response, I say this. Just as the slope varies from point to
## point, the standard error also changes from point to point.  If you
## choose a point in the data, say where x1 = 48 and x2 = 47.2, and then
## find the equivalent mean-centered point, which is x1c = -1.24 and
## x2c = -1.305, you will find the slopes and the standard errors are
## exactly the same.  The slope at x1,x2 in the not centered model is
## exactly the same as the slope in the centered model at x1c,x2c.
## More importantly, from either model, we can calculate the same
## values of not only the slope, but also the standard error, for any
## particular point in the data.

## What is the easiest way to see that the slopes and standard errors
## are the same?  It seems obvious to me that, since the predicted
## value plane is identical in the two models, then the slopes have
## to be the same as well. But perhaps that relys on intuition that
## is unique to me.

## Take the the noncentered model and calculate the slope and standard
## error at the mean of x1,x2.  The slope at that point in the x1
## dimension, in terms of the noncentered fit, is

## b1 + b3*x2
## Which we estimate from the noncentered model as:
coef(m2)["x1"] + coef(m2)["x1:x2"] * mean(dat3$x2, na.rm = TRUE)
##       x1
## -8.383827
##
## And the slope in the x2 dimension is
## b2 + b3*x1

coef(m2)["x2"] + coef(m2)["x1:x2"] * mean(dat3$x1, na.rm = TRUE)
##      x2
## -9.479689
##
## Please note, those estimates of the slopes EXACTLY match the
## coefficient estimates reported by the centered model. That is to
## say, if you restrict your attention to a particular value of
## (x1,x2), the centered and noncentered models produce EXATCLY the same
## slope estimates.

## And the standard errors are the same as well.
## Reproduce the standard errors in centered model from noncentered model

V <- vcov(m2)

sqrt(V["x1","x1"] + mean(dat3$x2)^2 * V["x1:x2","x1:x2"] +
     2 * mean(dat3$x2) * V["x1","x1:x2"])

## [1] 0.8400474

## That's the SAME number reported in the Std.Err. column for the centered
## model.

sqrt(V["x2","x2"] + mean(dat3$x1)^2 * V["x1:x2","x1:x2"] +
     2 * mean(dat3$x1) * V["x2","x1:x2"])

##[1]  0.7791997
##
## Bingo, Fits exactly. The estimates of the centered model are reproduced
## exactly from the notcentered model once the correct co-ordinate
## translation is put in place. The EXACT same t ratios, etc.

## Centering has NO EFFECT whatsoever on multicollinearity. It does
## not affect predicted values, it does not affect our confidence
## in estimates of slopes or predicted values.

## In short, if you understand what a multiple regression with
## interaction "really does," it is impossible to reach any conclusion
## except the following: mean centering makes no difference at all in
## the estimation or interpretation of regression models with
## interaction effects.

## Mean centering only aids the interpretation if one is too lazy to
## understand the curvature of the estimated plane and work with the
## predicted values that are relevant to a problem.  It is not more
## easier to interpret the EXACT SAME NUMBER from one model or the
## other.




m2rc <- residualCenter(m2)
summary(m2rc)

op <- par(no.readonly = TRUE)
par(mar = c(2,2,2,1), mfcol = c(2,2))

plotPlane(m1, "x1", "x2", plotPoints = TRUE, theta = -25,
          main = "No Interaction", ticktype = "detailed")

plotPlane(m2, "x1", "x2", plotPoints = TRUE, theta = -25,
          main = "Not Centered", ticktype = "detailed")

plotPlane(m2rc, "x1", "x2", plotPoints = TRUE, theta = -25,
          main = "Residual Centered", ticktype = "detailed")

plotPlane(m2mc, "x1", "x2", plotPoints = TRUE, theta = -25,
          main = "Mean Centered", ticktype = "detailed")

par(mfcol = c(1,1))


## I had trouble believing those plots. Is it possible that the
## predicted values from the residual centered regression are exactly the
## same as the predictions from the orinary non-centered m2 as well as
## the mean centered m2mc?  That seemed crazy, so I decided to
## double-check by manually calculating a residual-centered regression
## and extracting predicted values.

rcreg <-  lm ( I(x1*x2)  ~ x1 + x2, data = dat3)
rcx1x2 <- resid(rcreg)
m2manualrc <- lm(y ~ x1 + x2 + rcx1x2, data = dat3)
predm2manualrc <- predict(m2manualrc, newdata = dat3)



m2 <- lm( y ~ x1 * x2, data = dat3)
m2mc <- meanCenter(m2)
m2rc <- residualCenter(m2)

m2pred <- predict(m2, newdata = dat3)
dat3mc <- dat3
dat3mc$x1 <- dat3mc$x1c
dat3mc$x2 <- dat3mc$x2c
m2mcpred <- predict(m2mc, newdata = dat3mc)
m2rcpred <- predict(m2rc, newdata = dat3)

dat4 <- data.frame("m2pred"= m2pred, "m2mcpred" = m2mcpred,
                   "m2mcmaunal"= m2mcmanualpred, "m2rcpred" = m2rcpred,
                   "m2rcpred2" = predm2manualrc)

head(dat4)

cor(dat4)

## Well, the predicted values are the same.
## and the coefficients on the interaction terms are all the same.
## We already saw that the ordinary  interaction and the
## mean-centered regressions are identical.  Is it possible
## the residual centered version is just another identical
## model viewed from yet another point in the x1,x2 plane?
## Notice:
## m2
## lm(formula = y ~ x1 * x2, data = dat3)

## Coefficients:
## (Intercept)           x1           x2        x1:x2
##    -74.4572       1.5726       0.6462      -0.2053


## From rockchalk:::residualCenter(m2):
## m2rc
## lm(formula = y ~ x1 + x2 + x1.X.x2, data = dat)

## Coefficients:
## (Intercept)           x1           x2      x1.X.x2
##    424.1485      -8.9042      -9.2701      -0.2053



## From rockchalk:::meanCenter(m2)
## > m2mc
## The centered variables are:
## [1] "x1" "x2"
## The call that requested centering was:
## meanCenter.default(model = m2)

## Call:
## lm(formula = y ~ x1 * x2, data = stddat)

## Coefficients:
## (Intercept)           x1           x2        x1:x2
##   -456.6985      -8.3838      -9.4797      -0.2053



## If mean-centering is identical to residual centering,
## there must be some matrix algebra to prove it.

## b-hat = (X'X)-1 X'y

##          X = |1 x1 x2 x1*x2|

## mean centered:
##          X = |1  x1-m(x1) x2-m(x2)  (x1-m(x1))(x2-m(x2))|


## let pred(x1*x2) be the predicted value from the regression
## (x1*x2) = b0 + b1*x1 + b2*x2

## Then the design matrix for
## residual centered:
##          X = |1  x1 x2  (x1*x2)-pred(x1*x2) |

## I think if I were good at matrix algebra with partitioned
## matrices, I could demonstrate why the estimated coefficients
## are actually equivalent.
