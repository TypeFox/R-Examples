#################################################################################
#                                                                               #
#   Demonstration of the use the functions of the package "georob" for robust   #
#   geostatistical analyses                                                     #
#                                                                               #
#################################################################################

library( georob )

data( meuse )
data( meuse.grid )

# creation of a data set with replicated observations for some locations

set.seed( 1 )
meuse2 <- rbind( meuse, meuse )
meuse2[1:nrow( meuse ), "zinc"] <- exp( 
  log( meuse2[1:nrow( meuse ), "zinc"] ) + 
  rnorm( nrow( meuse ), sd = sqrt( 0.05) ) 
)

## 
 # 1. Fitting models, "meuse" data
 ##

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1000,
    control = control.georob(cov.bhat = TRUE, cov.ehat.p.bhat = TRUE))
summary(r.logzn.reml, correlation = TRUE)

logLik(r.logzn.reml)

waldtest(r.logzn.reml, .~. + ffreq)

## robust REML fit 
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)
    
summary(r.logzn.rob, correlation = TRUE)

plot(r.logzn.reml, lag.dist.def = seq( 0, 2000, by = 100 ))
lines(r.logzn.rob, col = "red")


# Gaussian REML Fit, data set with replicated observations
my.formula <- log( zinc ) ~ sqrt( dist )
r.logzn.reml2 <- georob(
  formula = my.formula,
  data = meuse2, 
  locations = ~ x + y,
  variogram.model = "RMexp",
  param = c( variance = 0.16, nugget = 0.03, scale = 208, snugget = 0.1 ),
  fit.param = c( variance = TRUE, nugget = TRUE, scale = TRUE, snugget = TRUE ),
  tuning.psi = 1000,
  verbose = 4
)
summary( r.logzn.reml2 )
plot( r.logzn.reml2, lag.dist.def = seq( 0, 2000, by = 100 ) )


# robust REML Fit, data set with replicated observations
r.logzn.rob2 <- georob(
  formula = my.formula,
  data = meuse2, 
  locations = ~ x + y,
  variogram.model = "RMexp",
  param = c( variance = 0.15, nugget = 0.02, scale = 208, snugget = 0.05 ),
  fit.param = c( variance = TRUE, nugget = TRUE, scale = TRUE, snugget = TRUE ),
  tuning.psi = 2,
  control = control.georob( initial.param = FALSE ),
  verbose = 2
)
summary( r.logzn.rob2 )
lines( r.logzn.rob2, col = "blue" )


## 
 # 2. Testing and residual diagnostics, "meuse" data set
 ##


# Wald-Test
t.georob <- update( r.logzn.rob, .~. + ffreq )

summary( t.georob )

waldtest( r.logzn.rob )
waldtest( t.georob, r.logzn.rob  )
waldtest( t.georob, .~.-ffreq )
waldtest( t.georob, "ffreq" )

waldtest( update( t.georob, .~.+soil ), t.georob )

# termplots
termplot( r.logzn.rob, partial = TRUE, se = TRUE )

## residual diagnostics
old.par <- par(mfrow = c(2,3))

plot(fitted(r.logzn.reml), rstandard(r.logzn.reml))
abline(h = 0, lty = "dotted")
qqnorm(rstandard(r.logzn.reml))
abline(0, 1)
qqnorm(ranef(r.logzn.reml, standard = TRUE))
abline(0, 1)
plot(fitted(r.logzn.rob), rstandard(r.logzn.rob))
abline(h = 0, lty = "dotted")
qqnorm(rstandard(r.logzn.rob))
abline(0, 1)
qqnorm(ranef(r.logzn.rob, standard = TRUE))
abline(0, 1)

par(old.par)

# display of robustness weights
plot( 
    y~x, meuse, 
    cex = sqrt( r.logzn.rob$rweights ) , asp = 1 
)


## 
 # 3. Cross-validation, "meuse" data set
 ##

r.cv.georob.reml<- cv( 
    r.logzn.reml, 
    seed = 1,
    lgn = TRUE, 
    return.fit = TRUE,
    verbose = 0
)
summary( r.cv.georob.reml )

r.cv.georob.rob <- cv( 
    r.logzn.rob, 
    seed = 1,
    lgn = TRUE, 
    return.fit = TRUE,
    verbose = 1
)
summary( r.cv.georob.rob )

# display of measured values vs. cross-validation predictions

# log scale
with( r.cv.georob.rob$pred, plot( data~pred ) ); abline( 0, 1)

# original scale
with( r.cv.georob.rob$pred, plot( lgn.data~lgn.pred ) ); abline( 0, 1)

#  Brier Score
plot( r.cv.georob.reml, "bs" )
plot( r.cv.georob.rob, "bs", col= "red", add = TRUE )



## 
 # 4. Kriging, "meuse" data set
 ##

# point kriging
r.luk.punkt <- predict(
    r.logzn.rob,
    type = "response",
    newdata = meuse.grid,
    control = control.predict.georob( extended.output = TRUE ),
    mmax = 1000,
    verbose = 1
)
str( r.luk.punkt )

# back-transformation
r.luk.punkt <- lgnpp( r.luk.punkt )
summary( r.luk.punkt )

# display
library( lattice )

levelplot( lgn.pred~x+y, r.luk.punkt )
f.colors <- colorRampPalette(c("darkblue",  "cyan", "yellow", "orange", "magenta"))
t.breaks <- c( seq( 0, 2000, by = 200 ), 2500, 3000, 3500 )

# predictions
levelplot( 
    lgn.pred ~ x + y, r.luk.punkt, 
    col.regions = f.colors(100), 
    at = t.breaks,
    aspect = "iso", 
    colorkey = list( 
        at = t.breaks, col = f.colors( length( t.breaks ) - 1 ),
        labels = list( at = t.breaks, labels = as.character( t.breaks ) )
    ),
    main = "lognormal kriging prediction zn content",
    panel = function( x, y, z, ..., xp, yp, zp, colp ){
        panel.levelplot( x, y, z, ... )
        panel.points(xp, yp, cex = zp, col = colp, lwd = 0.7 )
    },
    xp = meuse$x,
    yp = meuse$y,
    zp = sqrt( meuse$zinc )/10,
    colp = "grey"
)


# limits of prediction intervals
levelplot( 
    lgn.upper ~ x + y, r.luk.punkt, 
    col.regions = f.colors(100), 
    aspect = "iso", 
    at = t.breaks,
    colorkey = list( 
        at = t.breaks, col = f.colors( length( t.breaks ) - 1 ),
        labels = list( at = t.breaks, labels = as.character( t.breaks ) )
    ),
    main = "upper limit 95% prediction interval zn content"
)


levelplot( 
    lgn.lower ~ x + y, r.luk.punkt, 
    col.regions = f.colors(100), 
    aspect = "iso", 
    at = t.breaks,
    colorkey = list( 
        at = t.breaks, col = f.colors( length( t.breaks ) - 1 ),
        labels = list( at = t.breaks, labels = as.character( t.breaks ) )
    ),
    main = "lower limit 95% prediction interval zn content"
)

# standard error
t.breaks.se <- c( seq( 0, 300, by = 30 ), 400, 500, 600 )
levelplot( 
    lgn.lower ~ x + y, r.luk.punkt, 
    col.regions = f.colors(100), 
    aspect = "iso", 
    at = t.breaks.se,
    colorkey = list( 
        at = t.breaks, col = f.colors( length( t.breaks ) - 1 ),
        labels = list( at = t.breaks, labels = as.character( t.breaks.se ) )
    ),
    main = "prediction standard error zn content"
)



# block kriging
library(constrainedKriging)

r.luk.block <- predict(
  r.logzn.rob, newdata = meuse.blocks, 
  control = control.predict.georob( extended.output = TRUE, pwidth = 75, pheight = 75 )
)
str( r.luk.block, max = 2 )
str( r.luk.block@data, max = 2 )


# back-transformation under assumption of permanence of lognormality
r.luk.block <- lgnpp(r.luk.block, newdata = meuse.grid)
str( r.luk.block@data, max = 2 )


# display

# predictions
spplot( 
    r.luk.block, zcol = "lgn.pred", col.regions = f.colors(100), 
    at = t.breaks, main = "lognormal block prediction zn content"
)

# standard error
spplot( 
    r.luk.block, zcol = "lgn.se", col.regions = f.colors(100), 
    at = t.breaks.se, main = "standard error block prediction zn content"
)


## 
 # 5. Fitting models to "wolfcamp" data 
 ##

library(geoR)
data(wolfcamp)
d.wolfcamp <- data.frame(x = wolfcamp[[1]][,1], y = wolfcamp[[1]][,2],
    pressure = wolfcamp[[2]])

# fitting an isotropic IRF(0) model
r.irf0.iso <- georob(pressure ~ 1, data = d.wolfcamp, locations = ~ x + y, 
    variogram.model = "RMfbm",
    param = c( variance = 10, nugget = 1500, scale = 1, alpha = 1.5 ),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    tuning.psi = 1000)
  
summary(r.irf0.iso)

# fitting an isotropic IRF(0) model
r.irf0.aniso <- georob(pressure ~ 1, data = d.wolfcamp, locations = ~ x + y, 
    variogram.model = "RMfbm",
    param = c( variance = 5.9, nugget = 1450, scale = 1, alpha = 1 ),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    aniso = c( f1 = 0.51, f2 = 1, omega = 148, phi = 90, zeta = 0 ),
    fit.aniso = c( f1 = TRUE, f2 = FALSE, omega = TRUE, phi = FALSE, zeta = FALSE ), 
    tuning.psi = 1000)
summary(r.irf0.aniso)

plot(r.irf0.iso, lag.dist.def = seq(0, 200, by = 7.5))
plot(r.irf0.aniso, lag.dist.def = seq(0, 200, by = 7.5), 
    xy.angle.def = c(0, 22.5, 67.5, 112.5, 157.5, 180.), 
    add = TRUE, col = 2:5)
    
pchisq( 2*(r.irf0.aniso$loglik - r.irf0.iso$loglik), 2, lower = FALSE )



## 
 # 6.  Computing sample variogram and fitting variogram model to it,
 #     "wolfcamp" data
 ##

data(wolfcamp)

# fitting an isotropic IRF(0) model
r.sv.iso <- sample.variogram(wolfcamp$data,
    locations = wolfcamp[[1]], lag.dist.def = seq(0, 200, by = 15))

r.irf0.iso <- fit.variogram.model(r.sv.iso, variogram.model = "RMfbm",
    param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE ),
    method = "Nelder-Mead", hessian = FALSE, 
    control = list(maxit = 5000), verbose = 0)  
summary(r.irf0.iso, correlation = TRUE)

plot( r.sv.iso, type = "l")
lines( r.irf0.iso, line.col = "red")

# fitting an anisotropic IRF(0) model
r.sv.aniso <- sample.variogram(wolfcamp$data,
    locations = wolfcamp[[1]], lag.dist.def = seq(0, 200, by = 15),
    xy.angle.def = c(0., 22.5, 67.5, 112.5, 157.5, 180.))
summary(r.sv.aniso)

r.irf0.aniso <- fit.variogram.model(r.sv.aniso, variogram.model = "RMfbm",
    param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1.5),
    fit.param = c(variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE ),
    aniso = c(f1 = 0.4, f2 = 1., omega = 135, phi = 90., zeta = 0.),
    fit.aniso = c(f1 = TRUE, f2 = FALSE, omega = TRUE, phi = FALSE, zeta = FALSE),
    method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 5000), verbose = 0)
summary(r.irf0.aniso, correlation = TRUE)

plot(r.sv.aniso, type = "l")
lines(r.irf0.aniso, xy.angle = seq( 0, 135, by = 45))
