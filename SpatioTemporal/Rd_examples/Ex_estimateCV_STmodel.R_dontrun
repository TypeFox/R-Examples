##load data
data(mesa.model)
data(est.mesa.model)

################
## estimateCV ##
################
##create the CV structure defining 10 different CV-groups
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)

##use the best parameters and there starting values as
x.init <- coef(est.mesa.model, pars="cov")[,c("par","init")]

\dontrun{
  ##estimate different parameters for each CV-group
  est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)
}
##lets load precomputed results instead
data(est.cv.mesa)

##examine the estimation results
print( est.cv.mesa )
##estimated parameters for each CV-group
coef(est.cv.mesa, pars="cov")

###############
## predictCV ##
###############
\dontrun{
  ##Do cross-validated predictions using the just estimated parameters
  ##Ind.cv is infered from est.cv.mesa as est.cv.mesa$Ind.cv
  pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, LTA=TRUE)
}
##lets load precomputed results instead
data(pred.cv.mesa)

##prediction results
print( pred.cv.mesa )

##and CV-statistics
print( summary( pred.cv.mesa, LTA=TRUE) )


\dontrun{
  ##A faster option is to only consider the observations and not to compute
  ##variances
  pred.cv.fast <- predictCV(mesa.model, est.cv.mesa, only.obs=TRUE,
                            pred.var=FALSE)
  print( pred.cv.fast )
  summary( pred.cv.fast )
}
