transfH.fun <-
function(mlePP)
{
inddat<-mlePP@inddat
posE<-mlePP@posE
lambdafit<-mlePP@lambdafit
lambdafitc<-lambdafit*inddat

Ilambda <- inddat*cumsum(lambdafitc)
posEH <- Ilambda[posE]-lambdafitc[posE]

return(list(posEH=posEH,posE=posE, lambdafit=lambdafit, inddat=inddat))
}
