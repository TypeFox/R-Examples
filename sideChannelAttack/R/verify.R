verify.loo <-
function (model, filter, X, Y, nbreVarX, param.model=list(), param.fs=list(), ...) 
{
    UseMethod("verify.loo")
}
verify.ho <-
function (model, filter, Xlearn, Ylearn, Xval, Yval, nbreVarX, param.model=list(), param.fs=list(), ...) 
{
    UseMethod("verify.ho")
}
verify.cv <-
function (model, filter, X, Y, nbreVarX, k, param.model=list(), param.fs=list(), ...) 
{
    UseMethod("verify.cv")
}