filter.NULL <-
function (...) 
{
    UseMethod("filter.NULL")
}
filter.PCA <-
function (X, nbreVarX_, ...) 
{
    UseMethod("filter.PCA")
}
filter.RegressionTreeFilter <-
function (X, nbreVarX_, ...) 
{
    UseMethod("filter.RegressionTreeFilter")
}
filter.mRMR <-
function (X, Y, nbreVarX_, ...) 
{
    UseMethod("filter.mRMR")
}
filter.MAX <-
function (nbreVarX_, ...) 
{
    UseMethod("filter.MAX")
}