`ensembleVerifObs.ensembleData` <-
function (x) 
{
 class(x) <- "data.frame" 
 x$obs
}

