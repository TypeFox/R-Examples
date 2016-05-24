"examineFit" <-
function (resultfitModel, opt=vector()) 
{
     if(length(opt) == 0)
	plotoptions <-    resultfitModel$currModel@optlist[[1]]
     else 
	plotoptions <- opt

     plotter(resultfitModel$currModel@modellist[[1]],
     resultfitModel$currModel, 
     resultfitModel$currTheta,
     plotoptions)

}
