`print.reliability` <-
function(x, ...){   # define the characteristics of the output of the class "Reliability"
                   cat("\n Number of Items \n", unlist(x[1]),"\n")
                   cat("\n Number of Examinees \n",unlist(x[2]),"\n")
                   cat("\n Coefficient Alpha \n", round(unlist(x[3]),3),"\n")
                   invisible(x)}

