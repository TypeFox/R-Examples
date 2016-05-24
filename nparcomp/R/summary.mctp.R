summary.mctp <-
function(object,...)
{           
            cat("\n", "#----------------Nonparametric Multiple Comparisons for relative effects---------------#", "\n","\n",
        "-", "Alternative Hypothesis: ", object$text.Output,"\n",
        "-", "Estimation Method: Global Pseudo ranks","\n",
        "-", "Type of Contrast", ":", object$input$type, "\n", "-", "Confidence Level:",
            object$input$conf.level*100,"%", "\n", "-", "Method", "=", object$AsyMethod, "\n","\n",
                  "#--------------------------------------------------------------------------------------#","\n",

            "\n")
            cat( " #----Data Info-------------------------------------------------------------------------#","\n")
            print(object$Data.Info)
            cat("\n", "#----Contrast--------------------------------------------------------------------------#","\n")
            print(object$Contrast)
            cat("\n", "#----Analysis--------------------------------------------------------------------------#","\n")
            print(object$Analysis)
            cat("\n", "#----Overall---------------------------------------------------------------------------#","\n")
            print(object$Overall)
            if (object$input$correlation == TRUE){
            cat("\n", "#----Covariance------------------------------------------------------------------------#","\n")
            print(object$Covariance)
            cat("\n", "#----Correlation-----------------------------------------------------------------------#","\n")
            print(object$Correlation)
            }
            cat("\n", "#--------------------------------------------------------------------------------------#","\n")
}
