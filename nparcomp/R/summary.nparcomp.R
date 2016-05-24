summary.nparcomp <-
function(object,...)
{   
            cat("\n", "#------------Nonparametric Multiple Comparisons for relative contrast effects----------#", "\n","\n",
        "-", "Alternative Hypothesis: ", object$text.Output,"\n",
        "-", "Estimation Method: Global Pseudo ranks","\n",
        "-", "Type of Contrast", ":", object$input$type, "\n", "-", "Confidence Level:",
            object$input$conf.level*100,"%", "\n", "-", "Method", "=", object$AsyMethod, "\n","\n",
        "-", "Estimation Method: Pairwise rankings","\n", "\n",
                    "#---------------------------Interpretation--------------------------------------------#",
            "\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a","\n",
            "#-------------------------------------------------------------------------------------#", "\n",

            "\n")
            cat( " #----Data Info-------------------------------------------------------------------------#","\n")
            print(object$Data.Info)
            cat("\n", "#----Contrast--------------------------------------------------------------------------#","\n")
            print(object$Contrast)
            cat("\n", "#----Analysis--------------------------------------------------------------------------#","\n")
            print(object$Analysis)
            cat("\n", "#----Overall---------------------------------------------------------------------------#","\n")
            print(object$Overall)
            if (object$input$weight.matrix == TRUE){
            cat("\n", "#----Weight Matrix---------------------------------------------------------------------#","\n")
            print(object$Weight.Matrix)
            cat("\n", "#----All Pairs-------------------------------------------------------------------------#","\n")
            print(object$AllPairs)
            }
            if (object$input$correlation == TRUE){
            cat("\n", "#----Covariance------------------------------------------------------------------------#","\n")
            print(object$Covariance)
            cat("\n", "#----Correlation-----------------------------------------------------------------------#","\n")
            print(object$Correlation)
            }
            cat("\n", "#--------------------------------------------------------------------------------------#","\n")
}
