summary.nparttestpaired <-
function(object,...)
{         
            cat("\n", "# Nonparametric Paired t Test Procedures and Confidence Intervals for the relative  effect #", "\n","\n",
        "-", "Alternative Hypothesis: ", object$text.Output,"\n",
        "-", "Confidence level:", object$input$conf.level*100,"%", "\n", "-", "Method", "=", object$Method,
                           "\n","#---------------------------Interpretation-------------------------------------------------#",
            "\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a","\n",
                  "#------------------------------------------------------------------------------------------#","\n",

            "\n")
            cat( " #----Data Info-----------------------------------------------------------------------------#","\n")
            print(object$Info)
            cat("\n", "#----Analysis------------------------------------------------------------------------------#","\n")
            print(object$Analysis)
            cat("\n", "#------------------------------------------------------------------------------------------#","\n")

}
