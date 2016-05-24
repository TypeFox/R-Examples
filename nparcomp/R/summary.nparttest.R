summary.nparttest <-
function(object,...)
{         
            cat("\n", "#-----Nonparametric Test Procedures and Confidence Intervals for relative  effects-----#", "\n","\n",
        "-", "Alternative Hypothesis: ", object$text.Output,"\n",
        "-", "Confidence level:", object$input$conf.level*100,"%", "\n", "-", "Method", "=", object$AsyMethod,
                           "\n","#---------------------------Interpretation---------------------------------------------#",
            "\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a","\n",
                  "#--------------------------------------------------------------------------------------#","\n",

            "\n")
            cat( " #----Data Info-------------------------------------------------------------------------#","\n")
            print(object$Info)
            cat("\n", "#----Analysis--------------------------------------------------------------------------#","\n")
            print(object$Analysis)

}
