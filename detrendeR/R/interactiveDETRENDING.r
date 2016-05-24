interactiveDETRENDING =
function (input = NULL, output = NULL, TwoSteps = TRUE, ...) 
{
    if (makeFirstDetrending) {
        eval(parse(text = paste("assign('", output, ".cv1', apply(", 
            input, ", 2, detrendeR:::RemoveTrend, method=method1, BandwidthPerc=nPerc1, Bandwidth = n1, P=p1))", 
            sep = "")), envir = .GlobalEnv)
        cat("\nDETRENDING [", input, "]\n\nFirst detrending [", 
            first.detrending.method, "]\n", sep = "")
        if (interactive.detrend) {
            eval(parse(text = paste("detrendeR:::InteractiveDetrending(", 
                input, ", Detrend =", output, ".cv1,  method=method1,  n=n1, nPerc=nPerc1, p=p1)", 
                sep = "")), envir = .GlobalEnv)
            if (detrendeR:::.get("DETRENDING_INTERACTIVE_FLAG")) {
                eval(parse(text = paste(output, ".cv1<<-detrendeR:::.get('DETRENDING_INTERACTIVE_OUTPUT')", 
                  sep = "")), envir = .GlobalEnv)
            }
        }
        eval(parse(text = paste(output, ".in1<<-", input, "/", 
            output, ".cv1", sep = "")), envir = .GlobalEnv)
        eval(parse(text = paste("detrendeR:::RwlInfo(", output, 
            ".in1, print=TRUE)", sep = "")), envir = .GlobalEnv)
        if (TwoSteps) {
            cat("\nSecond detrending [", second.detrending.method, 
                "]\n", sep = "")
            eval(parse(text = paste("assign('", output, ".cv2', apply(", 
#                input, ", 2, detrendeR:::RemoveTrend, method=method2, BandwidthPerc=nPerc2, Bandwidth = n2, P=p2))", 
				output, ".in1, 2, detrendeR:::RemoveTrend, method=method2, BandwidthPerc=nPerc2, Bandwidth = n2, P=p2))", 

                sep = "")), envir = .GlobalEnv)
            if (interactive.detrend) {
                eval(parse(text = paste("detrendeR:::InteractiveDetrending(", 
                  output, ".in1, Detrend =", output, ".cv2, folderName = \"SecondDetrend\", method=method2,  n=n2, nPerc=nPerc2, p=p2)", 
                  sep = "")), envir = .GlobalEnv)
                if (detrendeR:::.get("DETRENDING_INTERACTIVE_FLAG")) {
                  eval(parse(text = paste(output, ".cv2<<-detrendeR:::.get('DETRENDING_INTERACTIVE_OUTPUT')", 
                    sep = "")), envir = .GlobalEnv)
                }
            }
            eval(parse(text = paste(output, ".in2<<-", output, 
                ".in1/", output, ".cv2", sep = "")), envir = .GlobalEnv)
            eval(parse(text = paste("detrendeR:::RwlInfo(", output, 
                ".in2, print=TRUE)", sep = "")), envir = .GlobalEnv)
        }
    }
}


