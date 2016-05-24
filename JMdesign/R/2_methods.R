setMethod("show",
          signature = "powerLongSurv",
          definition = 
            function(object){
              cat(object@title,"\n", object@subtitle,"\n")
              cat("Covariance matrix SigmaTheta :\n")
              print(object@SigmaTheta)
              cat("Measurement times: ", "t[1:",length(object@t),"]=", object@t, "\n",
                  "Subject proportions: ",  "p[1:",length(object@p),"]=",object@p, "\n",
                  "Number of subjects (N): ",object@N, "\n",
                  "Number of events (nevents): ",object@nevents, "\n",
                  "Percent censored (censr): ", object@censr, "\n",
                  "Median survival time (tmedian): ", object@tmedian, "\n",
                  "Mean follow-up time (meantf): ", object@meantf, "\n",
                  "Trajectory effect (beta): ", object@beta, "\n",
                  "Type I Error (alpha): ", object@alpha, "\n",
                  "Power (power): ", object@power, "\n")
            }
)


