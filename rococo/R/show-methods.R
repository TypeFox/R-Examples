setMethod("show", signature(object="RococoTestResults"),
    function(object)
    {
        cat("\n\tRobust Gamma Rank Correlation:\n\n")
        cat("data: ", object@input, " (length = ", object@length, ")\n", sep="")

        if (length(object@similarity) == 1 ||
            object@similarity[1] == object@similarity[2])
        {
            cat("similarity:", object@similarity[1], "\n")

            if (object@similarity[1] != "classical")
                cat("rx =", object@r.values[1], "/ ry =", object@r.values[2],
                    "\n")
        }
        else
        {
            cat(paste("similarity for x:", object@similarity[1]))
            if (object@similarity[1] != "classical")
                cat(" (rx = ", object@r.values[1], ")\n", sep="")
            else
                cat("\n")

            cat(paste("similarity for y:", object@similarity[2]))
            if (object@similarity[2] != "classical")
                cat(" (ry = ", object@r.values[2], ")\n", sep="")
            else
                cat("\n")
        }

        cat("t-norm:", object@tnorm$name, "\n")
        cat("alternative hypothesis:",
            switch(object@alternative,
                   two.sided="true gamma is not equal to 0",
                   greater="true gamma is greater than 0",
                   less="true gamma is less than 0"), "\n")
        cat("sample gamma =", object@sample.gamma, "\n")

        if (object@exact)
            cat("exact p-value = ",
                ifelse(object@p.value < .Machine$double.eps,
                       paste("<", format(.Machine$double.eps, digits=2)),
                       object@p.value),
                " (", object@count, " of ",
                object@numtests, " values)\n\n", sep = "")
        else
            cat("estimated p-value = ",
                ifelse(object@p.value < .Machine$double.eps,
                       paste("<", format(.Machine$double.eps, digits=2)),
                       object@p.value),
                " (", object@count, " of ",
                object@numtests, " values)\n\n", sep = "")
     }
)
