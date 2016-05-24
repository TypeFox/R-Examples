plotSimScat <-
function(obj, sigma=NULL, layout=c(4,1), type=c("p","r"),
             show=c("points","residuals"), ...){
        nsim <- prod(layout)
        if(is.null(sigma))sigma <- summary(obj)[["sigma"]]
        hat <- fitted(obj)
        xnam <- all.vars(formula(obj))[2]
        ynam <- all.vars(formula(obj))[1]
        df <- data.frame(sapply(1:nsim, function(x)rnorm(length(hat), sd=sigma)))
        if(show[1]=="points")df <- df + hat
        simnam <- names(df) <- paste("Simulation", 1:nsim, sep="")
        df[, c(xnam, ynam)] <- model.frame(obj)[, c(xnam, ynam)]
        if(show[1]!="points"){df[, "Residuals"] <- df[, ynam] - hat
                              ynam <- "Residuals"
                              legadd <- "residuals"
                          } else legadd <- "data"
        leg <- list(text=paste(c("Simulated", "Actual"), legadd), columns=2)
        formula <- formula(paste(paste(simnam, collapse="+"), "~", xnam))
        parset <- simpleTheme(pch=c(16,16), lty=2, col=c("black","gray"))
        gph <- xyplot(formula, data=df, outer=TRUE, par.settings=parset,
                      auto.key=leg, lty=2, layout=layout, type=type, ...)
        formxy <- formula(paste(ynam, "~", xnam))
        addgph <- xyplot(formxy, data=df, pch=16, col="gray")
        gph+as.layer(addgph, under=TRUE)
    }
