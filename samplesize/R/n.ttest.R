n.ttest <-
function(power = 0.80, alpha = 0.05, mean.diff = 0.80, sd1 = 0.83, sd2 = 2.65, k = 1, design = "unpaired", fraction = "balanced", variance = "equal")
{
    if(fraction == "unbalanced" & k == 1){
        warning("Groups are chosen unbalanced, but fraction argument k is 1")
    }
    if(design == "paired" & fraction == "unbalanced"){
        warning("Argument -unbalanced- is not used. Paired design is balanced")
    }
    if(design == "paired" & k != 1){
        warning("Argument -k- is set to 1. Paired design is balanced")
    }
    if(design == "paired" & variance == "unequal"){
        warning("Paired design assumes and uses equal variances")
    }
    if(design == "paired"){
        fraction = "balanced"
    }
    if(design == "paired"){
        variance = "equal"
    }
    if(design == "unpaired" & variance == "unequal"){
        warning("Arguments -fraction- and -k- are not used, when variances are unequal")
    }
    if(power > 1 | power < 0){
        stop("Power must be between 0 and 1.0")
    }
    if(power < 0.5){
        warning("Are you sure that Power should be lower than 50 % ?")
    }
    if(alpha > 1 | alpha < 0){
        stop("Type-I-error must be between 0 and 1.0")
    }
    if(alpha > 0.1){
        warning("Are you sure that the two-sided Type-I-Error should be larger than 10 % ?")
    }
    if(k < 0){
        stop("Fraction k must be greater than zero")
    }
    conf.level <- 1 - alpha / 2
    n.start <- 4
    switch(variance,

           "unequal" = {
               k <- sd2/sd1
               n1.pri <- n.start/(1 + k)
               n2.pri <- (k * n.start)/(1 + k)
               n1 <- max(n1.pri, 2)
               n2 <- max(n2.pri, 2)
               gamma <- sd1/(sd1 + sd2)
               c <- mean.diff/(sd1 + sd2)
               df_approx <- 1/ ((gamma)^2/ (n1 - 1) + (1 - gamma)^2/ (n2-1))
               tkrit.alpha <- qt(conf.level, df = df_approx)
               tkrit.beta <- qt(power, df = df_approx)
               n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
               while(n.start <= n.temp){
                   n.start <- n1 + n2 + 1
                   n1 <- n.start/(1 + k)
                   n2 <- (k * n.start)/(1 + k)
                   df_approx <- 1/ ((gamma)^2/ (n1 - 1) + (1 - gamma)^2/ (n2-1))
                   tkrit.alpha <- qt(conf.level, df = df_approx)
                   tkrit.beta <- qt(power, df = df_approx)
                   n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
               }
               output <- list("Total sample size" = ceiling(n1) + ceiling(k * n1), "Sample size group 1" = ceiling(n1), "sample size group 2" = ceiling(n2))
               return(output)
           },

        "equal" = {
        {
            switch(
            design,

            "paired" = {
                n.start <- 2
                c <- mean.diff / sd1
                tkrit.alpha <- qt(conf.level, df = n.start -1)
                tkrit.beta <- qt(power, df = n.start - 1)
                n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
                while(n.start <= n.temp){
                    n.start <- n.start + 1
                    tkrit.alpha <- qt(conf.level, df = n.start - 1)
                    tkrit.beta <- qt(power, df = n.start - 1 )
                    n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
                }
                output <- list("Total sample size" = n.start)
                return(output)
            },

            "unpaired" = {
                switch(
                fraction,

                "balanced" = {
                    c <- mean.diff / (2*sd1)
                    tkrit.alpha <- qt(conf.level, df = n.start - 1)
                    tkrit.beta <- qt(power, df = n.start - 1)
                    n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)

                    while(n.start <= n.temp){
                        n.start <- n.start + 1
                        tkrit.alpha <- qt(conf.level, df = n.start - 1)
                        tkrit.beta <- qt(power, df = n.start - 1)
                        n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
                    }
                    n1 <- ceiling(n.start / 2)
                    n2 <- ceiling(n.start / 2)
                    output <- list("Total sample size" = 2 * n1, "Sample size group 1" = n1, "Sample size group 2" = n2)
                    return(output)

                },

                "unbalanced" = {
                    df <- n.start - 2
                    c <- (mean.diff/sd1)*(sqrt(k)/(1+k))
                    tkrit.alpha <- qt(conf.level, df = df)
                    tkrit.beta <- qt(power, df = df)
                    n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)

                    while(n.start <= n.temp){
                        n.start <- n.start + 1
                        tkrit.alpha <- qt(conf.level, df = n.start - 2)
                        tkrit.beta <- qt(power, df = n.start - 2)
                        n.temp <- ((tkrit.alpha + tkrit.beta)^2)/(c^2)
                    }
                    n1 <- n.start / (1 + k)
                    n2 <- k * n1
                    output <- list("Total sample size" = ceiling(n1) + ceiling(n2), "Sample size group 1" = ceiling(n1), "Sample size group 2" = ceiling(n2), "Fraction" = k)
                    return(output)
                })}
                )}
            return(output)
        }
        )
}
