wald2ci <-
    function(x1, n1, x2, n2, conf.level, adjust = c("AC", "Wald")){
        switch(adjust,
               "AC" = {
                   n1m <- n1 + 2
                   n2m <- n2 + 2
                   p1hat <- (x1 + 1)/ (n1m)
                   p2hat <- (x2 + 1)/ (n2m)
               },
               "Wald" = {
                   n1m <- n1
                   n2m <- n2
                   p1hat = x1/ n1m
                   p2hat = x2/ n2m
               }
               )
        z = abs(qnorm((1-conf.level)/2))
        diff <- (p1hat - p2hat)
        ll = (p1hat - p2hat)  - z*sqrt((p1hat*(1-p1hat))/n1m + (p2hat*(1-p2hat))/n2m)
        ul = (p1hat - p2hat)  + z*sqrt((p1hat*(1-p1hat))/n1m + (p2hat*(1-p2hat))/n2m)
        cint <- c(ll, ul)
        attr(cint, "conf.level") <- conf.level
        rval <- list(conf.int = cint, estimate = diff)
        class(rval) <- "htest"
        return(rval)
    }

