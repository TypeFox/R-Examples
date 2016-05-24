contrMat <- function(n, type = c("Dunnett", "Tukey", "Sequen", "AVE",
                                 "Changepoint", "Williams", "Marcus",
                                 "McDermott", "UmbrellaWilliams", "GrandMean"),
                     base = 1) {

    if (length(n) < 2) stop("less than two groups")
    if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
    m <- NULL
    type <- match.arg(type)
    if (type %in% c("AVE", "Williams", "McDermott") && length(n) < 3)
        stop("less than three groups")
    k <- length(n)
    if (base < 1 || base > k) stop("base is not between 1 and ", k)
    CM <- c()
    rnames <- c()
    if (!is.null(names(n)))
        varnames <- names(n)
    else
        varnames <- 1:length(n)

    kindx <- 1:k

    switch(type, "Dunnett" = {
        for(i in kindx[-base])
            CM <- rbind(CM, as.numeric(kindx == i) - as.numeric( kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base])
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                CM  <- rbind(CM, as.numeric(kindx==j)-as.numeric(kindx==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i]))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            CM  <- rbind(CM, as.numeric(kindx==i)-as.numeric(kindx==i-1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i-1]))
        }
    }, "AVE" = {
        help <- c(1,  -n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        for (i in 2:(k-1)) {
            x <- sum(n[1:(i-1)])+sum(n[(i+1):k])
            help <- c(-n[1:(i-1)]/x, 1, -n[(i+1):k]/x)
            CM <- rbind(CM, help)
        }
        help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
        CM  <- rbind(CM, help)
        rnames <- paste("C", 1:nrow(CM))
    }, "Changepoint" = {
        for (i in 1:(k-1)) {
            help <- c(-n[1:i]/sum(n[1:i]), n[(i+1):k]/sum(n[(i+1):k]))
            CM <- rbind(CM, help)
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Williams" = {
        for (i in 1:(k-2)) {
            help <-  c(-1, rep(0, k-i-1), n[(k-i+1):k]/sum(n[(k-i+1):k]))
            CM <- rbind(CM, help)
        }
        help <- c(-1, n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Marcus" = {
        cm1 <- matrix(0, nrow=k-1, ncol=k)
        cm2 <- cm1
        for (i in 1:(k-1)) {
            cm1[i,(i+1):k] <- n[(i+1):k]/sum(n[(i+1):k])
            cm2[i,1:i] <- n[1:i]/sum(n[1:i])
        }
        ### row <- k*(k-1)/2
        index <- 1
        for (i in 1:(k-1)) {
            for (j in 1:i) {
                help <- cm1[i,]-cm2[j,]
                CM <- rbind(CM, help)
                index <- index+1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
     }, "McDermott" = {
         for(i in 1:(k-2)) {
             help  <- c(-n[1:i]/sum(n[1:i]), 1, rep(0, k-i-1))
             CM <- rbind(CM, help)
         }
         help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
         CM  <- rbind(CM, help)
         rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Tetrade" = {
        if (is.null(m)) stop(sQuote("m"), " is missing")
        a <- length(n)
        b <- length(m)
        if (!is.null(names(m)))
            varnamesm <- names(m)
        else
            varnamesm <- 1:length(m)
	idi <- 1:a
	idj <- 1:b
        for (i1 in 1:(a-1)) {
            for (i2 in (i1+1):a) {
	        for (j1 in 1:(b-1)) {
        	    for (j2 in (j1+1):b) {
                	CM <- rbind(CM, kronecker( ( as.numeric(idi==i1)-as.numeric(idi==i2) ),
                                                   ( as.numeric(idj==j1)-as.numeric(idj==j2) ) ) )
		        rnames <- c(rnames, paste( "(", paste(varnames[i1], varnamesm[j1], sep = ":"), "-",
                                                        paste(varnames[i1], varnamesm[j2], sep = ":"), ")", "-",
                                                   "(", paste(varnames[i2], varnamesm[j1], sep = ":"), "-",
                                                        paste(varnames[i2], varnamesm[j2], sep = ":"), ")",  sep=""))
            	    }
        	}
	    }
        }
    }, "UmbrellaWilliams" = {
        for (j in 1:(k-1)) {
            for (i in 1:(k - j)) {
                helper <- c(-1, rep(0, k - i - j),
                    n[((k - i + 1):k)-(j-1)]/sum(n[((k - i + 1):k)-(j-1)]),
                    rep(0, j-1))
                CM <- rbind(CM, helper)
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "GrandMean" = {
        CM <- matrix(rep(-n/sum(n), k), nrow = k, byrow = TRUE)
        diag(CM) <- diag(CM) + 1
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    })

    rownames(CM) <- rnames
    if (type == "Tetrade")
      colnames(CM) <- NULL ###levels(interaction(varnames, varnamesm))
    else
      colnames(CM) <- varnames
    attr(CM, "type") <- type
    class(CM) <- c("contrMat", "matrix")
    CM
}
