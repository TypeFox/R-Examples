# plotregimeid() : clustering and plotting function for msbvar permuted
#                   sample output
#
# Iliyan Iliev and Patrick T. Brandt
#
# 20101115 : Initial version
# 20101124 : Added more plots and first logical switches for types of
#            plots for posterior identification
# 20101130 : Added labelings for plots; begin adding logical
#            conditions for the plot types.
# 20110118 : Updated logical conditions for plots and added section
#            for plotting versus values of Q.
# 20120426 : Updated to add AR(1) density and traceplots for the
#            permuted, clustered regimes.
# 20140609 : Added a par() reset so future plot calls have the right setup

plotregimeid <- function(x,
                         type=c("all", "intercepts", "AR1", "Sigma", "Q"),
                         ask = TRUE, ...)
{
    # Plot setup stuff
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    devAskNewPage()
    # Get the constants we need for setting up the matrices for the
    # clustering steps.

    m <- x$m                     # Number of equations
    p <- x$p                     # Lag length
    h <- x$h                     # Number of regimes
    N2 <- length(x$ss.sample)   # Number of Gibbs draws

    mpplus1 <- m*p + 1  # Number of coefficients in one equation of
                        # VAR(p)

    nc <- m*(mpplus1)   # Number of total coefs in one regime

    ii <- seq(mpplus1, by=mpplus1, length=m)  # Intercept indices
    aii <- rep(1:m,m) + (rep(0:(m-1), each=m))*mpplus1 # AR indices

    # Get eqn names
    eqnnames <- colnames(x$init.model$y)

    if(type=="all" || type=="intercepts")
    {
    # Now stack the regression coefficients for clustering
        Beta <- matrix(t(x$Beta.sample), byrow=TRUE, ncol=nc)

    # Now do the clustering for the regimes.
    # Probably should let the users choose the number of starting
    # points for the clustering centers = nstart, later.

        intercepts <- as.data.frame(Beta[,ii])
    # Label things so the plots make sense
        colnames(intercepts) <- eqnnames

        cl.int <- kmeans(intercepts, centers=h, nstart=10)


    # Plot the intercepts based on the clustering
        pairs(intercepts, pch=".", col=cl.int$cluster)
        title("Intercepts pairs by regime", line=3)

    # Note that the lattice / coda plotting function (the first two),
    # need to wrapped in a print to be displayed.

	form <- as.formula(paste("~",
                                 paste(lapply(names(intercepts), as.name),
                                       collapse = "+")))

        devAskNewPage(ask=ask)
        print(densityplot(form, data=intercepts, outer=TRUE,
                    groups=cl.int$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Intercept densities by regime"))

        form <- eval(parse(text =
                           paste(paste(lapply(names(intercepts), as.name),
                           collapse = "+"), "~ idx")))

        idx <- 1:nrow(intercepts)

        devAskNewPage(ask=ask)
        print(xyplot(form, data=intercepts,
                     groups=cl.int$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Intercept traceplots by regime"))
    }


    if(type=="all" || type=="AR1")
    {
    # Now stack the regression coefficients for clustering
        Beta <- matrix(t(x$Beta.sample), byrow=TRUE, ncol=nc)

    # Now do the clustering for the regimes.
    # Probably should let the users choose the number of starting
    # points for the clustering centers = nstart, later.

        ar1 <- as.data.frame(Beta[,aii])

        # Set up names
        idx <- expand.grid(1:m, 1:m)
        tmp <- cbind(rep(eqnnames,m), idx)

        colnames(ar1) <- paste(tmp[,1], "(", tmp[,2], ",", tmp[,3],
                               ")", sep="")

        # Do the clustering
        cl.ar1 <- kmeans(ar1, centers=h, nstart=10)


        form <- as.formula(paste("~",
                                 paste(lapply(names(ar1), as.name),
                                       collapse = "+")))
        devAskNewPage(ask=ask)
        print(densityplot(form, data=ar1, outer=TRUE,
                    groups=cl.ar1$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="AR(1) densities by regime"))

        form <- eval(parse(text =
                           paste(paste(lapply(names(ar1), as.name),
                           collapse = "+"), "~ idx")))

        idx <- 1:nrow(ar1)

        devAskNewPage(ask=ask)
        print(xyplot(form, data=ar1,
                     groups=cl.ar1$cluster, type="l", ylab=NULL,
                     main="AR(1) coefficients traceplots by regime"))

    }


    if(type=="all" || type=="Sigma")
    {

    # Now replicate this for the variance elements in
    # x$Sigma.sample.  Note that the unique elements of the Sigma(h)
    # matrices are in this and make up a different dimension than the
    # regression effects.
        nvar <- (m*(m+1)*0.5) # Number of variance terms
        Sigmaout <- matrix(t(x$Sigma.sample),
                           byrow=TRUE, ncol=nvar)

    # Get variance indices
        tmp <- m:2
        tmp <- rep(c(1,tmp), h)
        ssidx <- tmp
        for(i in 2:(m*h)) ssidx[i] <- ssidx[i-1] + tmp[i]

        Sigmaout <- as.data.frame(Sigmaout[,ssidx[1:m]])
        colnames(Sigmaout) <- eqnnames

        cl.sigma <- kmeans(Sigmaout, centers=h, nstart=10)

    # Plots
        devAskNewPage(ask=ask)
        pairs(Sigmaout, pch=".", col=cl.sigma$cluster, ...)
        title("Variances pairs plot by regime", line=3)

        form <- as.formula(paste("~",
                          paste(lapply(names(Sigmaout), as.name),
                                collapse = "+")))

        devAskNewPage(ask=ask)
        print(densityplot(form, data=Sigmaout, outer=TRUE,
                    groups=cl.sigma$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Variance densities by regime"), ...)

        form <- eval(parse(text = paste(paste(lapply(names(Sigmaout), as.name),
                           collapse = "+"), "~ idx")))
        idx <- 1:nrow(Sigmaout)

        devAskNewPage(ask=ask)
        print(xyplot(form, data=Sigmaout,
                     groups=cl.sigma$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Variance traceplots by regime"), ...)
    }

    if(type=="all" || type=="Q")
    {
    # Do the same as the above for the elements of x$Q.sample
    # These do not need to be stacked like the other elements (why?)

        Q <- as.data.frame(x$Q.sample)

    # Cluster

        cl.Q <- kmeans(Q, centers=h, nstart=10)

        idx <- cbind(rep(1:h, each=h), rep(1:h, times=h))
        Qnames <- paste("Q_", idx[,1], idx[,2], sep="")
        colnames(Q) <- Qnames

    # Plots
        devAskNewPage(ask=ask)
        pairs(Q, pch=".", col=cl.Q$cluster, ...)
        title("Transitions pairs plot by regime", line=3)

	form <- as.formula(paste("~",
                          paste(lapply(names(Q), as.name),
                                collapse = "+")))

        devAskNewPage(ask=ask)
        print(densityplot(form, data=Q, outer=TRUE,
                    groups=cl.Q$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Transition densities by regime"), ...)

        form <- eval(parse(text = paste(paste(lapply(names(Q), as.name),
                           collapse = "+"), "~ idx")))
        idx <- 1:nrow(Q)

        devAskNewPage(ask=ask)
        print(xyplot(form, data=Q,
                     groups=cl.Q$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Transition traceplots by regime"), ...)


    # Combine all of the above into one set of plots / analyses across
    # intercepts -- like SFS 2001, Figure 8

    # Plots based on the posterior clustering of Q for the intercepts
       # Now stack the regression coefficients for clustering
        Beta <- matrix(t(x$Beta.sample), byrow=TRUE, ncol=nc)

    # Now do the clustering for the regimes.
    # Probably should let the users choose the number of starting
    # points for the clustering centers = nstart, later.

        intercepts <- as.data.frame(Beta[,ii])
    # Label things so the plots make sense
        colnames(intercepts) <- eqnnames

    # Get diagonal of Q
        qdiag <- diag(matrix(1:h^2, h, h))

        devAskNewPage(ask=ask)
        par(mfrow=c(2, round(m/2)), omi=c(0.5, 0.75, 0.75, 0.25))
        for(i in 1:m)
        {
            plot(intercepts[,i], matrix(unlist(Q[, qdiag]), ncol=1), pch=".",
                 col=cl.Q$cluster, xlab=names(intercepts)[i],
                 ylab="Transition Probability Regimes")
        }
        title("Intercepts by transition probability regimes",
              outer=TRUE, line=1)

    }
    # Clean up after plotting
    devAskNewPage(FALSE)
    invisible()
}
