"gibbs.A0" <- function(varobj, N1, N2, thin=1, normalization="DistanceMLA")
{
    # Leave the character strings on the R side and pass an int flag
    # to the normalization routine norm_svar (fully documented in
    # MSBVARfun.cpp) in the order below

    tmp <- sanity.check.gibbs(list(N1=N1, N2=N2, thin=thin,
                                   normalization=normalization))

    methodlist <- c("DistanceMLA", "DistanceMLAhat", "Euclidean",
                    "PositiveDiagA", "PositiveDiagAinv")

    if(tmp){ method <- which(methodlist==normalization)-1 } else
    {method <- which(methodlist==tmp)-1}

    cat("Normalization Method: ", normalization, "(", method,")\n")

    tmp2 <- .Call("gibbsA0.cpp", varobj, as.integer(N1), as.integer(N2),
                  as.integer(thin), as.integer(method),
                  gibbs.setup.bsvar(varobj)$UT)

    # Memory cleanup
    gc(); gc();

    # Set classing
    class(tmp2) <- c("gibbs.A0")
    return(tmp2)
}


# Converts the A0 object into something that coda can understand.
"A02mcmc" <- function(x)
{
    return(mcmc(matrix(x$A0.posterior$A0,
                       nrow=x$N2,
                       ncol=length(x$A0.posterior$struct),
                       byrow=T),
                thin=x$thin))
}
