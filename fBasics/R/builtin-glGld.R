
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:
#  .pgld
#  .qgld
#   .qgl.fmkl
#   .qgl.fm5
#   .qgl.rs
#   .qdgld
#   .qdgl.rs
#   .qdgl.fmkl
#   .qdgl.fm5
#  .rgld
#  .gl.parameter.tidy
#  .gl.check.lambda
################################################################################


# Code borrowed from
#   R's contributed package "gld" written by Robert King


# Rmetrics:
#   Note that gmm is not available on Debian as of 2009-11-11.
#   To run these functions under Debian/Rmetrics we have them
#   implemented here as a builtin.


# Package: gld
# Version: 1.8.4
# Date: 2008/10/01
# Title: Estimation and use of the generalised (Tukey) lambda distribution
# Author: Robert King <Robert.King@newcastle.edu.au>
# Maintainer: Robert King <Robert.King@newcastle.edu.au>
# Description: The generalised lambda distribution, or Tukey lambda
#   distribution, provides a wide variety of shapes with one
#   functional form.  This package provides random numbers,
#   quantiles, probabilities, densities and plots.  It also
#   includes an implementation of the starship estimation method
#   for the distribution.
# License: GPL (>= 2)
# URL: http://tolstoy.newcastle.edu.au/~rking/publ/rprogs/information.html
# Packaged: 2009-10-14


################################################################################


.dgld <-
function(x, lambda1=0, lambda2=NULL, lambda3=NULL, lambda4=NULL, param="fmkl",
    lambda5=NULL, inverse.eps=1e-8, max.iterations=500)
{
    # Tidy the parameters so gl.check.lambda will work
    lambdas <-
        .gl.parameter.tidy(lambda1, lambda2, lambda3, lambda4, param, lambda5)

    # Check the parameters
    if(!.gl.check.lambda(lambdas, param=param, vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
            "\ndo not produce a proper distribution with the", param,
            "parameterisation - see \ndocumentation for .gl.check.lambda"))
    }

    # calculate u=F(x) numerically, then use qdgl
    # Unless x is outside the range, then density should be zero
    extreme <- .qgld(c(0, 1), lambda1=lambdas, param=param)

    # It may be better to change this to simply
    # (x <= extreme[2])*(x >= extreme[1])
    outside.range <- !as.logical((x<=extreme[2])*(x>=extreme[1]))
    u <- .pgld(x, lambdas, param=param, inverse.eps=inverse.eps,
        max.iterations=max.iterations)
    dens <- .qdgld(u, lambda1=lambdas, param=param)
    dens[outside.range] <- 0

    dens
}


################################################################################


.pgld <-
function(q, lambda1=0, lambda2=NULL, lambda3=NULL, lambda4=NULL,
    param="fmkl", lambda5=NULL, inverse.eps=1e-8, max.iterations=500)
{
    # Thanks to Steve Su, <s.su@qut.edu.au>, for improvements to this code
    # If lambda1 is a vector, the default value for lambda2 will cause a
    # problem.
    # I did have a warning about this, but it will occur too often to make
    # up for the benefit, so I've deleted it.
    # Tidy the parameters so gl.check.lambda will work

    lambdas <-
        .gl.parameter.tidy(lambda1, lambda2, lambda3, lambda4, param, lambda5)

    # Check the parameters
    if(!.gl.check.lambda(lambdas, param=param, vect=TRUE)) {
        stop(paste("The parameter values", lambda1, lambda2, lambda3, lambda4,
            "\ndo not produce a proper distribution with the", param,
            "parameterisation - see \ndocumentation for gl.check.lambda"))
    }
    jr <- q; jr[sort.list(q)] <- seq(along=q)
    order.x<-order(q)
    xx<-sort(q)

    # xx has the sorted data, and jr & order.x the information to get back
    # to theoriginal order.
    extreme <- .qgld(c(inverse.eps, 1-inverse.eps), lambda1=lambdas,
        param=param)
    max.e<-extreme[2]
    min.e<-extreme[1]
    ind.min<-xx<=min.e
    ind.max<-xx>=max.e

    # This simpler comparison works here because we are using inverse.eps as our
    # tolerance
    q<-xx[as.logical((xx<max.e)*(xx>min.e))]

    # We only want to calculate the probabilities for q values inside
    # the support
    length.of.vector <- length(q)

    # Need a blank u to send to C
    u <- 0*q
    result <- switch(param,
        freimer=, # allows for alternate expressions
        frm=, # allows for alternate expressions
        FMKL=, # Notes on .C call - the "numerics", lambdas and inverse.eps
               # don't need the as.??? call as they are implicitly double
        fmkl=.C(
            "gl_fmkl_distfunc", lambdas[1], lambdas[2], lambdas[3], lambdas[4],
            as.double(0), as.double(1), inverse.eps,
            as.integer(max.iterations), as.double(q), as.double(u),
            as.integer(length.of.vector), PACKAGE="fBasics"),
            ramberg=, # Ramberg & Schmeiser
            ram=,
            RS=,
            rs=.C(
                "gl_rs_distfunc", lambdas[1], lambdas[2], lambdas[3], lambdas[4],
                as.double(0), as.double(1), inverse.eps,
            as.integer(max.iterations), as.double(q), as.double(u),
            as.integer(length.of.vector),
            PACKAGE="fBasics"),
        fm5=.C("gl_fm5_distfunc", lambdas[1], lambdas[2], lambdas[3],
            lambdas[4], lambdas[5],
                as.double(0), as.double(1), inverse.eps,
            as.integer(max.iterations), as.double(q), as.double(u),
            as.integer(length.of.vector),
            PACKAGE="fBasics"),
            stop("Error: Parameterisation must be one of fmkl, fm5 or rs")
            ) # closes "switch"
    if (!(is.numeric(result[[1]]))){
        stop("Values for quantiles outside range. This shouldn't happen")
    }
    if (param=="fm5") {
        u <- result[[11]]
    } else {
        u <- result[[10]]
    }
    xx[as.logical((xx<max.e)*(xx>min.e))]<-u
    xx[ind.min]<-0
    xx[ind.max]<-1

    # Revert to the original order of the dataset:
    xx[jr]
}



################################################################################


.qgld <-
function(p, lambda1, lambda2=NULL, lambda3=NULL, lambda4=NULL,
    param="fmkl", lambda5=NULL)
{
    lambdas <-
        .gl.parameter.tidy(lambda1, lambda2, lambda3, lambda4, param, lambda5)

    # Check the values are OK
    if(!.gl.check.lambda(lambdas, param=param, vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas, collapse=" "),
            "\ndo not produce a proper distribution for the", param,
            "parameterisation \n - see documentation for gl.check.lambda"))
    }
    result <- switch(param,
        freimer=,  # allows for alternate expressions
        frm=,  # allows for alternate expressions
        FMKL=,
        fmkl=.qgl.fmkl(p, lambdas),
        ramberg=, # Ramberg & Schmeiser
        ram=,
        RS=,
        rs=.qgl.rs(p, lambdas),
        fm5 = .qgl.fm5(p, lambdas),
        stop("Error: Parameterisation must be fmkl, fm5 or rs")
        ) # closes "switch"
    result
}


# ------------------------------------------------------------------------------


.qgl.fmkl <-
function(p, lambdas)
{
    # No checking - use qgl if you want that
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]
    p <- as.double(p)

    # abandoned this for the simpler one below
    # outside.range <- !as.logical(((p<1)*(p>0))|
    #   (sapply(p, all.equal, 1)=="TRUE")| (sapply(p, all.equal, 0)=="TRUE"))
    outside.range <- !as.logical((p<=1)*(p>=0))

    # u gets only the probabilities in [0, 1]
    u <- p[!outside.range]

    # If OK, determine special cases
    if (lambda3 == 0) {
        if (lambda4 == 0) { # both log
            quants <- lambda1 + (log(u) - log(1 - u))/lambda2
            }
        else    { # l3 zero, l4 non-zero
            quants <- lambda1 +
                (log(u) - ((1 - u)^lambda4 - 1)/lambda4)/lambda2
            }
    } else     { # lambda3 non-zero
        if (lambda4 == 0) { # non-zero, l4 zero
            quants <- lambda1 +
                ((u^lambda3 - 1)/lambda3 - log(1 - u))/lambda2
        } else    { # both non-zero - use usual formula
            quants <- lambda1 + ( ( u ^ lambda3 - 1)    / lambda3 -
                ( (1-u)^lambda4 - 1) /lambda4 ) / lambda2
        }
    }

    # Now we have the quantiles for p values inside [0, 1], put them in
    # the right place in the result vector
    result <- p*NaN
    result[!outside.range] <- quants

    # The remaining "quantiles" are NaN
    result
}


# ------------------------------------------------------------------------------


.qgl.fm5 <-
function(p, lambdas)
{
    # No parameter value checking. If you want that, use qgl!
    lambda5 = as.double(lambdas[5])
    lambda4 = as.double(lambdas[4])
    lambda3 = as.double(lambdas[3])
    lambda2 = as.double(lambdas[2])
    lambda1 = as.double(lambdas[1])
    p <- as.double(p)

    # abandoned this for the simpler
    # outside.range <- !as.logical(((p<1)*(p>0))|
    #   (sapply(p, all.equal, 1)=="TRUE")| (sapply(p, all.equal, 0)=="TRUE"))
    outside.range <- !as.logical((p<=1)*(p>=0))

    # u gets only the probabilities in [0, 1]
    u <- p[!outside.range]

    # If OK, determine special cases
    if (lambda3 == 0) {
        if (lambda4 == 0) { # both log
            quants <- lambda1 + ((1-lambda5)*log(u) -
                (1+lambda5)*log(1 - u))/lambda2
        } else { # l3 zero, l4 non-zero
            quants <- lambda1 +
                ((1-lambda5)*log(u) -
                (1+lambda5)*((1 - u)^lambda4 - 1)/lambda4)/lambda2
        }
    } else { # lambda3 non-zero
        if (lambda4 == 0) { # non-zero, l4 zero
            quants <- lambda1 +
                ((1-lambda5)*(u^lambda3 - 1)/lambda3 -
                (1+lambda5)*log(1 - u))/lambda2
        } else { # both non-zero - use usual formula
            quants <- lambda1 + ((1-lambda5)* ( u ^ lambda3 - 1) / lambda3
                - (1+lambda5)*( (1-u)^lambda4 - 1) /lambda4 ) / lambda2
        }
    }

    # Now we have the quantiles for p values inside [0, 1], put them in
    # the right place in the result vector
    result <- p*NaN
    result[!outside.range] <- quants

    # The remaining "quantiles" are NaN
    result
}


# ------------------------------------------------------------------------------


.qgl.rs <-
function(p, lambdas)
{
    u <- p
    # No parameter value checking - use qgl!
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]

    # At present, I'm rejecting zero values for l3 and l4, though I think
    # there are limit results, so one functional form.
    quants <- lambda1 + ( u ^ lambda3 - (1-u)^lambda4 ) / lambda2
    quants
}


# ------------------------------------------------------------------------------


.qdgld <-
function(p, lambda1, lambda2=NULL, lambda3=NULL, lambda4=NULL,
    param="fmkl", lambda5=NULL)
{
    # Don't allow characters in lambda5 -
    # common error with parameterisation stuff
    if(is.character(lambda5)) {
        stop(paste("lambda5=", lambda5,
            "It should be a number between -1 and 1"))
    }
    lambdas <-
        .gl.parameter.tidy(lambda1, lambda2, lambda3, lambda4, param, lambda5)

    # Check the values are OK
    if(!.gl.check.lambda(lambdas, param=param, vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas, collapse=" "),
        "\ndo not produce a proper distribution for the", param,
        "parameterisation \n - see documentation for gl.check.lambda"))
    }
    result <- switch(param,

    # Different tests apply for each parameterisation
        freimer=,  # allows for alternate expressions
        frm=,  # allows for alternate expressions
        FMKL=,
        fmkl=.qdgl.fmkl(p, lambdas),
        ramberg=, # Ramberg & Schmeiser
        ram=,
        RS=,
        rs=.qdgl.rs(p, lambdas),
        fm5 = .qdgl.fm5(p, lambdas),
        stop("Error: Parameterisation must be fmkl, fm5 or rs")
        ) # closes "switch"
    result
}


# ------------------------------------------------------------------------------


.qdgl.rs <-
function(p, lambdas)
{
    # Check the values are OK)
    if(!.gl.check.lambda(lambdas, param="rs", vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas, collapse=" "),
        "\ndo not produce a proper distribution with the RS parameterisation",
        " - see \ndocumentation for gl.check.lambda"))
    }
    outside.range <- !as.logical((p<=1)*(p>=0))

    # u gets only the probabilities in [0, 1]
    u <- p[!outside.range]
    dens <-  lambdas[2]/(lambdas[3] * (p^(lambdas[3] -1)) +
        lambdas[4] * ((1 - p)^(lambdas[4] -1)))

    dens
}


# ------------------------------------------------------------------------------


.qdgl.fmkl <-
function(p, lambdas)
{
    # Check the values are OK)
    if(!.gl.check.lambda(lambdas, param="fmkl", vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas, collapse=" "),
        "\ndo not produce a proper distribution with the FMKL",
        "parameterisation - see \ndocumentation for gl.check.lambda"))
    }
    outside.range <- !as.logical((p<=1)*(p>=0))

    # u gets only the probabilities in [0, 1]
    u <- p[!outside.range]

    # The density is given by 1/Q'(u)
    dens <- lambdas[2]/(p^(lambdas[3] - 1) + (1 - p)^(lambdas[4] - 1))

    dens
}


# ------------------------------------------------------------------------------


.qdgl.fm5 <-
function(p, lambdas)
{
    # Check the values are OK)
    if(!.gl.check.lambda(lambdas, param="fm5", vect=TRUE)) {
        stop(paste("The parameter values", paste(lambdas, collapse=" "),
            "\ndo not produce a proper distribution with the FM5",
            "parameterisation - see \ndocumentation for gl.check.lambda"))
    }
    outside.range <- !as.logical((p<=1)*(p>=0))

    # u gets only the probabilities in [0, 1]
    u <- p[!outside.range]

    # The density is given by 1/Q'(u)
    dens <- lambdas[2]/((1-lambdas[5])*(u^(lambdas[3] - 1)) +
        (1+lambdas[5])*((1 - u)^(lambdas[4] - 1)) )

    dens
}


################################################################################


.rgld <-
function(n, lambda1=0, lambda2=NULL, lambda3=NULL, lambda4=NULL,
    param="fmkl", lambda5=NULL)
{
    # Check the parameters
    lambdas <-
        .gl.parameter.tidy(lambda1, lambda2, lambda3, lambda4, param, lambda5)

    # Check the values are OK
    if(!.gl.check.lambda(lambdas, param=param, vect=TRUE)) {
            stop(paste("The parameter values",  lambdas,
            "\ndo not produce a proper distribution for the", param,
            "parameterisation \n - see documentation for gl.check.lambda"))
    }

    # Produce the uniform data
    p <- runif(n)

    # convert to gl
    res <- .qgld(p, lambda1=lambdas, param=param)
    res
}


################################################################################


.gl.parameter.tidy <-
function(lambda1, lambda2=NULL, lambda3=NULL, lambda4=NULL, param="fmkl",
    lambda5=NULL)
{
    # Don't allow characters in lambda5 - common error with parameterisation stuff
    if(is.character(lambda5)) {stop(paste("lambda5=", lambda5,
        "It should be a number between -1 and 1"))}

    # Don't allow numbers in parameterisation -
    #   included as a warning here, so the main one is a stop.
    if(!is.character(param)) {warning(paste("param=", param,
        "It shouldn't be a number, ",
        "it should be a string describing the parameterisation"))
    }
    if(length(lambda1) > 1) {
        # using a vector for the parameters.
        # Check that there aren't values in the individual lambda arguments
        if (!(is.null(lambda2) & is.null(lambda3)& is.null(lambda4) &
            is.null(lambda5)) )
            { stop("Call includes vector version of the lambda parameters",
            "as well as the \nscalar version")
        }
        if ((length(lambda1) < 4) | (length(lambda1) > 5 ) )
            { stop(paste("argument lambda1 has length", length(lambda1),
            "\nThis should be 1 (lambda parameters given as seperate ",
            "arguments), 4 (vector argument \n for RS or FMKL ",
            "parameterisation) or 5 (vector argument for fm5",
            "parameterisation"))
        }
        if (length(lambda1)== 5)
            { if (param != "fm5") {
                stop(paste("argument lambda1 has length", length(lambda1),
                "which is not valid for the", param, "\nparameterisation"))
            }
            # else --- fm5, in vector form, ready for gl.check.lambda
        }
        if (length(lambda1)== 4) {
            if (param == "fm5" )
                { stop(paste("argument lambda1 has length 4, which is not",
                " valid for the fm5 \nparameterisation"))
            }
        # else --- 4 parameter versions in vector form, ready for gl.check.lambda
        }
    } else {
        # single parameter arguments - check they are there,
        # then collect them together
        if (is.null(lambda2)) { stop("No value for lambda2") }
        if (is.null(lambda3)) { stop("No value for lambda3") }
        if (is.null(lambda4)) { stop("No value for lambda4") }
        if ((is.null(lambda5)) & param=="fm5" ) {
            stop("No value for lambda5") }
        if (!(is.null(lambda5)) & param!="fm5") {
            stop(paste("lambda5=", lambda5, " but there is no lambda 5 for the\n",
            param, "parameterisation")) }
        if (param != "fm5") { # A 4 parameter version
            lambda1 <- c(lambda1, lambda2, lambda3, lambda4)
        } else { # fm5
            lambda1 <- c(lambda1, lambda2, lambda3, lambda4, lambda5)
        }
    }

    # There is now an error if there is the wrong number of parameters, and
    # lambda1 returned as a vector with 4 or 5 elements
    # as.double is needed to remove data.frame attributes if lambda1 was
    # extracted from a data.frame

    as.double(lambda1)
}


# ------------------------------------------------------------------------------


.gl.check.lambda <-
function(lambdas, lambda2=NULL, lambda3=NULL, lambda4=NULL,
    param="fmkl", lambda5=NULL, vect=FALSE)
{
    # Checks to see that the lambda values given are allowed.
    # There is a function called .gl.parameter.tidy that does the tidying
    # around of parameters.  It return a single vector, which contains the
    # parameters.

    # If you call this after .gl.parameter.tidy, let it know with the vect=TRUE
    # argument

    # If vect=TRUE, we don't need to tidy

    if (vect) {
        if (!is.null(lambda3)) {
            warning("lambda3 should be null because you claim the",
            " parameters are in a vector")
        }
    } else {
        lambdas <-
            .gl.parameter.tidy(lambdas, lambda2, lambda3, lambda4, param, lambda5)
    }
    if(param=="fm5"){lambda5 = lambdas[5]}
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]

    # I did have a check for finite lambdas, but that caused a
    # problem with data frames,
    # so I removed it - still need to include the limit results
    param <- switch(param,
        # Different tests apply for each parameterisation
        freimer=,  # allows for alternate expressions
        frm=,  # allows for alternate expressions
        FMKL=,
        fmkl={
        if (lambda2<=0) {return(FALSE)}
        else {return(TRUE)}
        },
        ramberg=, # Ramberg & Schmeiser
        ram=,
        RS=,
        rs={
        if (lambda3*lambda4>0) { # regions 3 and 4
            # all values of lambda 3 and lambda 4 OK
            # check lambda 2
            if ((lambda3>0)&(lambda4>0)) { # region 3 - l2 >0
                if (lambda2<=0) {
                    ret <- FALSE
                } else {
                    ret <- TRUE
                }
            }
            if ((lambda3<0)&(lambda4<0)) { # region 4 - l2 <0
                if (lambda2>=0) {
                    ret <- FALSE
                } else {
                    ret <- TRUE
                }
            }
        } else { # other quadrants - lambda 2 must be negative, and lambda3
            # lambda 4 values need checking.
            if (lambda2>=0) {return(FALSE)}
            # Rectangular regions where RS is not defined
            if ((lambda3>0)&(lambda3<1)&(lambda4<0)) {return(FALSE)}
            if ((lambda4>0)&(lambda4<1)&(lambda3<0)) {return(FALSE)}
            # Different here because there are a
            # number of ways in which the parameters can fail.
            #
            # Curved regions where RS is not defined
            # change to shorter var names
            lc <- lambda3
            ld <- lambda4
            if ((lambda3>-1)&(lambda3<0)&(lambda4>1)) {  # region 5 or not?
                if ( ((1-lc)^(1-lc)*(ld-1)^(ld-1))/((ld-lc)^(ld-lc)) > -lc/ld )
                    {return(FALSE)}
                else     {return(TRUE)}
                }
            # Second curved region
            if ((lambda4>-1)&(lambda4<0)&(lambda3>1)) {  # region 6 or not?
                if ( ((1-ld)^(1-ld)*(lc-1)^(lc-1))/((lc-ld)^(lc-ld)) > -ld/lc )
                    {return(FALSE)}
                else     {return(TRUE)}
                }
            # There may be some limit results that mean these are not correct, but
            # I'll check that later
            # This is not the place where the possible l3, l4 zero values should appear
            if (lambda3 == 0) {
                warning('lambda 3 = 0 with RS parameterisation - possible problem')
                if (lambda4 == 0) {return(FALSE)}
                else {return(TRUE)}
                }
            if (lambda4 == 0) {
                warning('lambda 5 = 0 with RS parameterisation - possible problem')
                if (lambda4 == 0) {return(FALSE)}
                else {return(TRUE)}
                }
            # If we get here, then the parameters are OK.
            ret <- TRUE
            }
        },
        fm5={
            # make lambda5 - in here so it doesn't stuff up the
            # other parameterisations
            lambda5 <- lambdas[5]
            if (lambda2<=0) {ret <- FALSE}
            else { #Check lambda5 - should be between -1 and 1, but
                # I haven't checked this against a piece of paper
                if ((lambda5 >= -1) & (lambda5 <= 1)) {
                    ret <- TRUE
                } else {
                    ret <- FALSE
                }
            }
        },
        stop("Error when checking validity of parameters.\n",
            " Parameterisation must be fmkl, rs or fm5")
        ) # closes "switch"
    ret
}


################################################################################

