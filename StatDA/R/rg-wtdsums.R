"rg.wtdsums" <-
function(x, ri, xcentr = NULL, xdisp = NULL)
{
        # Function to compute weighted sums for a matrix based on computed
        # quantiles and user-defined relative importances, see Garrett, RG &
        # Grunsky, EC, 2001: Weighted Sums - knowledge based empirical indices
        # for use in exploration geochemistry, Geochem: Explor, Env, Anal, 1(4):135-141.
        # Optionally estimates of centres (mean) and dispersion (sd) may be provided,
        # e.g., with non-robust or other robust estimates.
        #
        # Example using the Howarth and Sinding-Larsen Norwegian stream sediment
        # data, and with sind.dat attached:
        #
        # rg.wtdsums(cbind(Zn,Cd,Fe,Mn),c(2,1,-1,-1))
        #
        ncolx <- length(x[1,  ])
        nleni <- length(ri)
        if(ncolx != nleni)
                stop("\n  Number of variables and importances do not match")
        w <- ri/sum(abs(ri))
        a <- w/sqrt(sum(w * w))
        if(is.null(xcentr) & is.null(xdisp)) {
                xsumm <- matrix(nrow = ncolx, ncol = 3)
                xsumm <- apply(x, 2, quantile, c(0.25, 0.5, 0.75), na.rm = TRUE)
                xcentr <- xsumm[2,  ]
                xdisp <- 0.74129999999999996 * (xsumm[3,  ] - xsumm[1,  ])
        }
        else {
                if(length(xcentr) != nleni)
                        stop("\n  Numbers of variables and centre measures do not match")
                if(length(xdisp) != nleni)
                        stop("\n  Numbers of variables and dispersions do not match")
        }
        xx <- sweep(x, 2, xcentr, "-")
        z <- sweep(xx, 2, xdisp, "/")
        ws <- as.matrix(z) %*% a
        invisible(list(input = deparse(substitute(x)), centr = xcentr, disp = xdisp, ri = ri,
                w = w, a = a, ws = ws))
}

