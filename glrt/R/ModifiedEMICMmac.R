ModifiedEMICMmac <-
function (A, EMstep = TRUE, ICMstep = TRUE, keepiter = FALSE, tol = 1e-06, tolbis = 1e-07, maxiter = 1000)

{

    if (!EMstep && !ICMstep)
    {

        print("One of EMstep or ICMstep must be true.")

        return(NULL)

    }

    Meps <- .Machine$double.eps

    m <- dim(A)[1]

    n <- dim(A)[2]

    tA <- t(A)

    if (m == 1)
    {

        ret <- NULL

        ret$sigma <- 1

        ret$weights <- n

        ret$lastchange <- 0

        ret$numiter <- 0

        return(ret)

    }

    WW <- matrix(0, m, n)

    for (i in 1:(m - 1))
 WW[i, ] <- A[i, ] - A[i + 1, ]

    WW[m, ] <- A[m, ]

    sigma <- cumsum(apply(A, 1, sum)/sum(A))

    numiter <- 0

    oldsigma <- rep(-1, m)

    if (keepiter)
 
        iter <- sigma

    while (max(abs(oldsigma - sigma)) > tol && numiter <= maxiter)
    {

        oldsigma <- sigma

        if (EMstep) {

            pvec <- diff(c(0, sigma))

            temp <- sweep(A, 1, pvec, FUN = "*")

            if (sum(apply(temp, 2, sum) == 0) == 0) {

                pvec <- apply(sweep(temp, 2, apply(temp, 2, sum), 
FUN = "/"), 1, sum)

                sigma <- cumsum(pvec)/sum(pvec)

            }

            if (keepiter)
 
                iter <- rbind(iter, sigma)

        }

        if (ICMstep) {

            Wps <- 1/(t(WW) %*% sigma)

            weights <- (abs(WW) %*% Wps^2)

            increment <- as.vector((WW %*% Wps)/weights)

            sigma <- sigma + increment

            sigma[m] <- 1

            if (keepiter)
 
                iter <- rbind(iter, sigma)

            temp <- PMGA(sigma[-m], weights[-m])

            poolnum <- c(0, cumsum(temp$poolnum))

            for (i in 2:length(poolnum))
 for (j in (poolnum[i - 1] + 1):poolnum[i])
 sigma[j] <- temp$est[i - 
1]

            if (keepiter)
 
                iter <- rbind(iter, sigma)

            temp <- c(0, sigma)

            pvec <- diff(c(0, oldsigma))

            ndir <- diff(c(0, temp[2:(m + 1)])) - pvec

            pvec <- Bisect(tA, pvec, ndir, Meps, tolbis)

            sigma <- cumsum(pvec)

            if (keepiter)
 
                iter <- rbind(iter, sigma)

        }

        numiter <- numiter + 1

    }

    if (numiter == maxiter)
 
        warning("EM/ICM may have failed to converge.")

    pf <- diff(c(0, sigma))

    if(ICMstep)

        ret <- list(sigma = sigma, pf = pf, llk = sum(log(t(A) %*% pf)), weights = as.vector(weights), lastchange = sigma - oldsigma, numiter = numiter, eps = tol)

    else

        ret <- list(sigma = sigma, pf = pf, llk = sum(log(t(A) %*% pf)), lastchange = sigma - oldsigma, numiter = numiter, eps = tol)

    if (keepiter) {

        if (EMstep && ICMstep)
 
            dimnames(iter) <- list(c("Seed", rep(c("EM", "Fisher", "PMGA", "Bisect"), numiter)), NULL)

        else if (EMstep)
 
            dimnames(iter) <- list(rep("EM", numiter + 1), NULL)

        else dimnames(iter) <- list(c("Seed", rep(c("Fisher", 
"PMGA"), numiter)), NULL)

        ret$iter <- iter

    }

    ret

}

