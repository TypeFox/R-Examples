
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTIONS:            FRACTIONAL BROWNIAN MOTION:
#  fbmSim                Generates fractional Brownian motion
#  .fbmSim.mvn            Numerical approximation of the stochastic integral
#  .fbmSim.chol           Choleki's decomposition of the covariance matrix
#  .fbmSim.lev            Method of Levinson
#  .fbmSim.circ           method of Wood and Chan
#  .fbmSim.wave           Wavelet synthesis
#  .convol                Internal Convolution
#  .fbmSlider             Displays fbm simulated time Series
################################################################################


fbmSim =
function(n = 100, H = 0.7, method = c("mvn", "chol", "lev", "circ", "wave"),
waveJ = 7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulation of fractional Brownian motion by five different methods

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion

    # Details:
    #   The underlying functions were ported from SPlus code written
    #   by J.F. Couerjolly. They are documented in the reference given
    #   below.

    # Reference:
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian
    #       Motion: A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Initialization:
    method = match.arg(method)

    # Match Function:
    fun = paste(".fbmSim.", method, sep = "")
    funFBM = match.fun(fun)

    # Simulate:
    if (method == "wave") {
        ans = funFBM(n, H, waveJ, doplot, fgn)
    } else {
        ans = funFBM(n, H, doplot, fgn)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.fbmSim.mvn =
function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n }
    #   by numerical approximation of stochastic integral

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Mandelbrot B. and Van Ness,
    #       Fractional brownian motions,
    #       fractional noises and applications, SIAM Review, 10, n.4, 1968.
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Initialization:
    dB1 = rnorm(n)
    borne = trunc(n^1.5)
    dB2 = rnorm(borne)
    fBm = rep(0, n)
    CH = sqrt(gamma(2 * H + 1) * sin(pi * H))/gamma(H + 1/2)    ##

    # Main Program:
    ind1 = (1:n)^(H - 1/2)
    ind2 = (1:(borne + n))^(H - 1/2)
    for(i in (2:n)) {
        I1 = dB1[(i - 1):1] * ind1[1:(i - 1)]
        I2 = (ind2[(i + 1):(i + borne)] - ind2[1:borne]) * dB2
        fBm[i] = sum(I1) + sum(I2)
    }
    fBm = fBm * n^( - H) * CH
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "mvnFBM Path"
    if (fgn) {
        ans = c(fBm[1], diff(fBm))
        Title = "mvnFGN Path"
    }

    # Plot of fBm
    if (doplot) {
        time = 1:n
        Nchar = as.character(n)
        Nleg = paste("N=", Nchar, sep = "")
        Hchar = as.character(round(H, 3))
        Hleg = paste(", H=", Hchar, sep = "")
        NHleg = paste(c(Nleg, Hleg), collapse = "")
        leg = paste(c(Title, NHleg), collapse = " - ")
        plot(time, ans, type = "l", main = leg, col = "steelblue")
        grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "mvn", H = H)
    ans
}


# ------------------------------------------------------------------------------


.fbmSim.wave =
function(n = 1000, H = 0.7, J = 7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   J : resolution
    #   doplot : = TRUE ----> plot of path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by wavelet synthesis

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Abry P. and Sellan F.,
    #       The wavelet-based synthesis
    #       for fractional Brownian motion, Applied and computational
    #       harmonic analysis, 1996 + Matlab scripts from P. Abry
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Daubechies filter of length 20
    Db20 = c(0.026670057901000001, 0.188176800078)
    Db20 = c(Db20, 0.52720118932000004, 0.688459039454)
    Db20 = c(Db20, 0.28117234366100002, -0.24984642432699999)
    Db20 = c(Db20, -0.19594627437699999, 0.127369340336)
    Db20 = c(Db20, 0.093057364604000006, -0.071394147165999997)
    Db20 = c(Db20, -0.029457536821999999, 0.033212674058999997)
    Db20 = c(Db20, 0.0036065535670000001, -0.010733175483)
    Db20 = c(Db20, 0.001395351747, 0.0019924052950000002)
    Db20 = c(Db20, -0.00068585669500000003, -0.000116466855)
    Db20 = c(Db20, 9.3588670000000005e-05, -1.3264203000000001e-05)
    secu = 2 * length(Db20)

    # Quadrature mirror filters of Db20
    Db20qmf = (-1)^(0:19) * Db20
    Db20qmf = Db20qmf[20:1]
    nqmf = -18

    # Truncated fractional coefficients appearing in fractional integration
    # of the multiresolution analysis
    prec = 0.0060000000000000001
    hmoy = c(1, 1)
    s = H + 1/2
    d = H - 1/2
    if (H == 1/2) {
        ckbeta = c(1, 0)
        ckalpha = c(1, 0)
    } else {
        # Truncature at order prec
        ckalpha = 1
        ckbeta = 1
        ka = 2
        kb = 2
        while(abs(ckalpha[ka - 1]) > prec) {
            g = gamma(1 + d)/gamma(ka)/gamma(d + 2 - ka)
            if (is.na(g))
                g = 0
            ckalpha = c(ckalpha, g)
            ka = ka + 1
        }
        while(abs(ckbeta[kb - 1]) > prec) {
            g = gamma(kb - 1 + d)/gamma(kb)/gamma(d)
            if (is.na(g))
                g = 0
            ckbeta = c(ckbeta, g)
            kb = kb + 1
        }
    }
    lckbeta = length(ckbeta)
    lckalpha = length(ckalpha)  ##

    # Number of starting points
    nbmax = max(length(ckbeta), length(ckalpha))
    nb0 = n/(2^(J)) + 2 * secu  ##

    # Sequence fs1:
    fs1 = .convol(ckalpha, Db20)
    fs1 = .convol(fs1, hmoy)
    fs1 = 2^( - s) * fs1
    fs1 = fs1 * sqrt(2) #

    # Sequence gs1:
    gs12 = .convol(ckbeta, Db20qmf)
    gs1 = cumsum(gs12)
    gs1 = 2^(s) * gs1
    gs1 = gs1 * sqrt(2) ##

    # Initialization:
    nb1 = nb0 + nbmax
    bb = rnorm(nb1)
    b1 = .convol(bb, ckbeta)
    bh = cumsum(b1)
    bh = bh[c(nbmax:(nb0 + nbmax - 1))]
    appro = bh
    tappro = length(appro)  ##

    # Useful function:
    dilatation = function(vect) {
        # dilates one time vector vect
        ldil = 2 * length(vect) - 1
        dil = rep(0, ldil)
        dil[seq(1, ldil, by = 2)] = vect
        drop(dil)
    }

    # Synthese's algorithm:
    for(j in 0:(J - 1)) {
        appro = dilatation(appro)
        appro = .convol(appro, fs1)
        appro = appro[1:(2 * tappro)]
        detail = rnorm(tappro) * 2^(j/2) * 4^( - s) * 2^( - j * s)
        detail = dilatation(detail)
        detail = .convol(detail, gs1)
        detail = detail[( - nqmf + 1):( - nqmf + 2 * tappro)]
        appro = appro + detail
        tappro = length(appro)
    }
    debut = (tappro - n)/2
    fBm = appro[c((debut + 1):(debut + n))]
    fBm = fBm - fBm[1]
    fGn = c(fBm[1], diff(fBm))
    fGn = fGn * 2^(J * H) * n^( - H)    # path on [0,1]
    fBm = cumsum(fGn)
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "waveFBM Path"
    if (fgn) {
        ans = c(fBm[1], diff(fBm))
        Title = "waveFGN Path"
    }

    # Plot of fBM/FGN:
    if (doplot) {
        time = 1:n
        Nchar = as.character(n)
        Nleg = paste("N=", Nchar, sep = "")
        Hchar = as.character(round(H, 3))
        Hleg = paste(", H=", Hchar, sep = "")
        NHleg = paste(c(Nleg, Hleg), collapse = "")
        leg = paste(c(Title, NHleg), collapse = " - ")
        plot(time, ans, type = "l", main = leg, col = "steelblue")
        grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "wave", H = H)
    ans
}


# ------------------------------------------------------------------------------


.convol =
function(x, y)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x,y : vectors

    # Value:
    #   convolution of vectors x and y

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # FUNCTION:

    # Convolution:
    if (missing(x) | missing(y)) {
        break
    } else {
        a = c(x, rep(0, (length(y) - 1)))
        b = c(y, rep(0, (length(x) - 1)))
        a = fft(a, inverse = FALSE)
        b = fft(b, inverse = FALSE)
        conv = a * b
        conv = Re(fft(conv, inverse = TRUE))
        conv = conv/length(conv)
        drop(conv)
    }
}


# ------------------------------------------------------------------------------


.fbmSim.chol =
function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ----> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n }
    #   by Choleki's decomposition of the covariance matrix of the fBm

    # Author:
    #   Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Construction of covariance matrix of fBm
    H2 = 2 * H
    matcov = matrix(0, n - 1, n - 1)
    for(i in (1:(n - 1))) {
        j = i:(n - 1)
        r = 0.5 * (abs(i)^H2 + abs(j)^H2 - abs(j - i)^H2)
        r = r/n^H2
        matcov[i, j] = r
        matcov[j, i] = matcov[i, j]
    }
    L = chol(matcov)
    Z = rnorm(n - 1)
    fBm = t(L) %*% Z
    fBm = c(0, fBm)

    # Result:
    ans = drop(fBm)
    Title = "cholFBM Path"
    if (fgn) {
        ans = c(fBm[1], diff(fBm))
        Title = "cholFGN Path"
    }

    # Plot of fBm:
    if (doplot) {
        time = 1:n
        Nchar = as.character(n)
        Nleg = paste("N=", Nchar, sep = "")
        Hchar = as.character(round(H, 3))
        Hleg = paste(", H=", Hchar, sep = "")
        NHleg = paste(c(Nleg, Hleg), collapse = "")
        leg = paste(c(Title, NHleg), collapse = " - ")
        plot(time, ans, type = "l", main = leg, col = "steelblue")
        grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "chol", H = H)
    ans
}


# ------------------------------------------------------------------------------


.fbmSim.lev =
function(n = 1000, H = 0.7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   plotfBm :  =1 ---> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by Levinson's method

    # Author:
    #   Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Peltier R.F.,
    #       Processus stochastiques fractals avec
    #       applications en finance, these de doctorat, p.42, 28.12.1997
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # Covariances of fGn:
    k = 0:(n - 1)
    H2 = 2 * H
    r = (abs((k - 1)/n)^H2 - 2 * (k/n)^H2 + ((k + 1)/n)^H2)/2

    # Initialization of algorithm:
    y = rnorm(n)
    fGn = rep(0, n)
    v1 = r
    v2 = c(0, r[c(2:n)], 0)
    k =  - v2[2]
    aa = sqrt(r[1]) #

    # Levinson's algorithm:
    for(j in (2:n)) {
        aa = aa * sqrt(1 - k * k)
        v = k * v2[c(j:n)] + v1[c((j - 1):(n - 1))]
        v2[c(j:n)] = v2[c(j:n)] + k * v1[c((j - 1):(n - 1))]
        v1[c(j:n)] = v
        bb = y[j]/aa
        fGn[c(j:n)] = fGn[c(j:n)] + bb * v1[c(j:n)]
        k =  - v2[j + 1]/(aa * aa)
    }
    fBm = cumsum(fGn)
    fBm[1] = 0

    # Result:
    ans = drop(fBm)
    Title = "levFBM Path"
    if (fgn) {
        ans = c(fBm[1], diff(fBm))
        Title = "levFGN Path"
    }

    # Plot of fBm:
    if (doplot) {
        time = 1:n
        Nchar = as.character(n)
        Nleg = paste("N=", Nchar, sep = "")
        Hchar = as.character(round(H, 3))
        Hleg = paste(", H=", Hchar, sep = "")
        NHleg = paste(c(Nleg, Hleg), collapse = "")
        leg = paste(c(Title, NHleg), collapse = " - ")
        plot(time, ans, type = "l", main = leg, col = "steelblue")
        grid()
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "lev", H = H)
    ans
}


# ------------------------------------------------------------------------------


.fbmSim.circ =
function(n = 100, H = 0.7, doplot = TRUE, fgn = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n : length of the desired sample
    #   H : self-similarity parameter
    #   doplot : = TRUE ---> plot path of fBm

    # Value:
    #   Simulation of a standard fractional Brownian motion
    #   at times { 0, 1/n,..., n-1/n } by Wood-Chan's method

    # Author:
    #    Coeurjolly 06/2000 of the original SPlus port

    # Reference:
    #   Wood A. and Chan G.,
    #       Simulation of stationnary Gaussian processes,
    #       Journal of computational and grahical statistics, Vol.3, 1994.
    #   Couerjolly J.F.,
    #       Simulation and Identification of the Fractional Brownian motion:
    #       A Bibliographical and Comparative Study,
    #       Journal of Statistical Software 5, 2000

    # FUNCTION:

    # First line of the circulant matrix, C, built via covariances of fGn
    lineC =
    function(n, H, m)
    {
        k = 0:(m - 1)
        H2 = 2 * H
        v = (abs((k - 1)/n)^H2 - 2 * (k/n)^H2 + ((k + 1)/n)^H2)/2
        ind = c(0:(m/2 - 1), m/2, (m/2 - 1):1)
        v = v[ind + 1]
        drop(v)
    }

    # Next power of two > n:
    m = 2
    repeat {
        m = 2 * m
        if (m >= (n - 1)) break
    }
    stockm = m  ##

    # Research of the power of two (<2^18) such that C is definite positive:
    repeat {
        m = 2 * m
        eigenvalC = lineC(n, H, m)
        eigenvalC = fft(c(eigenvalC), inverse = FALSE)
        ### DW: That doesn't work on complex vectors !
        ### if ((all(eigenvalC > 0)) | (m > 2^17)) break
        ### We use:
        if ((all(Re(eigenvalC) > 0)) | (m > 2^17)) break
    }
    if (m > 2^17) {
        cat("----> exact method, impossible!!", fill = TRUE)
        cat("----> can't find m such that C is definite positive", fill = TRUE)
        break
    } else {
        # Simulation of W=(Q)^t Z, where Z leads N(0,I_m)
        # and   (Q)_{jk} = m^(-1/2) exp(-2i pi jk/m):
        ar = rnorm(m/2 + 1)
        ai = rnorm(m/2 + 1)
        ar[1] = sqrt(2) * ar[1]
        ar[(m/2 + 1)] = sqrt(2) * ar[(m/2 + 1)]
        ai[1] = 0
        ai[(m/2 + 1)] = 0
        ar = c(ar[c(1:(m/2 + 1))], ar[c((m/2):2)])
        aic =  - ai
        ai = c(ai[c(1:(m/2 + 1))], aic[c((m/2):2)])
        W = complex(real = ar, imaginary = ai)  ##

        # Reconstruction of the fGn:
        W = (sqrt(eigenvalC)) * W
        fGn = fft(W, inverse = FALSE)
        fGn = (1/(sqrt(2 * m))) * fGn
        fGn = Re(fGn[c(1:n)])
        fBm = cumsum(fGn)
        fBm[1] = 0

        # Result:
        ans = drop(fBm)
        Title = "circFBM Path"
        if (fgn) {
            ans = c(fBm[1], diff(fBm))
            Title = "circFGN Path"
        }

        # Plot of fBm:
        if (doplot) {
            time = 1:n
            Nchar = as.character(n)
            Nleg = paste("N=", Nchar, sep = "")
            Hchar = as.character(round(H, 3))
            Hleg = paste(", H=", Hchar, sep = "")
            NHleg = paste(c(Nleg, Hleg), collapse = "")
            leg = paste(c(Title, NHleg), collapse = " - ")
            plot(time, ans, type = "l", main = leg, col = "steelblue")
            grid()
        }
    }

    # Return Value:
    ans = as.ts(ans)
    attr(ans, "control") <- c(method = "circ", H = H)
    ans
}


# ------------------------------------------------------------------------------


.fbmSlider =
function()
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Displays fbm simulated time Series

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        n      = fBasics:::.sliderMenu(no = 1)
        H      = fBasics:::.sliderMenu(no = 2)
        method = fBasics:::.sliderMenu(no = 3)

        # Select Method:
        Method = c("mvn", "chol", "lev", "circ", "wave")
        Method = Method[method]

        # Frame:
        par(mfrow = c(2, 1))

        # FBM TimeSeries:
        fbmSim(n = n, H = H, method = Method, doplot = TRUE)

        # FGN TimeSeries:
        fbmSim(n = n, H = H, method = Method, doplot = TRUE, fgn = TRUE)

        # Reset Frame:
        par(mfrow = c(1, 1))
    }

    # Open Slider Menu:
    fBasics:::.sliderMenu(refresh.code,
       names =       c(  "n",    "H", "method"),
       minima =      c(   10,   0.01,       1),
       maxima =      c(  200,   0.99,       5),
       resolutions = c(   10,   0.01,       1),
       starts =      c(  100,   0.70,       1))
}


################################################################################

