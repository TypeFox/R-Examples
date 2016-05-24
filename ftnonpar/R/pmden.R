"pmden" <-
function (x, DISCR = FALSE, verbose = FALSE, bandwidth = -1, 
    extrema.nr = -1, accuracy = mad(x)/1000, extrema.mean = TRUE, 
    maxkuipnr = 19, asympbounds = FALSE, tolerance = 1e-08, localsq = TRUE, 
    locsq.factor = 0.9) 
{
    nsamp <- length(x)
    if (asympbounds || nsamp > max(kuipdiffbounds.x)) 
        currbounds <- kuipdiffbounds[length(kuipdiffbounds.x), 
            ] * sqrt(max(kuipdiffbounds.x))/sqrt(nsamp)
    else {
        currbounds <- double(maxkuipnr)
        for (i in 1:maxkuipnr) currbounds[i] <- approx(kuipdiffbounds.x, 
            kuipdiffbounds[, i], nsamp, rule = 2)$y
    }
    if (maxkuipnr > dim(kuipdiffbounds)[2]) 
        stop("maxkuipnr is too large")
    if (DISCR) {
        datax <- as.double(levels(as.factor(x)))
        N <- length(datax)
        dataemp <- c(0, as.double(summary(as.factor(x), maxsum = N))/nsamp)
        fdist.y <- cumsum(dataemp)
        x <- 0:N
        nsamp <- N + 1
    }
    else {
        dataemp <- c(0, rep(1/(nsamp - 1), nsamp - 1))
        x <- sort(x)
        if (min((x[-1] - x[-nsamp])/(x[nsamp] - x[1])) < 1e-14) {
            newx <- x / accuracy
            j <- 1
            while(j < nsamp)
              {
              ind <- (j:nsamp)[newx[j:nsamp]<floor(newx[j])+1]
              if(length(ind) == 0) # this can happen if newx[j] is extremely large
                k <- j
              else
                k <- max((j:nsamp)[newx[j:nsamp]<floor(newx[j])+1])

              newx[j:k] <- (1:(k-j+1))/(k-j+2)+floor(newx[j])
              j <- k+1
              }
            x <- newx * accuracy
        }
        fdist.y <- c(seq(0, 1, length = nsamp))
    }
    if (bandwidth > 0) 
        eps <- rep(bandwidth, nsamp)
    else {
        currprecision <- 0.5
        eps <- rep(0.5, nsamp)
    }
    eps[1] <- 0
    eps[nsamp] <- 0
    repeat {
        lower <- fdist.y - eps
        upper <- fdist.y + eps
        fts <- tautstring(x, fdist.y, lower, upper, upper[1], 
            lower[length(lower)], extrmean = extrema.mean)
        x.string <- fts$string
        if (sum(x.string[-1] != x.string[-(nsamp - 1)]) > 0) {
            ind1 <- min((1:(nsamp - 2))[x.string[-1] != x.string[-(nsamp - 
                1)]])
            ind2 <- max((1:(nsamp - 2))[x.string[-1] != x.string[-(nsamp - 
                1)]])
            if (x.string[ind1] > x.string[ind1 + 1]) 
                fts$nmax <- fts$nmax + 1
            if (x.string[ind2] < x.string[ind2 + 1]) 
                fts$nmax <- fts$nmax + 1
        }
        lastunif <- approx(fts$knotst, fts$knotsy, x)$y
        if (verbose) {
            par(mfrow = c(2, 1))
            if (DISCR) {
                plot(datax, dataemp[-1], col = "grey")
                lines(datax, x.string, col = "red")
            }
            else {
                hist(x, 40, prob = TRUE)
                lines((rep(x, rep(2, length(x))))[-c(1, 2 * length(x))], 
                  rep(x.string, rep(2, length(x.string))), col = "red")
            }
            plot(x, upper, type = "l")
            lines(fts$knotst, fts$knotsy, col = "red")
            lines(fts$knotst, fdist.y[fts$knotsind], col = "green")
            lines(x, lower)
        }
        if (bandwidth > 0) 
            break
        if (extrema.nr > 0) {
            if (fts$nmax > extrema.nr) 
                eps <- eps + currprecision
            if (currprecision < tolerance) {
                if (fts$nmax <= extrema.nr) 
                  break
            }
            else {
                currprecision <- currprecision/2
                eps[eps > 0] <- eps[eps > 0] - currprecision
            }
        }
        else {
            diff <- cumsum(dataemp) - lastunif
            currkkuip <- kkuip(diff, maxkuipnr)$met
            kuipinds <- c(currkkuip[1], currkkuip[-1] - currkkuip[-maxkuipnr]) > 
                currbounds + 1e-08
            if (sum(kuipinds) != 0) 
                eps[eps > 0] <- (currbounds[kuipinds])[1]/2
            else if (localsq) {
                irmax <- floor(log(nsamp)) + 3
                icomax <- 1
                prp <- exp(-1)
                currsum <- 2 * exp(-1)
                while (log(currsum) < log(0.95)/nsamp) {
                  icomax <- icomax + 1
                  prp <- prp/icomax
                  currsum <- currsum + prp
                }
                kni <- rep(0, nsamp)
                ind1 <- diff(kni, lag = 1) > 0.04/(nsamp^2)
                if (sum(ind1) > 0) 
                  kni[c(FALSE, ind1) | c(ind1, FALSE)] <- 1
                if (nsamp >= 3) {
                  ind2 <- diff(kni, lag = 2) > 0.25/(nsamp^1.5)
                  if (sum(ind2) > 0) 
                    kni[c(FALSE, FALSE, ind2) | c(FALSE, ind2, 
                      FALSE) | c(ind2, FALSE, FALSE)] <- 1
                }
                tmp <- .Fortran("denlocal", as.double(lastunif), 
                  kni = as.integer(kni), as.integer(nsamp), icomax = as.integer(icomax), 
                  irmax = as.integer(irmax), PACKAGE = "ftnonpar")
                if (sum(tmp$kni) == 0) 
                  break
                else eps[tmp$kni == 1] <- eps[tmp$kni == 1] * 
                  locsq.factor
            }
            else break
        }
        if (verbose) {
            print("Press Enter")
            readline()
        }
    }
    list(y = x.string, widthes = upper - fdist.y, nmax = fts$nmax, 
        ind = fts$knotsind, trans = lastunif)
}

"kuipdiffbounds" <-
structure(c(0.307501980728653, 0.220314152964105, 0.157498224961261, 
0.101218629121418, 0.0717806990191302, 0.0508098522102511, 0.0321599217349666, 
0.260975723191012, 0.188679698351668, 0.137057382929534, 0.0881050582317245, 
0.0633510908581634, 0.0449285678974546, 0.0285894565096712, 0.148978663288177, 
0.110940438285337, 0.080922168300817, 0.0525782956802672, 0.0374876315772114, 
0.0269628371478656, 0.0170712581613634, 0.124382620372917, 0.0952923972152614, 
0.069502435660574, 0.0452402857283996, 0.0328627738951723, 0.0234861641266666, 
0.0150209787920118, 0.0931979549885054, 0.0721989123993246, 0.0533914989167988, 
0.0353843290154652, 0.0259823805766242, 0.0185174716975161, 0.0117425052969115, 
0.080181867346804, 0.0628086817228552, 0.0476545800307207, 0.0316471006759641, 
0.0230176879538813, 0.0164825579792217, 0.010551360533595, 0.0661650592780638, 
0.0526433140545621, 0.0404992596271707, 0.0270039151375531, 0.019671300684033, 
0.0141403144918236, 0.00911192652213594, 0.0576892057496199, 
0.047291120062897, 0.0363065131491498, 0.0244300390839729, 0.017840792812333, 
0.0128944894170847, 0.0083147603774937, 0.0497274415042684, 0.0412695182684986, 
0.0319415193967662, 0.0216317776162515, 0.0159815894358178, 0.0115404374512996, 
0.0074691950384487, 0.0440930205642902, 0.0371632683985669, 0.0291764170236889, 
0.0199711130191836, 0.0147851456503776, 0.0108183394040534, 0.00692649071583678, 
0.0386375955093702, 0.0335325558345777, 0.0264992890334411, 0.0183791050868901, 
0.0136438144386473, 0.00985636827251102, 0.00638630984182286, 
0.0350372845631238, 0.0306666019917607, 0.0246776701228097, 0.0171324304025882, 
0.0127252321205716, 0.00922867868345428, 0.00599548336578164, 
0.0316423784749728, 0.0280847488931776, 0.0227942157567389, 0.0158230403649484, 
0.011929053040925, 0.00864599572624527, 0.00564937236258158, 
0.0290331314985375, 0.0259967548049147, 0.0212965842285101, 0.0149711609910706, 
0.0112445050423623, 0.0082202319675742, 0.00535064138925489, 
0.025883090604231, 0.0240572813028095, 0.0200325810700139, 0.014186955977593, 
0.0105454708236573, 0.00776751548779649, 0.00506143758819203, 
0.0232576453467689, 0.0225364463917157, 0.0188249152247203, 0.0134862219744271, 
0.0100617604237537, 0.00742707159158747, 0.00481720923365059, 
0.0206750232753329, 0.0208361558801763, 0.017644365356675, 0.0128059279425379, 
0.00958895744708681, 0.00710294102891937, 0.00459511562898805, 
0.0191363489411528, 0.0194317742601465, 0.0167332157118633, 0.0121466576268726, 
0.00920060897330519, 0.00677157289599328, 0.00443133999240373, 
0.01789137940886, 0.0183482298482047, 0.0159440097924292, 0.0116206451168149, 
0.00879431813865428, 0.0064837615919975, 0.004265070149443), .Dim = c(7, 
19))
"kuipdiffbounds.x" <-
c(50, 100, 200, 500, 1000, 2000, 5000)
"kkuip" <- function (x, k = 1)
{
    tmp <- .C("kkuip", as.double(x), as.integer(length(x)), as.integer(k),
        norm = double(k), a = integer(k), b = integer(k),PACKAGE="ftnonpar")
    list(metric = tmp$norm, a = tmp$a, b = tmp$b)
}
"rtennormal" <-
function (n)
{
    rsamp <- sample(1:10, n, replace = TRUE, prob = rep(0.1, 10))
    mus <- 0.5*(10 * rsamp - 5)
    sigmas <- 1
    rnorm(n, mus, sigmas)
}
"rclaw" <- function (n)
{
    rsamp <- sample(0:5, n, replace = TRUE, prob = c(0.1, 0.1, 0.1,
        0.1, 0.1, 0.5))
    mus <- double(n)
    sigmas <- double(n)
    mus[rsamp != 5] <- rsamp[rsamp != 5]/2 - 1
    mus[rsamp == 5] <- 0
    sigmas[rsamp != 5] <- 0.1
    sigmas[rsamp == 5] <- 1
    rnorm(n, mus, sigmas)
}
"dclaw" <- function (x)
{
    out <- 0.5 * dnorm(x)
    for (i in 0:4) out <- out + 0.1 * dnorm(x, i/2 - 1, 0.1)
    out
}

