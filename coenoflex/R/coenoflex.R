coenoflex <- function (numgrd, numplt, numspc, grdtyp, grdlen, width, variab, 
    grdprd, alphad, pdist, sdist, skew, aacorr, cmpasy, cmpphy, 
    maxtot, noise, slack, autlin) 
{
    rndplt <- function(numplt, numgrd, grdlen, grdprd) {
        centrd <- matrix(0, nrow = numplt, ncol = numgrd)
        pltprd <- rep(0, numplt)
        grdpos <- 0
        tmp <- .Fortran("rndplt", as.integer(numplt), as.integer(numgrd), 
            centrd = as.double(centrd), as.double(grdlen), as.double(grdprd), 
            pltprd = as.double(pltprd), as.double(grdpos), PACKAGE = "coenoflex")
        centrd <- matrix(tmp$centrd, ncol = numgrd)
        pltprd <- tmp$pltprd
        plot.out <- list(centrd = centrd, pltprd = pltprd)
        invisible(plot.out)
    }
    fixplt <- function(numplt, numgrd, grdlen, grdprd) {
        centrd <- matrix(0, nrow = numplt, ncol = numgrd)
        pltprd <- rep(0, numplt)
        numpts <- rep(0, numgrd)
        index <- rep(0, numgrd)
        size <- 0
        expans <- 0
        grdpos <- 0
        totsam <- 0
        tmp <- .Fortran("fixplt", as.integer(numplt), as.integer(numgrd), 
            as.double(grdlen), as.double(grdprd), centrd = as.double(centrd), 
            pltprd = as.double(pltprd), as.double(size), as.double(expans), 
            as.double(grdpos), numpts = as.integer(numpts), as.integer(totsam), 
            as.integer(index), PACKAGE = "coenoflex")
        mat <- matrix(tmp$centrd, ncol = numgrd)
        mat <- mat[1:prod(tmp$numpts), ]
        plot.out <- list(centrd = mat, pltprd = tmp$pltprd[1:prod(tmp$numpts)])
        invisible(plot.out)
    }
    rndspc <- function(numspc, numgrd, grdlen, alphad, width, 
        variab, grdtyp, skew, aacorr) {
        maxabu <- rep(0, numspc)
        spcamp <- array(0, dim = c(numspc, numgrd, 5))
        fudge <- 0
        hcnadj <- 0
        maxval <- 0
        tmp <- .Fortran("rndspc", as.integer(numspc), as.integer(numgrd), 
            spcamp = as.double(spcamp), maxabu = as.double(maxabu), 
            as.double(grdlen), as.double(alphad), as.double(width), 
            as.double(variab), as.integer(grdtyp), as.double(skew), 
            as.double(aacorr), as.double(fudge), as.double(hcnadj), 
            as.double(maxval), PACKAGE = "coenoflex")
        spc.out <- list(spcamp = array(tmp$spcamp, dim = c(numspc, 
            numgrd, 5)), maxabu = tmp$maxabu)
        invisible(spc.out)
    }
    fixspc <- function(numspc, numgrd, grdlen, width, variab, 
        grdtyp, skew, aacorr) {
        spcamp <- array(0, dim = c(numspc, numgrd, 5))
        maxabu <- rep(0, numspc)
        numpts <- rep(0, numgrd)
        index <- rep(0, numgrd)
        size <- 0
        expans <- 0
        center <- 0
        fudge <- 0
        hcnadj <- 0
        tmp <- .Fortran("fixspc", as.integer(numspc), as.integer(numgrd), 
            spcamp = as.double(spcamp), maxabu = as.double(maxabu), 
            as.double(grdlen), as.double(width), as.double(variab), 
            as.integer(grdtyp), as.double(skew), as.double(aacorr), 
            as.double(size), as.double(expans), numpts = as.integer(numpts), 
            as.integer(index), as.double(center), as.double(fudge), 
            as.double(hcnadj), PACKAGE = "coenoflex")
        spcamp <- array(tmp$spcamp, dim = c(numspc, numgrd, 5))
        spcamp <- spcamp[1:prod(tmp$numpts), , ]
        maxabu <- tmp$maxabu[1:prod(tmp$numpts)]
        spc.out <- list(spcamp = spcamp, maxabu = maxabu)
        invisible(spc.out)
    }
    totphy <- function(numplt, numspc, numgrd, noise, slack, 
        cmpasy, cmpphy, centrd = plot$centrd, spcamp = spc$spcamp, 
        argmnt = autinfo$argmnt, grdlst = autinfo$grdlst, numper = autinfo$numper, 
        count = autinfo$count, maxabu = spc$maxabu, pltprd = plot$pltprd) {
        physio <- matrix(0, nrow = numspc, ncol = numgrd + 10)
        abunda <- matrix(0, nrow = numplt, ncol = numspc)
        diff <- rep(0, numspc)
        tmp <- .Fortran("totphy", as.integer(numplt), as.integer(numspc), 
            as.integer(numgrd), as.double(centrd), as.double(spcamp), 
            as.double(physio), as.integer(argmnt), as.integer(grdlst), 
            as.integer(numper), as.integer(count), as.double(maxabu), 
            abunda = as.double(abunda), as.double(pltprd), as.double(noise), 
            as.double(slack), as.double(maxtot), as.double(cmpasy), 
            as.double(cmpphy), as.double(diff), PACKAGE = "coenoflex")
        abunda <- matrix(tmp$abunda, nrow = numplt)
        abunda
    }

    autpar <- function (line) 
    {
        maxnst <- function(xline)
        {
            numlft <- 0
            maxnst <- 0
            for (i in 1:length(xline)) {
                if (xline[i] == '(') {
                    numlft <- numlft + 1
                    maxnst <- max(maxnst,numlft)
                }
                else if (xline[i] == ')') {
                    numlft <- numlft - 1
                }
            }
            if (numlft != 0) stop('unbalanced parenetheses')
            return(maxnst)
        }
    
        extract <- function(xline,maxnst,derived)
        {
            opers <- c('ave','min','max','geo','irm')
            numlft <- 0
            grads <- NULL
            for (i in 1:length(xline)) {
                if (xline[i] == '(') {
                    numlft <- numlft + 1
                    if (numlft == maxnst) {
                        for (j in (i+1):length(xline)) {
                            if (xline[j] == ')') {
                                oper <- paste(xline[i-3],xline[i-2],xline[i-1],sep='')
                                oper <- which(opers == oper)
                                for (k in (i+1):(j-1)) {
                                    if (xline[k] != ',') grads <- c(grads,as.numeric(xline[k]))
                                }
                                xline[i-3] <- as.character(derived)
                                xline <- xline[-c((i-2):j)]
                                out <- list(oper=oper,grads=grads,xline=xline) 
                                return(out)
                            }
                        }
                    }
                }
            }
        }
    
    
        xline <- strsplit(line,'')[[1]]
    
        parts <- 0
        derived <- 11
        numper <- rep(0,10)
        count <- 0
        grdlst <- matrix(0,nrow=10,ncol=10)
        argmnt <- rep(0,10)
        rownum <- 1
    
        while(length(grep('(',xline,fixed=TRUE)) > 0) {
            nest <- maxnst(xline)
            ext <- extract(xline,nest,derived)
            argmnt[rownum] <- ext$oper
            grdlst[rownum,1:length(ext$grads)] <- ext$grads
            count <- count + 1
            numper[rownum] <- length(ext$grads)
            xline <- ext$xline
            derived <- derived + 1
            rownum <- rownum + 1
       }
       out <- list(argmnt=argmnt,grdlst=grdlst,numper=numper,count=count)
       out
    }
    auteco <- function(nusmpc = numspc, numgrd = numgrd, argmnt = autinfo$argmnt, 
        grdlst = autinfo$grdlst, numper = autinfo$numper, count = autinfo$count, 
        physio = physio) {
        .Fortran("auteco", as.integer(numspc), as.integer(numgrd), 
            as.integer(argmnt), as.integer(grdlst), as.integer(numper), 
            as.integer(count), as.double(physio), PACKAGE = "coenoflex")
    }
    hlpgrd <- function() {
        cat("\n  Coenoflex simulates up to 10 gradients.  The\n", 
            " number of gradients should be few compared to\n", 
            " the number of plots, so that the distribution\n", 
            " of plots is not too sparse.\n\n")
        resp <- readline("Enter a number between 1 and 10 : ")
        return(resp)
    }
    hlpspc <- function() {
        cat("\n  Coenoflex simulates the distribution of species\n", 
            " along gradients.  The number of species should be\n", 
            " chosen to obtain a reasonable species richness\n", 
            " given the number of gradients.\n\n")
        resp <- readline(" Enter a maximum of 1000 : ")
        return(resp)
    }
    hlpplt <- function() {
        cat("\n  Coenoflex simulates the vegetation composition of\n", 
            " sample points along environmental gradients.  Enter\n", 
            " a sufficient number of points to generate a suitable\n", 
            " density given the number of gradients\n\n")
        resp <- readline(" Enter a maximum of 1000 : ")
        return(resp)
    }
    hlpgtp <- function() {
        cat("\n  Coenoflex uses two types of gradients, following Austin\n", 
            " and Smith (1987).  The first type of gradient is an\n", 
            " environmental gradient which includes both direct\n", 
            " environmental gradients and complex environmental\n", 
            " gradients sensu Austin and Smith.  Species have modal\n", 
            " responses along these gradients with minimum, optimum, and\n", 
            " maximum values along the gradient (although not necessarily\n", 
            " within the range of gradient being sampled).  Species have\n", 
            " independent optima chosen randomly from a sampling\n", 
            " distribution depending on the alpha-diversity parameter \n", 
            " for that gradient.\n\n", " Resource gradients differ in that all species have their\n", 
            " optimum at the high end of the gradient, with\n", 
            " independently varying minima.  Resource gradients should\n", 
            " have a productivity increase along them, specified with\n", 
            " the gradient productivity coefficient.\n\n")
        resp <- readline(" Enter e or r : ")
        return(resp)
    }
    hlpglt <- function() {
        cat("\n  Coenoflex calculates gradient lengths in arbitrary units.\n", 
            " These units are used to compare gradient lengths in\n", 
            " relative values, so that for example a gradient of 200\n", 
            " units is twice as long as a gradient of 100 units.\n", 
            " In addition, species amplitudes are measured in the same\n", 
            " units so that a species with an amplitude of 50 would\n", 
            " extend about one half the length of a gradient of 100\n", 
            " units.\n\n", " Thus, species mean amplitudes scaled according to gradient\n", 
            " lengths controls species turnover rates.\n")
        resp <- readline(" Enter a length up to 1000 : ")
        return(resp)
    }
    hlpamp <- function() {
        cat("\n  Coenoflex calculates the ecological amplitude of species in\n", 
            " gradient units.  Thus, if the gradient length is 100, and\n", 
            " a species amplitude is 50, it will extend about one half\n", 
            " the length of the gradient.  Species amplitude, along with\n", 
            " the number of species, determines the number of species\n", 
            " per plot and the rate of species turnover.\n\n")
        resp <- readline(" Enter a number in gradient units : ")
        return(resp)
    }
    hlpvar <- function() {
        cat("\n  Coenoflex varies individual species amplitudes according to\n", 
            " the variance parameter.  Previously you specified a mean\n", 
            " species amplitude along each gradient.  The variance\n", 
            " parameter modifies that value for each species by +/- %.\n", 
            " For example, if you specified a mean amplitude of 50 and a\n", 
            " variance of 50%, species amplitudes would vary randomly\n", 
            " from 25 units to 75 units.\n\n")
        resp <- readline(" Enter a variance in percent : ")
        return(resp)
    }
    hlpprd <- function() {
        cat("\n  Coenoflex can vary the productivity along gradients\n", 
            " systematically.  To increase productivity along the\n", 
            " gradient, enter the percent increase in productivity\n", 
            " from the low to high end of the gradient.  Resource\n", 
            " gradients should have a significant increase associated\n", 
            " with them; environmental gradients may or may not\n\n")
        resp <- readline(" Enter a percent increase from 0 to 200% : ")
        return(resp)
    }
    hlpald <- function() {
        cat("\n  Coenoflex can systematically vary the alpha-diversity along\n", 
            " gradients by modifying the location of species centroids.\n", 
            " Values from 0.5 to 1.0 increase the alpha-diversity at the\n", 
            " high end of the gradient, and values from 1.0 to 2.0\n", 
            " increase alpha-diversity at the low end of the gradient.\n", 
            " Values near 1.0 have small effects; the farther from 1.0\n", 
            " the more dramatic the effect.\n\n")
        resp <- readline(" Enter number between 0.5 and 2.0 : ")
        return(resp)
    }
    hlpamp <- function() {
        cat("\n  Coenoflex calculates the ecological amplitude of species in\n", 
            " gradient units.  Thus, if the gradient length is 100, and\n", 
            " a species amplitude is 50, it will extend about one half \n", 
            " the length of the gradient.  Species amplitude, along with\n", 
            " the number of species, determines the number of species\n", 
            " per plot and the rate of species turnover.\n\n")
        resp <- readline(" Enter a number in gradient units : ")
        return(resp)
    }
    hlpprf <- function() {
        cat("\n  Coenoflex can establish plot centers randomly or along\n", 
            " a grid.  If you select randomly, plot centers will have\n", 
            " random coordinates along each gradient with approximately\n", 
            " constant density throughout the sample space\n\n", 
            " If you select a grid, Coenoflex will attempt to space plots\n", 
            " along the gradients in proportion to their length, again\n", 
            " attempting to maintain approximately constant density.\n", 
            " If necessary, Coenoflex will change the number of plots\n", 
            " specified to achieve equal spacing.\n\n")
        resp <- readline(" Enter a g (for grids) or r (for random) : ")
        return(resp)
    }
    hlpskw <- function() {
        cat("\n  Coenoflex simulates the maximum abundance of organisms\n", 
            " according an equation similar to a log-random\n", 
            " distribution.  The skew coefficient controls how skewed\n", 
            " the distribution is toward low values.  A skew of 1.0\n", 
            " results in a normal distribution with a mean of 50% cover.\n", 
            " As the value increases, the mode moves left and the tail to\n", 
            " the right grows longer.\n\n")
        resp <- readline(" Enter value between 1.0 and 4.0 : ")
        return(resp)
    }
    hlphcn <- function() {
        cat("\n  Coenoflex allows species abundances to be correlated\n", 
            " with their ecological amplitude as specified in the\n", 
            " Hierarchical Continuum Model (Collins et al. 1993\n", 
            " J. Veg. Sci. 4:149-156.) \n", " Values greater than 0.0 result in positive correlation\n", 
            " between amplitude and abundancem while 0.0 results\n", 
            " in independence\n\n")
        resp <- readline(" Enter [0.00-1.00] : ")
        return(resp)
    }
    hlpsrf <- function() {
        cat("\n  Coenoflex can establish species modes randomly or along\n", 
            " a grid.  If you select randomly, species modes will have\n", 
            " random coordinates along each gradient with approximately\n", 
            " constant density throughout the sample space.\n", 
            " If you select a grid, Coenoflex will attempt to space\n", 
            " species modes along the gradients in proportion to their\n", 
            " length, again attempting to maintain approximately\n", 
            " constant density.  If necessary, Coenoflex will change\n", 
            " the number of species specified to achieve equal spacing.\n\n")
        resp <- readline(" Enter a g (for grids) or r (for random) : ")
        return(resp)
    }
    hlpasy <- function() {
        cat("\n  Coenoflex has several mechanisms to simulate\n", 
            " competition.  The competition asymmetry coefficient\n", 
            " simulates a competitive advantage to larger plants\n", 
            " intended to simulate asymmetric competition such as\n", 
            " for light.  The algorithm scales species abundance\n", 
            " proportional to raw abundance raised to the competition\n", 
            " asymmetry coefficient.  Accordingly, a value of 1.0 \n", 
            " confers no advantage, and increasingly larger values\n", 
            " give large species more advantage.\n\n")
        resp <- readline(" Enter a number between 1.0 and 3.0): ")
        return(resp)
    }
    hlpnoi <- function() {
        cat("\n  Coenoflex allows you to simulate noise in the species\n", 
            " response model.  Values are +/- in percent.  For example,\n", 
            " a response of 25 means plus or minus 25% of the actual\n", 
            " value.\n\n")
        resp <- readline(" Enter a number between 0 and 100 : ")
        return(resp)
    }
    hlpslk <- function() {
        cat("\n  Coenoflex includes a form of noise known as slack.  In\n", 
            " contrast to noise, slack causes a species to be absent\n", 
            " in a favorable environment, due perhaps to poor dispersal\n", 
            " (or incorrect identification in the field).  The values\n", 
            " are probabilities of absence in favorable conditions, for\n", 
            " example a value of 0.1 means a species will be absent in\n", 
            " 10% of favorble plots.\n\n")
        resp <- readline(" Enter a value between 0.00 and 1.00 : ")
        return(resp)
    }
    hlpaut <- function() {
        cat("\n  Coenoflex combines the physiological response of a species\n", 
            " along each gradient into an autecological response using\n\n", 
            " min(x,y,...) = minimum response of gradients x,y,... \n", 
            " max(x,y,...) = maximum response of gradients x,y,... \n", 
            " irm(x,y,...) = integrated rate methodology response \n", 
            "                                 of gradients x,y,... \n", 
            " geo(x,y,...) = geometric mean response of gradients \n", 
            "                                              x,y,... \n", 
            " ave(x,y,...) = arithmetic mean rsponse of gradients \n", 
            "                                 of gradients x,y,... \n\n", 
            " functions can be nested as e.g.  min(1,irm(2,3)) \n\n")
        resp <- readline(" Enter a function : ")
        return(resp)
    }
    hlpstd <- function() {
        cat("\n  Coenoflex includes the option to standardize the mean\n", 
            " total plot abundance to a specific value.  If the\n", 
            " gradients have productivity coefficients specified, this\n", 
            " mean value will appear at the center of the simulated\n", 
            " coenospace.\n\n", " In order for the competition asymmetry algorithm to work\n", 
            " a total abundance must be set.  However, in some cases it\n", 
            " may mask unrealistic vegetation models (too sparse or\n", 
            " too dense).\n\n")
        resp <- readline(" Enter Y or N : ")
        return(resp)
    }
    hlptot <- function() {
        cat("\n  Enter the total cover for samples in percent.  Values\n", 
            " of 100 simulate a typical standarization of plot total = 100\n", 
            " but other values are possible\n\n")
        resp <- readline(" Enter the value in percent : ")
        return(resp)
    }
    if (missing(numgrd)) 
        numgrd <- readline(" enter the number of gradients                   [1-10] : ")
    if (numgrd == "?") 
        numgrd <- hlpgrd()
    numgrd <- as.numeric(numgrd)
    if (numgrd <= 0 | numgrd > 10) {
        cat("\n You can only enter a number between 1 and 10\n")
        numgrd <- as.numeric(readline(" enter the number of gradients                   [1-10] : "))
    }
    if (missing(grdtyp)) 
        grdtyp <- rep(NA, numgrd)
    else if (length(grdtyp) != numgrd) 
        stop("You must specify a gradient type for each gradient")
    grdtyp[grdtyp == "E" | grdtyp == "e"] <- 1
    grdtyp[grdtyp == "R" | grdtyp == "r"] <- 2
    if (missing(grdlen)) 
        grdlen <- as.numeric(rep(NA, numgrd))
    else if (length(grdlen) != numgrd) 
        stop("You must specify a gradient length for each gradient")
    if (missing(width)) 
        width <- rep(NA, numgrd)
    else if (length(width) != numgrd) 
        stop("You must specify a mean species amplitude for each gradient")
    if (missing(variab)) 
        variab <- rep(NA, numgrd)
    else if (length(variab) != numgrd) 
        stop("You must specify species amplitude variability for each gradient")
    if (missing(grdprd)) 
        grdprd <- rep(NA, numgrd)
    else if (length(grdprd) != numgrd) 
        stop("You must specify the gradient productivity for each gradient")
    if (missing(alphad)) 
        alphad <- rep(NA, numgrd)
    else if (length(alphad) != numgrd) 
        stop("You must specify the gradient alpha diversity for each gradient")
    for (i in 1:as.numeric(numgrd)) {
        if (is.na(grdtyp[i])) {
            cat(paste("\n enter gradient type for gradient ", 
                i, "\n"))
            resp <- readline(" enter e for environment or r for resource     [e or r] : ")
            if (resp == "?") 
                resp <- hlpgtp()
            if (resp == "e" || resp == "E") 
                grdtyp[i] <- 1
            else if (resp == "r" || resp == "R") 
                grdtyp[i] <- 2
            else {
                cat("\n You can only enter e or r\n")
                cat(paste("\n enter gradient type for gradient ", 
                  i, "\n"))
                resp <- readline(" enter e for environment or r for resource     [e or r] : ")
                if (resp == "e" || resp == "E") 
                  grdtyp[i] <- 1
                else grdtyp[i] <- 2
            }
        }
        if (is.na(grdlen[i])) {
            grdlen[i] <- readline(paste(" enter gradient length for gradient ", 
                i, "       [10-1000] : "))
            if (grdlen[i] == "?") 
                grdlen[i] <- hlpglt()
            grdlen <- as.numeric(grdlen)
            if (grdlen[i] < 10 || grdlen[i] > 1000) {
                cat("\n Enter a number between 10 and 1000\n")
                grdlen[i] <- as.numeric(readline(paste(" enter gradient length for gradient ", 
                  i, "       [10-1000] : ")))
            }
        }
        if (is.na(width[i])) {
            cat("\n enter the mean width of species amplitudes\n")
            width[i] <- readline("      reasonable values are ~[0.25-3 x gradient length] : ")
            if (width[i] == "?") 
                width[i] <- hlpamp()
            width[i] <- as.numeric(width[i])
        }
        if (is.na(variab[i])) {
            variab[i] <- readline(" enter the variability in percent                [0-50] : ")
            if (variab[i] == "?") 
                variab[i] <- hlpvar()
            variab[i] <- as.numeric(variab[i])
        }
        if (is.na(grdprd[i])) {
            grdprd[i] <- readline(" enter the productivity response in percent     [0-200] : ")
            if (grdprd[i] == "?") 
                grdprd[i] <- hlpprd()
            variab[i] <- as.numeric(variab[i])
        }
        if (is.na(alphad[i])) {
            alphad[i] <- readline(" enter the trend in alpha-diversity           [0.5-2.0] : ")
            if (alphad[i] == "?") 
                alphad[i] <- hlpald()
            alphad[i] <- as.numeric(alphad[i])
        }
    }
    if (missing(numplt)) {
        numplt <- readline("\n enter the number of plots                      [1-500] : ")
        if (numplt == "?") 
            numplt <- hlpplt()
        numplt <- as.numeric(numplt)
        if (numplt <= 0 || numplt > 500) {
            cat("\n enter a number between 1 and 500\n")
            numplt <- as.numeric(readline("\n enter the number of plots                      [1-500] : "))
        }
    }
    if (missing(pdist)) {
        pdist <- readline(" random or grid plot locations                 [r or g] :  ")
        if (pdist == "?") 
            pdist <- hlpprf()
        if (pdist == "r" || pdist == "R") {
            plot <- rndplt(numplt, numgrd, grdlen, grdprd)
        }
        else {
            plot <- fixplt(numplt, numgrd, grdlen, grdprd)
        }
    }
    else {
        if (pdist != "r" && pdist != "R" && pdist != "g" && pdist != 
            "G") 
            stop("The plot distribution must be R (random) or G (grid)")
        if (pdist == "r" || pdist == "R") {
            plot <- rndplt(numplt, numgrd, grdlen, grdprd)
        }
        else {
            plot <- fixplt(numplt, numgrd, grdlen, grdprd)
        }
    }
    if (missing(numspc)) {
        numspc <- readline("\n enter the number of species                    [1-500] : ")
        if (numspc == "?") 
            numspc <- hlpspc()
        numspc <- as.numeric(numspc)
        if (numspc < 1 || numspc > 500) {
            cat("\n enter a number between 1 and 500\n")
            numspc <- as.numeric(readline("\n enter the number of species                    [1-500] : "))
        }
    }
    if (missing(skew)) {
        skew <- readline(" enter the skew                             [1.00-4.00] : ")
        if (skew == "?") 
            skew <- hlpskw()
        skew <- as.numeric(skew)
        if (skew < 1 || skew > 4) {
            cat("\n enter a number between 1.0 and 4.0\n")
            skew <- as.numeric(readline(" enter the skew                             [1.00-4.00] : "))
        }
    }
    if (missing(aacorr)) {
        aacorr <- readline(" enter the amplitude/abundance correlation  [0.00-1.00] : ")
        if (aacorr == "?") 
            aacorr <- hlphcn()
        aacorr <- as.numeric(aacorr)
        if (aacorr < 0 || aacorr > 1) {
            cat("\n enter a number between 0.0 and 1.0\n")
            aacorr <- as.numeric(readline(" enter the amplitude/abundance correlation  [0.00-1.00] : "))
        }
    }
    if (missing(sdist)) {
        sdist <- readline(" random or grid species modes                  [r or g] : ")
        if (sdist == "?") 
            sdist <- hlpsrf()
        if (sdist == "r" || sdist == "R") {
            spc <- rndspc(numspc, numgrd, grdlen, alphad, width, 
                variab, grdtyp, skew, aacorr)
        }
        else {
            spc <- fixspc(numspc, numgrd, grdlen, width, variab, 
                grdtyp, skew, aacorr)
        }
    }
    else {
        if (sdist != "r" && sdist != "R" && sdist != "G" && sdist != 
            "g") 
            stop("The species distribution must be R (random) or G (grid)")
        else if (sdist == "r" || sdist == "R") 
            spc <- rndspc(numspc, numgrd, grdlen, alphad, width, 
                variab, grdtyp, skew, aacorr)
        else spc <- fixspc(numspc, numgrd, grdlen, width, variab, 
            grdtyp, skew, aacorr)
    }
    if (missing(cmpasy)) {
        cmpasy <- readline(" enter the competition asymmetry            [1.00-3.00] : ")
        if (cmpasy == "?") 
            cmpasy <- hlpasy()
        cmpasy <- as.numeric(cmpasy)
        if (cmpasy < 1 || cmpasy > 3) {
            cat("\n enter a number between 1.0 and 3.0\n")
            cmpasy <- as.numeric(readline(" enter the competition asymmetry            [1.00-3.00] : "))
        }
    }
    cmpphy <- 0
    if (missing(maxtot)) {
        yorn <- readline(" standardize mean total plot abundance?        [y or n] : ")
        if (yorn == "?") 
            yorn <- hlpstd()
        if (yorn == "y" || yorn == "Y") {
            maxtot <- readline(" enter the standard total percent cover        [10-500] : ")
            if (maxtot == "?") 
                maxtot <- hlptot()
            maxtot <- as.numeric(maxtot)
            if (maxtot < 10 || maxtot > 500) {
                cat("\n enter a number between 10 and 500\n")
                maxtot <- as.numeric(readline(" enter the standard total percent cover        [10-500] : "))
            }
        }
        else {
            maxtot <- 0
        }
    }
    if (missing(noise)) {
        noise <- readline(" enter the noise in percent                     [0-100] : ")
        if (noise == "?") 
            noise <- hlpnoi()
        noise <- as.numeric(noise)
        if (noise < 0 | noise > 100) {
            cat("\n enter a number between 0 and 100\n")
            noise <- as.numeric(readline(" enter the noise in percent                     [0-100] : "))
        }
    }
    if (missing(slack)) {
        slack <- readline(" enter the slack as a fraction              [0.00-1.00] : ")
        if (slack == "?") 
            slack <- hlpslk()
        slack <- as.numeric(slack)
        if (slack < 0 | slack > 1) {
            cat("\n enter a number between 0.0 and 1.0\n")
            slack <- as.numeric(readline(" enter the slack as a fraction              [0.00-1.00] : "))
        }
    }
    if (missing(autlin)) {
        autlin <- readline(" enter the autecological function\n : ")
        if (autlin == "?") 
            autlin <- hlpaut()
    }
    autinfo <- autpar(autlin)
    abunda <- totphy(numplt, numspc, numgrd, noise, slack, cmpasy, 
        cmpphy)
    params <- list(numplt = numplt, numspc = numspc, numgrd = numgrd, 
        grdtyp = grdtyp, grdlen = as.numeric(grdlen), width = as.numeric(width), 
        variab = as.numeric(variab), grdprd = as.numeric(grdprd), 
        alphad = as.numeric(alphad), skew = skew, aacorr = aacorr, 
        cmpasy = cmpasy, cmpphy = cmpphy, pltdst = pdist, spcdst = sdist, 
        funct = autlin, noise=noise, slack = slack)
    site <- matrix(plot$centrd, nrow = numplt)
    veg <- matrix(abunda, ncol = numspc)
    spcamp <- array(spc$spcamp, dim = c(numspc, numgrd, 5))
    maxabu <- spc$maxabu
    pltprd <- plot$pltprd
    out <- list(params = params, site = site, veg = veg, spcamp = spcamp, 
        maxabu = maxabu, pltprd = pltprd)
    class(out) <- "coenoflex"
    invisible(out)
}

plot.coenoflex <- function (x, which = "all", ...)
{
    if (class(x) != "coenoflex")
        stop("must pass an object of class coenoflex")
    if (which == "plots" | which == "all") {
        if (x$params$numgrd == 1) {
            plot(x$site[,1],xlab='Sample',ylab='Gradient Location')
            readline("hit return to continue")
        }
        else if (x$params$numgrd > 1) {
            for (i in 1:(x$params$numgrd - 1)) {
                for (j in (i + 1):x$params$numgrd) {
                    plot(x$site[, i], x$site[, j], asp = 1, xlab = paste("Gradient",
                      i), ylab = paste("Gradient", j), main = "Plot Locations")
                    readline("hit return to continue")
                }
            }
        }
    }
    if (which == "species" | which == "all") {
        if (x$params$numgrd == 1) {
            plot(x$spcamp[,1,3],xlab='Species',ylab='Gradient Location')
            readline("hit return to continue")
        }
        else if (x$params$numgrd > 1) {
            for (i in 1:(x$params$numgrd - 1)) {
                for (j in (i + 1):x$params$numgrd) {
                    plot(x$spcamp[, i, 3], x$spcamp[, j, 3], asp = 1,
                      xlab = paste("Gradient", i), ylab = paste("Gradient",
                        j), main = "Species Centroids")
                    lines(c(0, x$params$grdlen[i], x$params$grdlen[i],
                      0, 0), c(0, 0, x$params$grdlen[j], x$params$grdlen[j],
                      0), col = 2)
                    readline("hit return to continue")
                }
            }
        }
    }
    if (which == "amplitude" | which == "all") {
        for (i in 1:x$params$numgrd) {
            axis <- 0:x$params$grdlen[i]
            plot(c(0, x$params$grdlen[i]), c(0, max(x$maxabu)),
                xlab = paste("Gradient", i), ylab = "Abundance",
                type = "n")
            for (j in 1:x$params$numspc) {
                val <- x$spcamp[j, i, ]
                curve <- rep(0, x$params$grdlen[i] + 1)
                for (k in 0:x$params$grdlen[i]) {
                  if (k < val[1])
                    curve[k + 1] <- 0
                  else if (val[1] < k & k <= val[2])
                    curve[k + 1] <- 2 * ((k - val[1])/(val[3] -
                      val[1]))^2
                  else if (val[2] < k & k <= val[3])
                    curve[k + 1] <- 1 - 2 * ((val[3] - k)/(val[3] -
                      val[1]))^2
                  else if (val[3] < k & k <= val[4])
                    curve[k + 1] <- 1 - 2 * ((k - val[3])/(val[5] -
                      val[3]))^2
                  else if (val[4] <- k & k <= val[5])
                    curve[k + 1] <- 2 * ((val[5] - k)/(val[5] -
                      val[3]))^2
                  else curve[k + 1] <- 0
                }
                curve <- curve * x$maxabu[j]
                lines(axis, curve, col = (j%%6) + 2, pty = "l")
            }
            readline("hit return to continue")
        }
    }
}

