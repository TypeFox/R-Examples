############################################################################################
## package 'secr'
## autoini.R
## last changed 2009 04 21, 2009 07 12 (failure R2.9)
## 2009 07 17 add f.lower, f.upper to uniroot call for small speed improvement
## 2009 10 27 changes to dbar for robustness with counts etc.
## 2009 11    RPSV
## 2009 12 11 update integralprw1 call for nmix=1
## 2010 02 14 update integralprw1 call for spherical (no longer used)
## 2010 08 28 traps renamed to trps
## 2010 11 01 naivesigma shifted outside autoini
## 2011 01 24 minor editing
## 2012 04 19 added bad data check in naivesigma
## 2012 12 14 computeD added but suppressed for now
## 2012 12 15 tol argument added
## 2012 12 24 binomN and adjust.g0 arguments added
## 2014-08-27 dist2 optional input to integralprwi set to -1
## 2014-08-29 linear adjustment (mask 'spacing') when no mask 'area'
## 2015-10-04 markocc argument for integralprw1
## 2015-11-24 fixed markocc-related bug and use RPSV CC = TRUE
############################################################################################

naivesigma <- function (obsRPSV, trps, mask, wt, detectfn, z, tol = 0.001) {
    naiveRPSVcall <- function (sigma)
    {
        temp <- .C ("naiveRPSV", PACKAGE = 'secr',
          as.double (sigma),               # Parameter : detection scale
          as.double (z),                   # Parameter : detection shape (fixed)
          as.integer(k),
          as.integer(m),
          as.integer (wt),
          as.double (unlist(trps)),        # x,y locations of traps (first x, then y)
          as.double (unlist(mask)),        # x,y points on mask (first x, then y)
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value
        )
        if (temp$value > 0)
            obsRPSV - temp$value
        else
            obsRPSV                        # dummy value
    }

    if (is.null(mask))
        mask <- make.mask(trps, buffer = 10 * obsRPSV, nx=32)
    k <- nrow(trps)
    m <- nrow(mask)

    if (obsRPSV <= 0)   ## added 2012-04-19
        temp <- NA
    else
        temp <- try(uniroot (naiveRPSVcall, lower = obsRPSV/10, upper = obsRPSV*10,
                             tol=tol)$root)
    ifelse (inherits(temp,'try-error'), NA, temp)
}

autoini <- function (capthist, mask, detectfn = 0, thin = 0.2, tol = 0.001,
                     binomN = 1, adjustg0 = TRUE, ignoreusage = FALSE)

# obtain approximate fit of HN SECR model
# for use as starting values in MLE
# uses external C code

{
    naivedcall <- function (sigma)
    {
        temp <- .C ("naived",  PACKAGE = 'secr',
          as.double (sigma),               # Parameter : detection scale
          as.integer(k),
          as.integer(m),
          as.integer (wt),
          as.double (unlist(trps)),        # x,y locations of traps (first x, then y)
          as.double (unlist(mask)),        # x,y points on mask (first x, then y)
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value
        )
        db-temp$value
    }
    naivecap2 <- function (g0, sigma, cap)
    # using new algorithm for number of captures per animal 2009 04 21
    {
        temp <- .C ("naivecap2",  PACKAGE = 'secr',
          as.integer(prox),
          as.double (g0),                  # Parameter : detection magnitude
          as.double (sigma),               # Parameter : detection scale
          as.integer(s),                   # number of occasions
          as.integer(k),                   # number of traps
          as.integer(m),                   # number of points in mask
          as.integer(wt),                  # trap weights
          as.double (unlist(trps)),        # x,y locations of traps (first x, then y)
          as.double (unlist(mask)),        # x,y points on mask (first x, then y)
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value
        )
        cap - temp$value
    }
    naiveesa <- function (g0, sigma)
    {
      nc <- 1
      g0sigma0 <- matrix(rep(c(g0,sigma), c(2,2)), nrow = 2)
      gs0 <- rep(1,2*s*k)
      #--------------------------------------------
      # 2012-12-15 allow for incomplete grid
      if (!is.null(usage(trps))) {
          ## only considers binary use/non-use
          gs0[rep(t(usage(trps)==0), each = 2)] <- -1
      }
      #--------------------------------------------
      area <- attr(mask, 'area')  # area of single mask cell
      if (is.null(area))   ## 2014-08-29 linear
          area <- attr(mask, 'spacing')/1000
      miscparm <- c(1,0,0,0)  ## dummy value
      temp <- try ( .C("integralprw1",  PACKAGE = 'secr',
          as.integer(dettype),
          as.double(g0sigma0),
          as.integer(rep(1,nc)),        # group number 2013-06-24
          as.integer(nc),
          as.integer(s),
          as.integer(k),
          as.integer(m),
          as.integer(1),                # number of groups 2012-11-13
          as.integer(1),
          as.integer(rep(1,nc)),        # knownclass 2013-04-12, 2013-06-04
          as.double(unlist(trps)),
          as.double(-1),                ## optional dist2 2014-08-27
          as.double(unlist(usge)),
          as.integer(markocc),          # 2015-10-04
          as.double(unlist(mask)),
          as.integer(nrow(g0sigma0)),   # rows in lookup
          as.integer(gs0),              # index of nc+1,S,K to rows in g0sigma0
          as.integer(1),
          as.double(area),
          as.double(miscparm),
          as.integer(detectfn),
          as.integer(0),                # binomN
          as.integer(0),                # useD
          a=double(nc),
          resultcode=integer(1))
      )

      if (temp$resultcode != 0)
          stop ("error in external function 'integralprw1'; ",
                "possibly the mask is too large")
      temp$a
    }

    #############
    ## main line
    computeD <- TRUE
    if (length(tol)==1) tol <- rep(tol,2)

    if (nrow(capthist)<5)
        stop ("too few values for autoini")  ## message changed 2015-01-06

    ## added 2010-07-01
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (! (detectfn %in% c(0)))
        stop ("only halfnormal detection function implemented in 'autoini'")

    trps <- traps(capthist)
    usge <- usage(trps)
    if (is.null(usge) | ignoreusage)
        ## assuming k = nk i.e. not polygon or transect detector
        usge <- matrix(1, nrow = nrow(trps), ncol = ncol(capthist))
    dettype <- detectorcode(trps)
    n       <- nrow(capthist)    # number of individuals
    s       <- ncol(capthist)    # number of occasions
    k       <- nrow(trps)
    markocc <- markocc(traps(capthist))
    if (is.null(markocc)) 
        markocc <- rep(1,s)     ## assume all occasions were simple marking occasions
    allsighting <- !any(markocc>0)

    if (!(dettype %in% c(-1:5,8)))
        list(D=NA, g0=NA, sigma=NA)
    else {
        prox     <- detector(trps) %in% c('proximity', 'count','signal','times')
        ## wt is the number of opportunities for capture given binary usage
        wt <- apply(usge>0, 1, sum)

        # optionally thin mask
        if ((nrow(mask)>100) & (thin>0) & (thin < 1))
            mask <- mask[runif(nrow(mask)) < thin,]
        else
            thin <- 1.0
        m        <- nrow(mask)                   # number of points in mask
        if (length(dim(capthist))>2)
            cpa     <- sum(abs(capthist))/n      # captures per animal
        else
            cpa     <- sum(abs(capthist)>0)/n    # captures per animal

        obsRPSV <- RPSV(capthist, CC = TRUE)

        if (is.na(obsRPSV) | (obsRPSV<1e-10)) {    ## try db
            db <- dbar(capthist)
            if (!is.null(attr(trps,'spacing'))) {
                if (is.na(db)) {
                    warning ("could not calculate 'dbar'; using detector spacing")
                    db <- attr(trps, 'spacing')
                }
                if (db < (attr(trps, 'spacing')/100)) {
                    warning ("'dbar' close to zero; using detector spacing instead")
                    db <- attr(trps, 'spacing')
                }
            }
            if (is.na(db) | is.nan(db) | (db<1e10) )
                return(list(D=NA,g0=NA,sigma=NA))
            else
                tempsigma <- uniroot (naivedcall, lower = db/10, upper = db*10, tol=tol[2])$root
        }
        else {
            ## 2015-11-24 simpler to use RPSV
            ## tempsigma <- naivesigma(obsRPSV = obsRPSV, trps = trps, mask = mask, wt = wt,
            ## detectfn=detectfn, z=1, tol = tol[2])
            tempsigma <- obsRPSV
        }
        if (is.null(usage(trps))) wt <- rep(s,k)
        low <- naivecap2(0.001, sigma=tempsigma, cap=cpa)
        upp <- naivecap2(0.999, sigma=tempsigma, cap=cpa)
        badinput <- FALSE
        if (is.na(low) | is.na(upp)) badinput <- TRUE
        else if (sign(low) == sign(upp)) badinput <- TRUE
        if (badinput) {
            # not sure what conditions cause this 28/4/2008
            # observed number cap more than expected when g0=1 28/8/2010
            # maybe better in future to set g0 = 0.9
            warning ("'autoini' failed to find g0; setting initial g0 = 0.1")
            tempg0 <- 0.1
        }
        else tempg0 <- uniroot (naivecap2, lower=0.001, upper=0.999,
                 f.lower = low, f.upper = upp, tol=tol[1],
                 sigma=tempsigma, cap=cpa)$root
        if (computeD) {
            if (allsighting)  ## includes all-zero rows
                tempD <- n / maskarea(mask)
            else {
                esa <- naiveesa (tempg0, tempsigma)
                tempD <-  n / esa * thin
            }
        }
        else tempD <- NA

        ## 2012-12-18,24 adjust for large effort and/or binomN

        if (adjustg0) {
            if (binomN == 0)
                adjusted.g0 <- tempg0 / usge[usge>0]
            else
                adjusted.g0 <- 1 - (1 - tempg0)^ ( 1 / (usge[usge>0] * binomN) )
            tempg0 <- mean(adjusted.g0)
        }

        list(D = tempD, g0 = tempg0, sigma = tempsigma)
    }
}
##################################################################################
