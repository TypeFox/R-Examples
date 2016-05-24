################################################################################
###
### TITLE:              pop_reconstruction_functions.R
###
### AUTHOR:             Mark Wheldon
###
### DESC:               Functions to do population reconstruction developed for
###                     "Reconstructing Past Populations with Uncertainty from
###                     Fragmentary Data" submitted to Journal of the American
###                     Statistical Assocation".
###
###-----------------------------------------------------------------------------
################################################################################


### ------------------------ CCMPP FUNCTION- ---------------------- ###
### --------------------------------------------------------------- ###

popRecon.ccmp.female <-
    function(pop, surv, fert, srb = 1.05, mig
             ,proj.steps, age.int = 5
             ,label.dims = FALSE, base.year = "1960"
             )

    #-- Based on Preston et al. (2001), Ch. 6 --#
    #
    #   pop     :  population count at baseline
    #   fert    :  matrix of age specific fertility rates NOT yet
    #                mulitplied by age.int
    #   srb     :  sex ratio at birth matrix
    #   surv    :  Survivorship probabilities: the probability of
    #                reaching the age at the start of the interval.
    #              The first row should be nL0/(n*l0).
    #              The last row is survival for age.int years in the open
    #                interval
    #   mig     :  Net number of migrants as a
    #                 _proportion_ of prev time period's population
    #   proj.steps
    #           :  number of time periods to project forward
    #                if missing, set to ncol(fert)
    #   age.int :  needed for correct interpretation of survival
    #                and fertility rates
    #   label.dims
    #           :  should output have dimnames set?
    #   base.year
    #           :  start year for projections (aesthetic)

{

    #-- Checks --#

    # If proj.steps is greater than the number of columns in surv,
    #   reduce or recycle fert, surv matrices


    #-- Constants --#

    n.age.grps <- length(pop)
    n.surv <- nrow(surv)


    #-- Derive proj.steps from ncol(fert) --#

    if(missing(proj.steps)) proj.steps <- ncol(fert)


    #-- Make SRB into matrix if not already --#

    if(is.null(dim(srb))) srb <- matrix(rep(srb, proj.steps))


    #-- Loop over number of required time periods --#
    #...............................................#

    #-- Initialise pop.matrix --#

    pop.mat <- matrix(0, nrow = n.age.grps, ncol = 1 + proj.steps)
    pop.mat[,1] <- pop


    #-- Initialize leslie matrix --#

    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)


    #-- Project --#
    #.............#

    for(i in 1:proj.steps)
    {

        #-- Make the leslie matrix --#
        #............................#

        #-- Rows 2:(n.age.grps) = survival ratios --#

        lesM[2:n.age.grps,1:(n.age.grps-1)] <-
            diag(surv[-c(1,n.surv),i])
        lesM[n.age.grps,n.age.grps] <- surv[(n.surv),i]

        #-- First row = fert and survival ratios --#

        k <- 1/(1+srb[i]) * surv[1,i] * 0.5

        dbl.fert <- age.int*fert[,i] + c(age.int*fert[-1,i], 0) *
            surv[-1,i]

        lesM[1,] <- k * dbl.fert

        #-- Migrants --#

       net.numb.mig <- mig[,i] * pop.mat[,i]


        #-- Project --#
        #.............#

        pop.mat[,i+1] <-
            lesM %*% (pop.mat[,i] + 0.5 * net.numb.mig) +
                0.5 * net.numb.mig

    }


    #-- Add dim names --#
    #...................#

    if(label.dims) {
        ages <- seq(from = 0
                    ,to = age.int*(nrow(as.matrix(pop))-1)
                    ,by = age.int)
        yrs <- (0:proj.steps) * age.int + as.numeric(base.year)
        dimnames(pop.mat) <- list(ages, yrs)
    }


    #-- Output --#
    #............#

    return(pop.mat)

}


### -------------------- SUMMARIZATION FUNCTIONS------------------- ###
### --------------------------------------------------------------- ###

life.expectancy.stationary <- function(z)
{
    x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
    5 * sum(cumprod(x))
}


make.leslie.matrix <-
    function(pop, surv, fert, srb = 1.05, age.int = 5, label.dims = FALSE)

    ##-- Make the leslie matrix for CCMPP --##
    ##
    ##   pop     :  population count at baseline
    ##   fert    :  matrix of age specific fertility rates NOT yet
    ##                mulitplied by age.int
    ##   srb     :  sex ratio at birth matrix
    ##   surv    :  Survivorship probabilities: the probability of
    ##                reaching the age at the start of the interval.
    ##              The first row should be nL0/(n*l0).
    ##              The last row is survival for age.int years in the open
    ##                interval
    ##   proj.steps
    ##   age.int :  needed for correct interpretation of survival
    ##                and fertility rates
    ##   label.dims
    ##           :  should output have dimnames set?

{
    ## Constants
    n.age.grps <- length(pop)
    n.surv <- length(surv)

    ## Make Leslie matrix
    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)

    ## first row = fert and birth survival
    k <- 1/(1+srb) * surv[1] * 0.5
    dbl.fert <- age.int*fert + c(age.int*fert[-1], 0) * surv[-1]
    lesM[1,] <- k * dbl.fert

    ## rows 2:(n.age.grps) = survival ratios
    lesM[2:n.age.grps,1:(n.age.grps-1)] <- diag(surv[-c(1,n.surv)])
    lesM[n.age.grps,n.age.grps] <- surv[n.surv]

    if(label.dims) {
        age.labs <- seq(from = 0, by = 5, length = n.age.grps)
        dimnames(lesM) <- list(age.labs, age.labs)
    }

    ## return
    return(lesM)
}

net.number.migrants <- function(n1, n2, L)
{
    ##-- Find net number of migrants in a CCMPP projection --##
    ##
    ## ARGUMENTS
    ##
    ##   n1      :  Population count vector at time t
    ##   n2      :  Population count vector at time t + delta
    ##   L       :  Leslie matrix used to get population at t + delta
    ##
    ##
    ## METHOD
    ##
    ## Invert n2 = L(n1 + 0.5 mig) + (0.5)*mig
    ## Can get proportions by pre-multiplying output by 'solve(diag(n1))'

    ## Make sure inputs are of correct form
    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    L <- as.matrix(L)

    return(2 * solve(L + diag(nrow(L))) %*% (n2 - L %*% n1))
}



### ----------------------- HELPER FUNCTIONS ---------------------- ###
### --------------------------------------------------------------- ###


## ........... Misc Functions .......... ##
## ..................................... ##

estMod.logit.mar29 <- function(p) log(p / (1 - p))

estMod.invlogit.mar29 <- function(x)
{
    if(any(is.infinite(exp(x)))) {
        y <- x
        y[is.infinite(exp(x))] <- 1
        y[!is.infinite(exp(x))] <-
            estMod.invlogit.mar29(y[!is.infinite(exp(x))])
        return(y)
    }
    else return(exp(x) / (1 + exp(x)))
}


##--- Generates random draws from inverse gamma ---##

estMod.rinvGamma.mar29 <- function(n, shape, scale)
{
    return(1/rgamma(n, shape = shape, rate = scale))
}


##--- Returns value of inverse gamma pdf ---##

estMod.dinvGamma.mar29 <- function(x, shape, scale, log = FALSE)
{
    if(log) d <-
        shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
    else d <- scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
    return(d)
}


##--- Creates column names for mcmc objects ---##

estMod.makeColNames.mar29 <- function(m)
{
    ## m:    matrix of input values

    e <- expand.grid(rownames(m), colnames(m))
    apply(e, 1, FUN = function(z) paste(z[2], z[1], sep = "."))

}


## .... Projection for census years .... ##
## ..................................... ##

proj.cen.yrs <-
    function(full.proj, bline.yr, vr.yrs, cen.yrs, proj.yrs
             ,labels = FALSE)

{
    bline.yr <- as.numeric(bline.yr)
    vr.yrs <- as.numeric(vr.yrs)
    cen.yrs <- as.numeric(cen.yrs)
    proj.yrs <- as.numeric(proj.yrs)

    interp.counts <- matrix(NA, nrow = nrow(full.proj), ncol = length(cen.yrs))

    if(labels) {
        if(!is.null(rownames(full.proj))) {
            rownames(interp.counts) <- rownames(full.proj)
        } else {
            rownames(intper.counts) <- seq(from = 0, by = 5
                                       ,length.out = nrow(full.proj)
                                       )
        }
        colnames(interp.counts) <- cen.yrs
    }

    cen.in.proj <- cen.yrs %in% proj.yrs
    cen.notin.proj <- !cen.in.proj
    proj.in.cen <- proj.yrs %in% cen.yrs


    ## Check to see if interpolation necessary

    if(all(cen.in.proj)) {
        interp.counts <- full.proj[,proj.in.cen]
    } else {

        ## Growth rates
        ## (Does this for all, including those not needing interpolation)

        ratio.mat <- matrix(NA,
                            ,nrow = nrow(full.proj)
                            ,ncol = ncol(full.proj) - 1)
        for(j in 1:ncol(ratio.mat)) {
            ratio.mat[,j] <- full.proj[,(j+1)] / full.proj[,j]
        }

        yr.diffs <- diff(proj.yrs)

        ## if(any(is.nan(log(ratio.mat)))) {
        ##     cat("\nratio.mat=\n");  print(ratio.mat)
        ##     stop()
        ## }
        r.mat <- log(ratio.mat) / yr.diffs


        ## Duplicate last column of r.mat so that "interpolation"
        ## can be done for years beyond the end of proj.yrs,
        ## assuming that the growth rate remains constant
        ## (although don't know why we'd do this)

        r.mat <- cbind(r.mat, r.mat[,ncol(r.mat)])


        ## Interpolate

        for(j in which(cen.notin.proj)) {
            base.yr.ind <-
                max(which(proj.yrs <= cen.yrs[j]))
            time.gap <- cen.yrs[j] - proj.yrs[base.yr.ind]
            growth <- r.mat[,base.yr.ind]
            interp.counts[,j] <-
                full.proj[,base.yr.ind] * exp(growth * time.gap)
        }
        interp.counts[,cen.in.proj] <- full.proj[,proj.in.cen]

    }

    return(interp.counts)

}


## ............. Likelihood ............ ##
## ..................................... ##

log.lhood.mar29 <-
    function(log.n.census, log.n.hat, ll.var)
{
    ##.. log.n.census and log.n.hat should already be logged

    ##.. log.n.hat should be log projected counts for census years
    ##   interpolated if necessary


    ##-- value of log likelihoods --##

    density <- dnorm(log.n.census,
                     mean = log.n.hat,
                     sd = sqrt(ll.var),
                     log = TRUE
                     )

    ##-- joint log likelihood --##

    return(sum(density))
}


## .............. Posterior ............ ##
## ..................................... ##

log.post.mar29 <- function(## estimated vitals
                           f, s, g, baseline.n
                           ## fixed prior means on vitals
                           ,prior.mean.f, prior.mean.s
                           ,prior.mean.g, prior.mean.b
                           ## fixed prior parameters on variance distns
                           ,alpha.f, beta.f, alpha.s, beta.s
                           ,alpha.g, beta.g
                           ,alpha.n, beta.n
                           ## updated variances on prior distns
                           ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n
                           ## value of the log likelihood
                           ,log.like
                           ## non zero rows of fertility matrix
                           ,non.zero.fert
                           )
{

    ##-- Values of prior densities for vitals --##

    ##.. Note that log densities are calculated for numerical stability.
    ##     f, baseline.n, prior.mean.f, prior.mean.b are logged coming
    ##     in, s, prior.mean.s is logit transformed coming in, g and
    ##     prior.mean.g are not transformed coming in.
    log.f.prior <- dnorm(as.vector(f[non.zero.fert,])
                         ,mean = as.vector(prior.mean.f[non.zero.fert,])
                         ,sd = sqrt(sigmasq.f)
                         ,log = TRUE)
    log.s.prior <- dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                         ,log = TRUE)
    log.g.prior <- dnorm(g, mean = prior.mean.g
                         ,sd = sqrt(sigmasq.g)
                         ,log = TRUE)
    log.b.prior <- dnorm(baseline.n, mean = prior.mean.b
                         ,sd = sqrt(sigmasq.n)
                         ,log = TRUE)

    ##-- Values of prior densities for variances --##

    log.sigmasq.f.prior <-
        log(estMod.dinvGamma.mar29(sigmasq.f, alpha.f, beta.f))
    log.sigmasq.s.prior <-
        log(estMod.dinvGamma.mar29(sigmasq.s, alpha.s, beta.s))
    log.sigmasq.g.prior <-
        log(estMod.dinvGamma.mar29(sigmasq.g, alpha.g, beta.g))
    log.sigmasq.n.prior <-
        log(estMod.dinvGamma.mar29(sigmasq.n, alpha.n, beta.n))


    ##-- The log posterior is the SUM of these with the log.like --##

    return(sum(log.f.prior, log.s.prior, log.g.prior, log.b.prior
               ,log.sigmasq.f.prior
               ,log.sigmasq.s.prior
               ,log.sigmasq.g.prior
               ,log.sigmasq.n.prior
               ,log.like))

}


## ......... Acceptance Ratio .......... ##
## ..................................... ##

acc.ra.mar29 <- function(log.prop, log.current)
{
    min(1, exp(log.prop - log.current))
}

acc.ra.var.mar29 <-
    function(log.prop.post, log.curr.post, log.prop.var, log.curr.var)
{
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post
               ))
}



### --------------------------- SAMPLER --------------------------- ###
### --------------------------------------------------------------- ###

popRecon.sampler <-
    function(#.. number of iterations and burn-in (not saved)
             n.iter, burn.in = 0, thin.by = 1

             #.. fixed variance hyper-parameters
             ,al.f = 1, be.f = 0.0109, al.s = 1, be.s = 0.0109
             , al.g = 1, be.g = 0.0436, al.n = 1, be.n = 0.0109

             #.. fixed prior means
             ,mean.f, mean.s, mean.g, mean.b

             #.. inital values for vitals and variances
             #   *vitals not transformed coming in*
             ,start.f = mean.f, start.s = mean.s, start.g = mean.g, start.b = mean.b
             ,start.sigmasq.f = 5, start.sigmasq.s = 5, start.sigmasq.g = 5
             ,start.sigmasq.n = 5

             #.. census data
             #   *not transformed coming in*
             ,pop.data

             #.. **variances** for proposal distributions used in M-H
             #   steps which update vital rates.
             ,prop.vars

             #.. ccmp function
             ,ccmp.function = popRecon.ccmp.female

             #.. number of periods to project forward over (e.g.,
             #     number of five-year steps)
             ,proj.periods = ncol(mean.f)

             #.. age group width
             ,age.size = 5

             #.. print algorithm progress
             ,verb = FALSE

             #.. tolerance defining allowable survival probabilities
             ,s.tol = 10^(-10)
             )
{

    ## .............. Sampler .............. ##
    ## ..................................... ##


    ## -------- Begin timing ------- ##

    ptm <- proc.time()


    ## ------ Match functions ------ ##


    ## ------- Check input dimensions ------- ##

    input.dims <- sapply(list(start.s = start.s[1:(nrow(start.s)-1),]
                              , start.g = start.g
             ,mean.f = mean.f, mean.s = mean.s[1:(nrow(mean.s)-1),]
             ,mean.g = mean.g
             )
           ,"dim")
    mismatch.dims <- apply(input.dims, 2, "identical", dim(start.f))
    if(!all(mismatch.dims))
        stop("Dims of these inputs do not match 'dim(start.f)'", "\n"
             ,paste(names(mismatch.dims)[!mismatch.dims], collapse = "  "))


    ## ------- Check Years --------- ##

    all.vr.years <-
        list(start.s = colnames(start.s), start.g = colnames(start.g)
             ,mean.f = colnames(mean.f), mean.s = colnames(mean.s)
             ,mean.g = colnames(mean.g)
             )
    mismatch.yrs <- sapply(all.vr.years, "identical", colnames(start.f))
    if(!all(mismatch.dims))
        stop("Years of these inputs do not match years of 'start.f'", "\n"
             ,paste(names(mismatch.yrs)[!mismatch.yrs], collapse = "  "))

    all.vr.years.eq <-
        sapply(all.vr.years, FUN = function(z) all.equal(colnames(start.f), z))
    if(!all(all.vr.years.eq)) {
        stop(paste("colnames(", names(all.vr.years)[min(which(!all.vr.years.eq))], ")"
                   ," != colnames(start.f). There may be more...", sep = "")
             )
    }

    proj.years <-
        seq(from = as.numeric(colnames(start.f)[1]), by = age.size
            ,length = proj.periods + 1)

    if(!all.equal(proj.years[1:ncol(start.f)], as.numeric(colnames(start.f))))
        stop("colnames(start.f) !=  seq(from = as.numeric(colnames(start.f)[1]), by = age.size, length = ncol(start.f))")

    vr.years <- as.numeric(colnames(start.f))
    baseline.year <- as.numeric(colnames(start.b))
    census.years <- as.numeric(colnames(pop.data))


    ## ------- Determine fert.rows --------- ##

    zero.elements <- mean.f == 0
    fert.rows <- as.logical(apply(zero.elements, 1, function(z) !all(z)))

    ## ------- Type of migration data ------- ##

    ## No reason for this anymore; remove later
    mig.string <- "prop"


    ## ---------- Storage ---------- ##

    #.. MCMC objects for posterior samples
    # Samples are stored as 2D arrays for compatibility with coda's
    # mcmc format with iterations as rows, year*age.group as columns.
    # Age.group cycles fastest across columns, e.g.,
    # _____________________________________________________
    #   1960  | 1960  | 1960  | ... | 1965  | 1965  | ...
    #   15.19 | 20.24 | 25.29 | ... | 15.19 | 20.24 | ...
    # 1  --   |  --   |  --   | ... |  --   |  --   | ...
    # 2  --   |  --   |  --   | ... |  --   |  --   | ...
    #   etc.
    # _____________________________________________________

    ## How many stored?
    n.stored <- ceiling(n.iter / thin.by)

      # Fertility
      fert.rate.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = sum(fert.rows) * ncol(start.f))
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(fert.rate.mcmc) <-
          estMod.makeColNames.mar29(start.f[fert.rows,])

      # Survival proportions
      surv.prop.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.s) * ncol(start.s))
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(surv.prop.mcmc) <-
          estMod.makeColNames.mar29(start.s)

      # lx
      lx.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.b) * (proj.periods))
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(lx.mcmc) <-
          estMod.makeColNames.mar29(matrix(0, nrow = nrow(start.b)
                                 ,ncol = proj.periods
          ,dimnames = list(rownames(start.b)
           ,seq(from = as.numeric(colnames(start.b)[1]) +
           age.size, by = age.size, length = proj.periods))
                                 )
                          )

      # migration proportions
      mig.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.g) * ncol(start.g))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(mig.mcmc) <-
          estMod.makeColNames.mar29(start.g)

      # baseline counts
      baseline.count.mcmc <-
          mcmc(matrix(nrow = n.stored, ncol = nrow(start.b))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(baseline.count.mcmc) <- estMod.makeColNames.mar29(start.b)

      # variances
      variances.mcmc <-
          mcmc(matrix(nrow = n.stored, ncol = 4)
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(variances.mcmc) <-
          c("fert.rate.var", "surv.prop.var", "mig.var"
            ,"population.count.var")

    #.. Record acceptance rate

    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.g = 0
             ,sigmasq.n = 0
             )


    #.. Count how often acceptance ratio missing or na

    ar.na <- acc.count


    #.. Count how often projection gives negative population

    pop.negative <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             )


    #.. Count how often surv probs are outside tolerance

    s.out.tol <- matrix(0, nrow = nrow(mean.s), ncol = ncol(mean.s)
                        ,dimnames = dimnames(mean.s))


    ## -------- Initialize -------- ##

    #.. Set current vitals and variances to inital values
    #   Take logs/logits here where required

    log.curr.f <- log(start.f) #<-- log(0) stored as "-Inf". Gets
    log.prop.f <- log(start.f) #    converted to 0 under exponentiation
    logit.curr.s <- estMod.logit.mar29(start.s)
    curr.g <- start.g
    log.curr.b <- log(start.b)

    curr.sigmasq.f <- start.sigmasq.f
    curr.sigmasq.s <- start.sigmasq.s
    curr.sigmasq.g <- start.sigmasq.g
    curr.sigmasq.n <- start.sigmasq.n


    #.. Fixed means for vitals and baseline
    #   Set these to inputs, take logs where required.

    log.mean.f <- log(mean.f)
    logit.mean.s <- estMod.logit.mar29(mean.s)
    mean.g <- mean.g
    log.mean.b <- log(mean.b)


    #.. Fixed census data
    #   Take logs here

    log.census.mat <- log(pop.data)


    #.. Set current projection: base on initial values

    log.curr.proj <-
        log(proj.cen.yrs(full.proj =
                         ccmp.function(pop = exp(log.curr.b),
                                       surv = estMod.invlogit.mar29(logit.curr.s),
                                       fert = exp(log.curr.f),
                                       mig = curr.g,
                                       proj.steps = proj.periods,
                                       age.int = age.size)
                         ,bline.yr = baseline.year
                         ,vr.yrs = vr.years
                         ,cen.yrs = census.years, proj.yrs = proj.years
                         ))


    #.. Current log posterior

    log.curr.posterior <-
        log.post.mar29(f = log.curr.f
                       ,s = logit.curr.s
                       ,g = curr.g
                       ,baseline.n = log.curr.b
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.g = mean.g
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.g = al.g, beta.g = be.g
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.g = curr.sigmasq.g
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood.mar29(
                        log.n.census = log.census.mat
                        ,log.n.hat = log.curr.proj
                        ,ll.var = curr.sigmasq.n)
                       ,non.zero.fert = fert.rows
                       )


    ## -------- Begin loop ------- ##
    #...............................#

    if(verb) {
        cat("\n\ntotal iterations = ", n.iter+burn.in
            ,"\nburn in = ", burn.in
            ,"\nthin = ", thin.by
            ,",\nnumber stored = ", n.stored, sep = "")
        cat("\n\nfert.rows = ", which(fert.rows)
            ,"\ncensus years = ", census.years
            ,"\nvital rate years = ", vr.years
            ,"\nprojection years = ", proj.years
            )
        cat("\n\n"
            ,"iter ", " quantity\n", "---- ", " --------"
            ,sep = "")
          }

    for(i in 1:(n.iter + burn.in)) {

      # k is the index into the storage objects
      k <- (i - burn.in - 1) / thin.by + 1


      ## -------- Vital Rate M-H Steps ------- ##

      ##...... Fertility .....##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Fertility")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(log.curr.f[fert.rows,])) {

        #.. make a matrix conformable w fertility rate matrix
        log.prop.f.mat <-
            matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
        log.prop.f.mat[fert.rows,][j] <-
            rnorm(1, 0, sqrt(prop.vars$fert.rate[j]))

        #.. make proposal
        log.prop.f <- log.curr.f + log.prop.f.mat

        # - Run CCMP (project on the original scale)
        #   ** Don't allow negative population
        full.proj <- ccmp.function(pop = exp(log.curr.b),
                            fert = exp(log.prop.f), #<-- use proposal
                            surv = estMod.invlogit.mar29(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                                   age.int = age.size)

        if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
           || is.nan(sum(full.proj))) {
            if(i > burn.in) {
                pop.negative$fert.rate[j] <-
                    pop.negative$fert.rate[j] + 1/n.iter
            }
        } else {
            prop.proj <-
                proj.cen.yrs(full.proj = full.proj
                             ,bline.yr = baseline.year
                             ,vr.yrs = vr.years
                             ,cen.yrs = census.years, proj.yrs = proj.years
                             )
            log.prop.proj <- log(prop.proj)

          # - Calculate log posterior of proposed vital under projection
          log.prop.posterior <-
              log.post.mar29(f = log.prop.f #<-- use proposal
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n)
                             ,non.zero.fert = fert.rows
                             )

          #- Acceptance ratio
          ar <- acc.ra.mar29(log.prop = log.prop.posterior,
                             log.current = log.curr.posterior)

          # - Move or stay
          #.. stay if acceptance ratio 0, missing, infinity, etc.
          if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$fert.rate[j] <-
                ar.na$fert.rate[j] + 1/n.iter
          } else {
            #.. if accept, update current fert rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
              if(i > burn.in) acc.count$fert.rate[j] <-
                  acc.count$fert.rate[j] + 1/n.iter
              log.curr.f <- log.prop.f
              log.curr.proj <- log.prop.proj
              log.curr.posterior <- log.prop.posterior
            }
            #.. if reject, leave current fert rates and projections
            #   alone

          } # close else after checking for ar=0, missing, inf

        } # close else after checking neg or zero population

      } # close loop over all age-spec fertility rates

      #.. Store proposed fertility rate matrix
      if(k %% 1 == 0 && k > 0) fert.rate.mcmc[k,] <-
          as.vector(exp(log.curr.f[fert.rows,]))


      ##...... Survival ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Survival")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.s)) {

        #.. make a matrix conformable w rate matrix
        logit.prop.s.mat <-
            matrix(0, nrow = nrow(logit.curr.s)
                   ,ncol = ncol(logit.curr.s))
        logit.prop.s.mat[j] <- rnorm(1, 0, sqrt(prop.vars$surv.prop[j]))

        #.. make proposal
        logit.prop.s <- logit.curr.s + logit.prop.s.mat

        #.. If proposal resulted in back-transformed s = 0 or 1, do
        #   nothing
        if(estMod.invlogit.mar29(logit.prop.s[j]) > 1 - s.tol ||
           estMod.invlogit.mar29(logit.prop.s[j]) < s.tol) {
          #.. leave current surv rates and projections
          #   alone (simply do not propose
          #   extreme survival probabilities)
          s.out.tol[j] <- s.out.tol[j] + 1/n.iter
        } else {

          # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
            full.proj <- ccmp.function(pop = exp(log.curr.b),
                                       fert = exp(log.curr.f),
                                       surv = estMod.invlogit.mar29(logit.prop.s), #<-- use prop
                                       mig = curr.g,
                                       proj.steps = proj.periods,
                                       age.int = age.size)

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] <-
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
                prop.proj <-
                    proj.cen.yrs(full.proj = full.proj
                                 ,bline.yr = baseline.year
                                 ,vr.yrs = vr.years
                                 ,cen.yrs = census.years, proj.yrs = proj.years
                                 )
            log.prop.proj <- log(prop.proj)

            # - Calculate log posterior of proposed vital under projection
            log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.prop.s #<-- use proposal
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n)
                             ,non.zero.fert = fert.rows
                             )

            #- Acceptance ratio
            ar <- acc.ra.mar29(log.prop = log.prop.posterior,
                               log.current = log.curr.posterior)

            # - Move or stay
            #.. stay if acceptance ratio 0, missing, infinity, etc.
            if(is.na(ar) || is.nan(ar) || ar < 0) {
              if(i > burn.in) ar.na$surv.prop[j] <-
                  ar.na$surv.prop[j] + 1/n.iter
            } else {
              #.. if accept, update current surv rates,
              #   update current projection and count acceptance
              if(runif(1) <= ar) {
                if(i > burn.in) acc.count$surv.prop[j] <-
                    acc.count$surv.prop[j] + 1/n.iter
                logit.curr.s <- logit.prop.s
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
              } #.. if reject, leave current surv rates and projections
                #   alone

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(k %% 1 == 0 && k > 0) surv.prop.mcmc[k,] <-
        as.vector(estMod.invlogit.mar29(logit.curr.s))


      ##...... Migration ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Migration")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(curr.g)) {

        #.. make a matrix conformable w rate matrix
        prop.g.mat <-
            matrix(0, nrow = nrow(curr.g), ncol = ncol(curr.g))
        prop.g.mat[j] <- rnorm(1, 0, sqrt(prop.vars$mig[j]))

        #.. make proposal
        prop.g <- curr.g + prop.g.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
        full.proj <- ccmp.function(pop = exp(log.curr.b),
                              fert = exp(log.curr.f),
                              surv = estMod.invlogit.mar29(logit.curr.s),
                              mig = prop.g, #<-- use proposal
                              proj.steps = proj.periods,
                              age.int = age.size)

      if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
         || is.nan(sum(full.proj))) {
        if(i > burn.in) {
          pop.negative$mig[j] <-
              pop.negative$mig[j] + 1/n.iter
        }
      } else {

        prop.proj <-
            proj.cen.yrs(full.proj = full.proj
                     ,bline.yr = baseline.year
                     ,vr.yrs = vr.years
                     ,cen.yrs = census.years, proj.yrs = proj.years
                     )
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = prop.g #<-- use proposal
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n)
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- acc.ra.mar29(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$mig[j] <-
                ar.na$mig[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$mig[j] <-
                    acc.count$mig[j] + 1/n.iter
                curr.g <- prop.g
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

    } # close loop over all age-specific migration proportions

      #.. Store proposed migration proportion matrix
      if(k %% 1 == 0 && k > 0) mig.mcmc[k,] <- as.vector(curr.g)


      ##...... Baseline population ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

      # - Proposal

      #.. cycle through components (never update last
      #   value as this affects years beyond the estimation period)
      for(j in 1:length(log.curr.b)) {

      #.. make a matrix conformable w rate matrix
      log.prop.b.mat <- matrix(0, nrow = nrow(log.curr.b), ncol = 1)
      log.prop.b.mat[j] <- rnorm(1, 0, sqrt(prop.vars$baseline.pop.count[j]))

      #.. make proposal
      log.prop.b <- log.curr.b + log.prop.b.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
      full.proj <- ccmp.function(pop = exp(log.prop.b), #<-- use proposal
                            fert = exp(log.curr.f),
                            surv = estMod.invlogit.mar29(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size)

      if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
         || is.nan(sum(full.proj))) {
        if(i > burn.in) {
          pop.negative$baseline.count[j] <-
              pop.negative$baseline.count[j] + 1/n.iter
        }
      } else {
      prop.proj <-
          proj.cen.yrs(full.proj = full.proj
                     ,bline.yr = baseline.year
                     ,vr.yrs = vr.years
                     ,cen.yrs = census.years, proj.yrs = proj.years
                     )
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.prop.b #<-- use proposal
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n)
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- acc.ra.mar29(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$baseline.count[j] <-
                ar.na$baseline.count[j] + 1/n.iter
        } else {
            #.. if accept, update current mig rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$baseline.count[j] <-
                    acc.count$baseline.count[j] + 1/n.iter
                log.curr.b <- log.prop.b
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

  } # close loop over all age-specific baseline counts

      #.. Store proposed baseline count matrix
      if(k %% 1 == 0 && k > 0) baseline.count.mcmc[k,] <-
          as.vector(exp(log.curr.b))


      ## ------- Variance Updates ------- ##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Variances")

      ##...... Fertility rate ......##

      prop.sigmasq.f <-
        estMod.rinvGamma.mar29(1, al.f +
                         length(mean.f[fert.rows,])/2,
                  be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                  log.mean.f[fert.rows,])^2)
                  )

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- acc.ra.var.mar29(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = estMod.dinvGamma.mar29(prop.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             ,log.curr.var = estMod.dinvGamma.mar29(curr.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.f <-
                ar.na$sigmasq.f + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.f <-
                    acc.count$sigmasq.f + 1/n.iter
                curr.sigmasq.f <- prop.sigmasq.f
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"fert.rate.var"] <- curr.sigmasq.f


      ##...... Survival Proportion ......##

      prop.sigmasq.s <-
        estMod.rinvGamma.mar29(1, al.s + length(mean.s)/2,
                  be.s +
                    0.5*sum((logit.curr.s - logit.mean.s)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = prop.sigmasq.s  #<-- use proposal
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- acc.ra.var.mar29(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = estMod.dinvGamma.mar29(prop.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             ,log.curr.var = estMod.dinvGamma.mar29(curr.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.s <-
                ar.na$sigmasq.s + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.s <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.s <- prop.sigmasq.s
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"surv.prop.var"] <- curr.sigmasq.s


      ##...... Migration Proportion ......##

      prop.sigmasq.g <-
        estMod.rinvGamma.mar29(1, al.g + length(mean.g)/2,
                  be.g +
                    0.5*sum((curr.g - mean.g)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = prop.sigmasq.g #<-- use proposal
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- acc.ra.var.mar29(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = estMod.dinvGamma.mar29(prop.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             ,log.curr.var = estMod.dinvGamma.mar29(curr.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.g <-
                ar.na$sigmasq.g + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.g <-
                    acc.count$sigmasq.g + 1/n.iter
                curr.sigmasq.g <- prop.sigmasq.g
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"mig.var"] <- curr.sigmasq.g


      ##...... Population Count ......##

      prop.sigmasq.n <-
        estMod.rinvGamma.mar29(1, al.n + (length(mean.b) +
                                    length(log.census.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2) +
                  sum((log.census.mat - log.curr.proj)^2)
                  )
                         )

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post.mar29(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                             ,log.like = log.lhood.mar29(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = prop.sigmasq.n #<-- use proposal
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- acc.ra.var.mar29(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = estMod.dinvGamma.mar29(prop.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.census.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.census.mat - log.curr.proj)^2))
                              ,log = TRUE)
                             ,log.curr.var = estMod.dinvGamma.mar29(curr.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.census.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.census.mat - log.curr.proj)^2))
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.n <-
                ar.na$sigmasq.n + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.n <-
                    acc.count$sigmasq.n + 1/n.iter
                curr.sigmasq.n <- prop.sigmasq.n
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) {
        variances.mcmc[k,"population.count.var"] <- curr.sigmasq.n
      }


      ## ------- Store current population ------- ##

      lx.mcmc[k,] <-
          as.vector(ccmp.function(pop = exp(log.curr.b),
                                  surv = estMod.invlogit.mar29(logit.curr.s),
                                  fert = exp(log.curr.f),
                                  mig = curr.g,
                                  proj.steps = proj.periods,
                                  age.int = age.size
                                  )[-(1:ncol(baseline.count.mcmc))]
                    )


      if(verb && identical(i%%1000, 0)) cat("\n\n")

  } # Ends outer-most loop

    ## ......... End Loop ........ ##
    #...............................#


    ## ---------- Output --------- ##

    #cat("inital values", "\n\n")
    #.. initial values
    start.vals <- list(fert.rate = start.f
                      ,surv.prop = start.s
                      ,mig.XXX = start.g
                      ,baseline.count = start.b
                      ,start.sigmasq.f = start.sigmasq.f
                      ,start.sigmasq.s = start.sigmasq.s
                      ,start.sigmasq.g = start.sigmasq.g
                      ,start.sigmasq.n = start.sigmasq.n
                      ,pop.data = pop.data
                      )

    #.. fixed parameters
    fixed.params <- list(alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.mig.XXX = al.g
                         ,beta.mig.XXX = be.g
                         ,alpha.population.count = al.n
                         ,beta.population.count = be.n
                         ,mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.mig.XXX = mean.g
                         ,mean.baseline.count = mean.b
                         ,mean.pop.data = pop.data
                         )

    #.. migration type
    new.names <-
        lapply(list(start.vals, fixed.params), FUN = function(z) {
        sub("XXX", mig.string, names(z))
    })
    names(start.vals) <- new.names[[1]]
    names(fixed.params) <- new.names[[2]]


    #cat("algorithm statistics", "\n\n")
    #.. algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             ,acc.prop.adj4neg = mapply(FUN = function(a, b, n) {
                 (a * n) / (n - b)
             },
              acc.count[1:4], pop.negative, MoreArgs = list(n = n.iter)
              )
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    #cat("algorithm parameters", "\n\n")
    #.. algorithm parameters
    alg.params <- list(prop.vars = prop.vars
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit", mig.XXX = "I"
                        ,baseline.count = "log"
                        ,population.count = "log")
                       ,projection.periods = proj.periods
                       ,age.gp.size = age.size
                       ,non.zero.fert.rows = fert.rows
                       ,surv.tolerance = s.tol
                       ,burn.in = burn.in
                       ,iters = n.iter
                       ,years = list(vital.rate.years = vr.years
                        ,baseline.year = baseline.year
                        ,census.years = census.years
                        ,projection.years = proj.years
                        )
                       )
    names(alg.params) <- sub("XXX", mig.string, names(alg.params))

    #.. results
    ret.list <- list(fert.rate.mcmc = fert.rate.mcmc
                  ,surv.prop.mcmc = surv.prop.mcmc
                  ,mig.XXX.mcmc = mig.mcmc
                  ,baseline.count.mcmc = baseline.count.mcmc
                  ,lx.mcmc = lx.mcmc
                  ,variances.mcmc = variances.mcmc
                  ,alg.stats = alg.stats
                  ,fixed.params = fixed.params
                  ,start.vals = start.vals
                  ,alg.params = alg.params
                  )
    names(ret.list) <- sub("XXX", mig.string, names(ret.list))

    return(ret.list)

  }
