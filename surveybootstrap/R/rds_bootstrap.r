

#####################################################
##' draw RDS bootstrap resamples for one chain
##'
##' this function uses the algorithm described in the
##' supporting online material for Weir et al 2012
##' (TODO PROPER CITE) to take bootstrap resamples
##' of one chain from an RDS dataset
##'
##' @param chain the chain to draw resamples for
##' @param mm the mixing model to use
##' @param dd the degree distns to use
##' @param parent.trait a vector whose length is the number
##' of bootstrap reps we want
##' @param idvar the name of the variable used to label the
##' columns of the output (presumably some id identifying the
##' row in the original dataset they come from -- see below)
##' @return a list of dataframes with one entry for each respondent in the chain.
##' each dataframe has one row for each bootstrap replicate. so if we take 10
##' bootstrap resamples of a chain of length 50, there will be 50 entries in
##' the list that is returned. each entry will be a dataframe with 10 rows.
rds.boot.draw.chain <- function(chain, mm, dd, parent.trait, idvar="uid") {

    thisid <- paste0(idvar, ".", chain$data[,idvar])

    ## choose the traits of the referrals from this set of parents
    trait.draws <- mm$choose.next.state.fn(parent.trait)

    ## choose the degree (and other vars) of each referral based on
    ## the drawn set of traits
    deg.draws <- dd$draw.degrees.fn(trait.draws)

    ## if no children, this is the last wave
    ## return a list whose only entry is dataframe of draws for this wave
    if (is.null(chain$children)) {
        return(setNames(list(deg.draws), thisid))
    }

    ## otherwise, if there are children, recursively get bootstrap
    ## draws for them
    child.draws <- plyr::llply(chain$children,
                         rds.boot.draw.chain,
                         mm=mm,
                         dd=dd,
                         parent.trait=trait.draws)
    child.draws <- unlist(child.draws, recursive=FALSE)

    
    return(setNames(c(list(deg.draws), child.draws), 
                    c(thisid, names(child.draws))))

}

#####################################################
##' draw RDS bootstrap resamples
##'
##' draw boostrap resamples for an RDS dataset, using
##' the algorithm described in the supporting online
##' material of Weir et al 2012 (TODO PROPER CITE)
##'
##' TODO -- consider constructing chains, mm from other args
##'
##' TODO be sure to comment the broken-out trait variables
##'      (ie these could all be different from the originals)
##'
##' @param chains a list whose entries are the chains
##' we want to resample
##' @param mm the mixing model
##' @param dd the degree distributions
##' @param num.reps the number of bootstrap resamples we want
##' @param keep.vars if not NULL, then the names of variables
##' from the original dataset we want appended to each bootstrap
##' resampled dataset (default is NULL)
##' @return a list of length \code{num.reps}; each entry in
##' the list has one bootstrap-resampled dataset
rds.chain.boot.draws <- function(chains,
                                 mm,
                                 dd,
                                 num.reps,
                                 keep.vars=NULL) {

    traits <- mm$traits

    ## get the bootstrap resamples for each chain
    res <- plyr::llply(chains,
                 function(this.chain) {

                     ## TODO -- should handle the case where there is
                     ## missingness in a seed trait?

                     seed.trait <- traits.to.string(this.chain$data[,traits,drop=FALSE],
                                                    traits)$traits

                     return(rds.boot.draw.chain(this.chain,
                                                mm,
                                                dd,
                                                rep(seed.trait,num.reps)))
                 })

    ## within each chain, we now have a list whose entries are dataframes,
    ## one for each repsondent, with one row for each bootstrap resample
    ## convert this into a list whose entries are dataframes, one for each
    ## bootstrap resample, and whose rows are respondents
    res.byboot <- plyr::llply(res,
                        function(this.chain.res) {
                            by.rep <- plyr::llply(1:num.reps,
                                            function(this.rep.id) {
                                                plyr::ldply(this.chain.res,
                                                      function(x) x[this.rep.id,])
                                            })
                            return(by.rep)
                        })

    ## assemble the bootstrap resamples from each chain together
    ## to end up with a list whose entries are datasets, one for each bootstrap
    ## resample (across all chains)
    res.dat <- plyr::llply(1:num.reps,
                     function(this.rep.id) {
                         plyr::ldply(res.byboot,
                              function(x) { 
                                  this.dat <- x[[this.rep.id]] 
                                  ## remove the .id column
                                  this.dat$.id <- NULL
                                  ## add the trait columns back in...
                                  return(cbind(this.dat,
                                               unparse.trait(this.dat$trait,
                                                             traits)))
                              })
                     })

    return(res.dat)

}

#####################################################
##' draw RDS bootstrap resamples using the
##' algorithm in Salganik 2006 (TODO PROPER CITE)
##'
##' this algorithm picks a respondent from the survey
##' to be a seed uniformly at random. it then generates
##' a bootstrap draw by simulating the markov process
##' forward for n steps, where n is the size of the draw
##' required.
##'
##' if you wish the bootstrap dataset to end up with
##' variables from the original dataset other than the
##' traits and degree, then you must specify this when
##' you construct \code{dd} using the 
##' '\code{estimate.degree.distns} function.
##'
##'
##' TODO be sure to comment the broken-out trait variables
##'      (ie these could all be different from the originals)
##'
##' @param chains a list with the chains constructed from the survey
##' using \code{make.chain}
##' @param mm the mixing model
##' @param dd the degree distributions
##' @param num.reps the number of bootstrap resamples we want
##' @return a list of length \code{num.reps}; each entry in
##' the list has one bootstrap-resampled dataset
rds.mc.boot.draws <- function(chains,
                              mm,
                              dd,
                              num.reps) {

    traits <- mm$traits

    all.data <- plyr::ldply(chains,
                      chain.data)

    ## NB: only using traits that are nonmissing
    all.traits <- traits.to.string(all.data[,traits,drop=FALSE], traits)
    all.traits.str <- all.traits$traits

    chain.sizes <- plyr::laply(chains,
                         chain.size)

    num.chains <- length(chain.sizes)

    res <- plyr::llply(1:num.reps,
                 function(this.rep) {

                     ## NB: an alternative would be do to the bootstrap this way,
                     ## but this is not what has been done in the literature so far
                     ## so it's commented out
                     ##
                     ## draw chains with the same length as the
                     ## empirically observed ones
                     ## draw seeds (with replacement, right?)
                     ## these.seeds <- sample(all.traits.str, size=num.chains, replace=TRUE)
                     ##
                     ## these.traits <- unlist(llply(1:num.chains,
                     ##                              function(x) {
                     ##                                  res <- mc.sim(mm,
                     ##                                                these.seeds[x],
                     ##                                                chain.sizes[x])
                     ##                                  return(res)
                     ##                              }))

                     ## use only one chain
                     this.seed <- sample(all.traits.str, size=1)

                     ## draw the chain of respondent traits, starting from our seed,
                     ## using the markov model we were given
                     these.traits <- mc.sim(mm,
                                            this.seed,
                                            sum(chain.sizes))

                     ## draw degrees for each of the respondents in the chain we
                     ## just drew
                     these.degs <- dd$draw.degrees.fn(these.traits)

                     return(cbind(these.degs,
                                  unparse.trait(these.traits, traits)))
                 })

    return(res)

}
