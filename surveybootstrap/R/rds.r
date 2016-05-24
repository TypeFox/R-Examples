#####################################################
## rds.r
##
## limited implementation of RDS estimator

#####################################################
##' take a set of traits and turn into a string
##'
##' this is a helper function that is useful when we wish
##' to make several traits into one variable
##'
##' @param data the respondent info
##' @param traits the names of the traits to build the model on
##' @param na.action for now, defaults to 'drop' (meaning all rows of data
##' with any missingness on the traits are dropped). anything else
##means
##' NAs are treated like any other value
##' @param sep the separator character used to combine values
##' @return a list whose entries are \code{used.idx}, which indicates
##' which rows from the original dataset were used (may not be all of
##them
##' if there is missingness); and \code{traits}, which has the string
##' version of the traits
##' @keywords internal
traits.to.string <- function(data, traits, na.action="drop", sep=".") {

  if (na.action == "drop") {
    touse.idx <- plyr::aaply(data[, traits],
                       1,
                       function(x) { return(! any(is.na(x))) },
                       .expand=FALSE)
    touse.idx <- which(touse.idx)
  } else {
    touse.idx <- 1:nrow(data)
  }

  traits.str <- plyr::aaply(data[touse.idx, traits],
                      1,
                      paste0,
                      collapse=sep,
                      .expand=FALSE)

  return(list(used.idx=touse.idx,
              traits=traits.str,
              sep=sep,
              names=traits))

}

#####################################################
##' unparse a collapsed trait string
##'
##' for a few of the RDS-related functions, it is useful
##' to combine several traits into one variable as a string;
##' for example, "male" and "young" might become
##' "male.young". this function takes a string with
##' combined traits and explodes it back into
##' several variables
##'
##' @param trait.string a vector whose values are collapsed
##' traits
##' @param names a vector with the names of each trait (in order)
##' @param sep the character used to separate the traits in their
##' collpased string representation
##' @return a dataframe whose rows correspond to the entries in
##' \code{trait.string}, with one column per trait
##' @keywords internal
unparse.trait <- function(trait.string, names, sep="\\.") {

    if (sep == ".") {
        sep <- "\\."
    }

    vals <- str_split(trait.string, sep)

    vals <- plyr::ldply(vals, as.numeric)

    colnames(vals) <- names

    return(vals)

}

#####################################################
##' estimate degree distributions by trait
##'
##' break down RDS degree distributions by trait,
##' and return an object which has the degrees
##' for each trait as well as functions to draw
##' degrees from each trait.
##'
##' @details one of the items returned as a result is a function,
##' \code{draw.degrees.fn}, which takes one argument,
##' \code{traits}. this is a vector of traits and,
##' for each entry in this vector, \code{draw.degress.fn}
##' returns a draw from the empirical distribution of
##' degrees among respondents with that trait. so,
##' \code{draw.degrees.fn(c("0.0", "0.1", "0.1")} would
##' return a degree drawn uniformly at random from among
##' the observed degrees of respondents with trait "0.0"
##' and then two degrees from respondents with trait "0.1"
##'
##' @param survey.data the respondent info
##' @param d.hat.vals the variable that contains
##' the degrees for each respondent
##' @param traits a vector of the names of the columns
##' of \code{survey.data} which refer to the traits
##' @param keep.vars additional vars to return along with degrees
##' @return an object with
##' \itemize{
##'   \item \code{distns} a list with one entry per trait value; each
##entry has a dataframe with all of the degrees from respondents with
##the given trait
##'   \item \code{draw.degrees.fn} a function which gets called with one
##argument, \code{traits}. See description above.
##'   \item \code{keep.vars} the name of the other vars that are kept (if any)
##' }
estimate.degree.distns <- function(survey.data,
                                   d.hat.vals,
                                   traits,
                                   keep.vars=NULL) {

  st <- traits.to.string(survey.data,
                         traits)

  degs <- get.var(survey.data, d.hat.vals)

  ## NOTE: we need to guarantee that the order of deg.dat's columns is
  ##      1: trait
  ##      2: degree
  ## [3...]: keep.vars
  if (! is.null(keep.vars)) {
      ## TODO -- should eventually make grabbing these others vars more robust
      other.vars <- survey.data[st$used.idx, keep.vars]
      deg.dat <- data.frame(trait=st$traits, 
                            degree=degs[st$used.idx],
                            other.vars)
  } else {
      deg.dat <- data.frame(trait=st$traits, degree=degs[st$used.idx])
  }

  ## TODO -- for now, we can just represent the degrees with duplicates
  ## and use SI sampling to pick one when we need to. if there are
  ## huge datasets, this might need to be changed later

  ## to placate R CMD CHECK
  trait <- NULL

  deg.distns <- plyr::dlply(deg.dat,
                      .(trait),
                      identity)

  deg.fns <- unlist(plyr::llply(deg.distns,
                          function(this.trait.deg) {
                              ## if we don't force evaluation here,
                              ## R's lazy evaluation implies that only the
                              ## last version of this.trait.deg will
                              ## get used
                              ## (see, eg,
                              ##  http://adv-r.had.co.nz/Functions.html)
                              force(this.trait.deg)
                              return(function(n=1) {
                                  idx <- sample(1:nrow(this.trait.deg), size=n, replace=TRUE)
                                  return(this.trait.deg[idx,])
                              })
                          }))

  draw.degrees.fn <- function(traits) {
      tocall <- deg.fns[traits]
      degs <- plyr::llply(tocall, do.call, args=list())
      degs <- do.call("rbind", degs)
      rownames(degs) <- NULL
      return(degs)
  }

  return(list(distns=deg.distns, draw.degrees.fn=draw.degrees.fn, keep.vars=keep.vars))

}

#####################################################
##' construct a mixing model from GoC/RDS data
##'
##' given a dataset with the respondents and a dataset
##' on the parents (in many cases the same individuals),
##' and a set of relevant traits,
##' estimate mixing parameters and return a markov model
##'
##' @param survey.data the respondent info
##' @param parent.data the parent info
##' @param traits the names of the traits to build the model on
##' @return a list with two entries:
##' \itemize{
##' \item \code{mixing.df} the data used to estimate the mixing
##function
##' \item \code{choose.next.state.fn} a function which can be passed
##' a vector of states and will return a draw of a subsequent state
##for
##' each entry in the vector
##' \item \code{mixing.df} a dataframe (long-form) representation of
##' the transition counts used to estimate the transition probabilities
##' \item \code{states} a list with an entry for each state. within
##' each state's entry are
##' \itemize{
##' \item \code{trans.probs} a vector of estimated
##' transition probabilities
##' \item \code{trans.fn} a function which,
##' when called, randomly chooses a next state with probabilities given
##' by the transition probs.
##' }}
estimate.mixing <- function(survey.data, parent.data, traits) {

  ## reduce to dataframe with [ child trait, parent trait ]
  ## then basically do a cross tab

  pkey <- attr(parent.data, "key")
  ckey <- attr(survey.data, "key")

  if (is.null(pkey) || is.null(ckey)) {
    stop("parent and survey datasets need to have an attribute which indicates what their keys are")
  }

  st <- traits.to.string(survey.data,
                         traits)
  pt <- traits.to.string(parent.data,
                         traits)

  parent.tomix <- data.frame(key=parent.data[pt$used.idx,pkey],
                             parent.trait=pt$traits)
  survey.tomix <- data.frame(key=survey.data[st$used.idx,ckey],
                             child.trait=st$traits)

  mix.data <- merge(survey.tomix,
                    parent.tomix,
                    by="key",
                    all.x=TRUE)

  res <- list()

  ## we'll reutrn a dataframe which has the parameters we use to estimate the mixing pattern
  res$mixing.df <- as.data.frame(xtabs(~ child.trait + parent.trait, data=mix.data))

  # to placate R CMD CHECK
  parent.trait <- NULL

  ## return a fn which will give us the next step in the chain from each state
  ## based on these transition probabilities
  res$states <- plyr::dlply(res$mixing.df,
                    .(parent.trait),
                    function(this.trait) {

                      if (all(this.trait$Freq == 0)) {
                        return(NULL)
                      }

                      probs <- this.trait$Freq / sum(this.trait$Freq)
                      names(probs) <- this.trait$child.trait

                      return(list(trans.probs=probs,
                                  trans.fn=function() {
                                    draw <- rmultinom(1, 1, probs)
                                    return(names(probs)[as.logical(draw)])
                                  }))
                    })

  ## given a list of preceding states, this function obtains a
  ## succeeding state for each one
  res$choose.next.state.fn <- function(prev.states) {
      tfns <- plyr::llply(res$states, function(x) { x$trans.fn })

      tocall <- tfns[prev.states]

      next.states <- unlist(plyr::llply(tocall, do.call, args=list()))
      return(next.states)

  }

  res$traits <- traits

  return(res)
}


#####################################################
##' run a markov model
##'
##' run a given markov model for n time steps, starting
##' at a specified state
##'
##' this uses the markov model produced by estimate.mixing
##'
##' @param mm the markov model object returned by \code{estimate.mixing}
##' @param start the name of the state to start in
##' @param n the number of time-steps to run through
##' @return a vector with the state visited at each time step. the first entry
##'         has the starting state
##' @export
mc.sim <- function(mm, start, n) {

  if (n <= 1) {
      return(start)
  }

  if (is.null(mm) || is.null(mm$state) || is.null(mm$state[[ start ]])) {
    stop(paste("there was a problem starting at state", start))
  }

  path <- rep(NA, n)

  path[1] <- start

  for(i in 2:n) {
    path[i] <- mm$state[[ path[i-1] ]]$trans.fn()
  }

  return(path)

}

###########################################################
##' determine whether or not one id is a parent of another
##'
##' this function allows us to determine which ids are
##' directly descended from which other ones. it is the only part
##' of the code that relies on the ID format used by the
##' Curitiba study (TODO CITE); by modifying this function,
##' it shold be possible to adapt this code to another study
##'
##' @param id the id of the potential child
##' @param seed.id the id of the potential parent
##' @return TRUE if \code{id} is the direct descendant of \code{seed.id}
##' and FALSE otherwise
is.child.ct <- function(id, seed.id) {

    res <- str_locate(paste(id), paste(seed.id))

    if (! any(is.na(res)) &&
        res[1] == 1 &&
        res[2] == (nchar(id)-1)) {
        return(TRUE)
    }

    return(FALSE)
}

#####################################################
##' build an RDS seed's chain from the dataset
##'
##' text
##' TODO assumes that the chain is a tree (no loops)
##'
##' @param seed.id the id of the seed whose chain we
##' wish to build from the dataset
##' @param survey.data the dataset
##' @param is.child.fn a function which takes two ids as arguments;
##' it is expected to return TRUE if the second argument is the parent of the
##' first, and FALSE otherwise. it defaults to \code{\link{is.child.ct}}
##' @return info
make.chain <- function(seed.id, survey.data, is.child.fn=is.child.ct) {

    key.var <- attr(survey.data, "key")
    keys <- survey.data[,key.var]

    is.child <- get.fn(is.child.fn)

    these.data <- survey.data[keys==seed.id,]

    child.ids <- keys[which(plyr::laply(keys, is.child, seed.id=seed.id))]

    ## if no children of this seed, finish
    if (length(child.ids)==0) {
        return (list(data=these.data,
                     children=NULL))
    }

    return(list(data=these.data,
                children=plyr::llply(child.ids,
                               make.chain,
                               survey.data=survey.data)))

}

#####################################################
##' get the height (maximum depth) of a chain
##'
##' get the height (maximum depth) of a chain
##'
##' @param chain the chain object
##' @return the maximum depth of the chain
max.depth <- function(chain) {
    if (is.null(chain$children)) {
        return(1)
    }

    return(1 + max(plyr::laply(chain$children,
                         max.depth)))
}

#####################################################
##' get the size of a chain
##'
##' count the total number of respondents in the chain
##' and return it
##'
##' @param chain the chain object
##' @return the number of respondents involved in
##' the chain
chain.size <- function(chain) {

    if (is.null(chain$children)) {
        return(1)
    }

    return(1 + sum(unlist(lapply(chain$children,
                                 chain.size))))
}

#####################################################
##' chain.vals
##'
##' get all of the values of the given variable
##' found among members of a chain
##'
##' @param chain the chain to get values from
##' @param qoi.var the name of the variable to
##' get from each member of the chain
##' @return a vector with all of the values of \code{qoi.var}
##' found in this chain. (currently, the order of the values
##' in the vector is not guaranteed)
chain.vals <- function(chain, qoi.var="uid") {

    if (is.null(chain$children)) {
        return(chain$data[,qoi.var])
    }

    return(c(chain$data[,qoi.var],
             unlist(lapply(chain$children,
                           chain.vals,
                           qoi.var=qoi.var))))
}

#####################################################
##' get a dataset from a chain
##'
##' take the data for each member of the given chain
##' and assemble it together in a dataset
##'
##' @param chain the chain to build a dataset from
##' @return a dataset comprised of all of the chain's
##' members' data put together. the order of the rows
##' in the datset is not specified.
chain.data <- function(chain) {
    if (is.null(chain$children)) {
        return(chain$data)
    }

    return(rbind(chain$data,
                 plyr::ldply(chain$children,
                       chain.data)))
}
