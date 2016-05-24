# Update a factor based on the levels of another factor.
#
# @param .factor Factor
# @param reference.factor Factor whose levels to use.
# @param levels New levels
#
update.levels <- function(.factor,
                          reference.factor=NULL,
                          .levels=levels(reference.factor),
                          ordered=is.ordered(reference.factor))
{
    factor(.factor, levels=.levels, ordered=ordered)
}


# Retrieve the estimated tau for the given sample.
#
# @param dl de.lorean object
# @param sample.iter Which sample to use, defaults to best sample
#
tau.for.sample <- function(dl, sample.iter=dl$best.sample) {
    (
        dl$samples.l$tau
        %>% filter(sample.iter == iter)  # Filter correct iteration
        %>% arrange(c)  # Sort by cell
    )$tau
}


# The levels of the gene factor
#
# @param dl The de.lorean object.
#
gene.levels <- function(dl) levels=levels(dl$gene.meta$gene)


# The levels of the cell factor
#
# @param dl The de.lorean object.
#
cell.levels <- function(dl) levels=levels(dl$cell.meta$cell)


#' Melt an expression matrix.
#'
#' @param dl The de.lorean object.
#' @param expr Matrix of expression values.
#'
#' @export melt.expr
#'
melt.expr <- function(dl, expr=dl$expr) (
    expr
    %>% melt(varnames=c("gene", "cell"), value.name="x")
    %>% mutate(gene=factor(gene, levels=gene.levels(dl)),
               cell=factor(cell, levels=cell.levels(dl)))
)


# Cast an expression matrix.
#
# @param dl The de.lorean object
# @param expr.l Expression values in long format
#
cast.expr <- function(expr.l) expr.l %>% acast(gene ~ cell, value.var="x")


#' Analyse variance of expression between and within capture times.
#'
#' @param dl de.lorean object
#' @param adjust.cell.sizes Choose whether to adjust the expression values by
#'        the cell size estimates
#'
#' @export
#'
analyse.variance <- function(dl, adjust.cell.sizes) {
    within(dl, {
        #
        # First melt expression data into long format
        #
        expr.l <- melt(expr, varnames=c("gene", "cell"), value.name="x")
        expr.l$gene <- factor(expr.l$gene, levels=levels(gene.meta$gene))
        expr.l$cell <- factor(expr.l$cell, levels=levels(cell.meta$cell))
        #
        # Estimate a pseudo reference mean for each gene
        #
        gene.expr <- (expr.l
            %>% group_by(gene)
            %>% dplyr::summarise(x.mean=mean(x))
        )
        stopifnot(! is.na(gene.expr))
        # Estimate the cell size by the median of the expression
        # adjust by the gene's mean
        cell.expr <- (expr.l
            %>% left_join(gene.expr)
            %>% group_by(cell)
            %>% dplyr::summarise(S.hat=median(x - x.mean))
        )
        stopifnot(! is.na(cell.expr))
        expr.l <- expr.l %>% left_join(cell.expr)
        if (adjust.cell.sizes) {
            # Adjust the expression by the cell size estimates
            expr.l <- expr.l %>% mutate(x.hat=x - S.hat)
        } else {
            # Use the raw expression values
            expr.l <- expr.l %>% mutate(x.hat=x)
        }
        stopifnot(! is.na(expr.l))
        # Resummarise the adjusted expression data
        gene.expr <- (expr.l
            %>% group_by(gene)
            %>% dplyr::summarise(x.mean=mean(x),
                          x.sd=sd(x),
                          phi.hat=mean(x.hat),
                          x.hat.sd=sd(x.hat))
            %>% filter(! is.na(x.sd), ! is.na(x.hat.sd))
        )
        stopifnot(! is.na(gene.expr))
        stopifnot(nrow(gene.expr) > 0)  # Must have some rows left
        # Examine the variation within and between times for each gene
        gene.time.expr <- (
            expr.l
            %>% left_join(cell.meta)
            %>% group_by(gene, capture)
            %>% dplyr::summarise(x.mean=mean(x.hat),
                          x.var=mean((x.hat-mean(x.hat))**2),
                          num.capture=n(),
                          Mgc=mean(x),
                          Vgc=var(x))
            %>% filter(! is.na(x.var), ! is.na(Vgc))
        )
        stopifnot(! is.na(gene.time.expr))
        stopifnot(nrow(gene.time.expr) > 0)  # Must have some rows left
        # Decomposition of variance within and between time.
        gene.var <- (
            gene.time.expr
            %>% group_by(gene)
            %>% dplyr::summarise(Omega=weighted.mean(x.var, num.capture, na.rm=TRUE),
                          Psi=weighted.mean((x.mean-mean(x.mean))**2,
                                            num.capture, na.rm=TRUE))
            %>% filter(! is.na(Psi), Omega > 0 | Psi > 0)
        )
        stopifnot(! is.na(gene.var))
        stopifnot(nrow(gene.var) > 0)  # Must have some rows left
        # No longer needed
        rm(expr.l)
    })
}


#' Estimate hyperparameters for model using empirical Bayes.
#'
#' @param dl de.lorean object
#' @param sigma.tau Noise s.d. in temporal dimension, that is prior s.d. for tau
#' @param length.scale Length scale for stationary GP covariance function.
#'   Defaults to the range of the observed capture times.
#' @param model.name The model's name:
#'   \itemize{
#'     \item 'simplest-model': The simplest model (does not estimate the gene
#'       means).
#'     \item 'simple-model': Like 'simplest-model' but estimates the gene
#'       means.
#'     \item 'lowrank': Low rank approximation to the 'simplest-model'.
#'   }
#'
#' @export
#'
estimate.hyper <- function(
    dl,
    sigma.tau = .5,
    length.scale = NULL,
    model.name = 'simplest-model'
) {
  dl <- within(dl, {
    #
    # Remember options that depend on the model
    #
    opts$model.name <- model.name
    opts$estimate.phi <- switch(
        opts$model.name,
        "simple-model" = TRUE,
        "simplest-model" = FALSE,
        "lowrank" = FALSE,
        NA)
    opts$adjust.cell.sizes <- switch(
        opts$model.name,
        "simple-model" = TRUE,
        "simplest-model" = FALSE,
        "lowrank" = FALSE,
        NA)
    #
    # Set up temporal hyper-parameters
    #
    opts$sigma.tau <- sigma.tau
    time.range <- range(cell.meta$obstime)
    time.width <- time.range[2] - time.range[1]
    if (is.null(length.scale)) {
        opts$length.scale <- time.width
    } else {
        opts$length.scale <- length.scale
    }
    # message("Length scale: ", opts$length.scale)
  })
  dl <- analyse.variance(dl, dl$opts$adjust.cell.sizes)
  # Expected variance of samples from zero-mean Gaussian with covariance K.obs
  V.obs <- with(dl,
    expected.sample.var(
      cov.matern.32(
        cov.calc.dists(unique(dl$cell.meta$obstime)),
        dl$opts$length.scale)))
  within(dl, {
    gene.var <- gene.var %>% left_join(
      gene.time.expr
      %>% group_by(gene)
      %>% dplyr::summarise(omega.hat=mean(Vgc),
                           psi.hat=var(Mgc)/V.obs)
      %>% filter(! is.na(psi.hat),
                 ! is.na(omega.hat),
                 omega.hat > 0,
                 psi.hat > 0)
    )
    hyper <- list(
      mu_S=mean(cell.expr$S.hat),
      sigma_S=sd(cell.expr$S.hat),
      mu_phi=mean(gene.expr$phi.hat),
      sigma_phi=sd(gene.expr$phi.hat),
      mu_psi=mean(log(gene.var$psi.hat), na.rm=TRUE),
      sigma_psi=sd(log(gene.var$psi.hat), na.rm=TRUE),
      mu_omega=mean(log(gene.var$omega.hat), na.rm=TRUE),
      sigma_omega=sd(log(gene.var$omega.hat), na.rm=TRUE),
      sigma_tau=opts$sigma.tau,
      l=opts$length.scale
    )
    stopifnot(all(! sapply(hyper, is.na)))
  })
}


#' Filter genes
#'
#' @param dl de.lorean object
#' @param number Number to sample if filter function or genes not supplied.
#' @param genes The genes to keep.
#' @param .filter Function that gakes a list of genes as input and returns
#'     a vector of TRUE/FALSE
#'
#' @export
#'
filter.genes <- function(dl,
                         .filter=function(x) x %in% genes,
                         number=NULL,
                         genes=sample(rownames(dl$expr), number))
{
    within(dl, {
        expr <- expr[.filter(rownames(expr)),]
        message("Have ", nrow(expr), " genes after filtering")
    })
}


#' Filter cells
#'
#' @param dl de.lorean object
#' @param number Number to sample if filter function or cells not supplied.
#' @param cells The cells to keep.
#' @param .filter Function that gakes a list of cells as input and returns
#'     a vector of TRUE/FALSE
#'
#' @export
#'
filter.cells <- function(dl,
                         .filter=function(x) x %in% cells,
                         number=NULL,
                         cells=sample(colnames(dl$expr), number))
{
    within(dl, {
        expr <- expr[,.filter(colnames(expr))]
        message("Have ", ncol(expr), " cells after filtering")
    })
}




# Sample so many cells per capture time.
#
# @param dl de.lorean object
# @param number Number to sample from each capture time
#
sample.per.capture <- function(dl, cells.per.capture) {
    sample.at.most <- function(.df, number) {
        sample_n(.df, min(number, nrow(.df)))
    }
    sampled.cells <- (
        dl$cell.meta
        %>% filter(cell %in% colnames(dl$expr))
        %>% group_by(capture)
        # %>% do(sample_n(., min(cells.per.capture, length(cell))))
        %>% do(sample.at.most(., cells.per.capture))
    )
    filter.cells(dl, cells=sampled.cells$cell)
}


# Sample genes and cells
#
# @param dl de.lorean object
#
sample.genes.and.cells <- function(
    dl,
    max.cells = 0,
    max.genes = 0)
{
    within(dl, {
        opts$max.cells <- max.cells
        opts$max.genes <- max.genes
        if (opts$max.cells && ncol(expr) > opts$max.cells) {
            expr <- expr[,sample(ncol(expr), opts$max.cells)]
            # Remove genes that are not expressed in at least two cells
            num.cells.expr <- rowSums(! is.na(expr))
            expr <- expr[num.cells.expr > 1,]
        }
        if (opts$max.genes && nrow(expr) > opts$max.genes) {
            expr <- expr[sample(nrow(expr), opts$max.genes),]
        }
    })
}


#' Calculate inducing pseudotimes for sparse approximation
#'
#' @param dl de.lorean object
#' @param num.inducing Number of inducing points
#' @param period Period of expression patterns
#'
#' @export
#'
calc.inducing.pseudotimes <- function(dl, num.inducing, period = 0) with(dl, {
  if (period > 0) {
    seq(from=0,
        to=period * (1 - 1/num.inducing),
        length.out=num.inducing)
  } else {
    seq(from=time.range[1] - 3 * hyper$sigma_tau,
        to=time.range[2] + 3 * hyper$sigma_tau,
        length.out=num.inducing)
  }
})


#' Prepare for Stan
#'
#' @param dl de.lorean object
#' @param num.test Number of test points to consider
#' @param num.inducing Number of inducing points
#' @param period Period of expression patterns
#' @param hold.out Number genes to hold out for generalisation tests
#'
#' @export prepare.for.stan
#'
prepare.for.stan <- function(
    dl,
    num.test = 101,
    num.inducing = 30,  # M
    period = 0,
    hold.out = 0
) {
    within(dl, {
        opts$num.test <- num.test
        opts$period <- period
        opts$periodic <- opts$period > 0
        stopifnot(hold.out < nrow(expr))
        .G <- nrow(expr) - hold.out
        .C <- ncol(expr)
        #
        # Permute genes to make held out genes random
        expr <- expr[sample(.G+hold.out),]
        #
        # Calculate the map from gene indices to genes and their meta data
        gene.map <- (
            data.frame(g=1:(.G+hold.out),
                       gene=factor(rownames(expr),
                                   levels=levels(gene.meta$gene)))
            %>% mutate(is.held.out=g>.G)
            %>% left_join(gene.expr)
            %>% left_join(gene.var)
            %>% left_join(gene.meta))
        stopifnot(! is.na(gene.map[
            c("g", "gene", "x.mean", "x.sd",
            "phi.hat", "x.hat.sd", "Omega", "Psi",
            "psi.hat", "omega.hat")
        ]))
        #
        # Rename phi.hat if we are not estimating phi
        if (! opts$estimate.phi) {
            gene.map <- gene.map %>% dplyr::rename(phi=phi.hat)
        }
        #
        # Calculate the map from cell indices to genes and their meta data
        cell.map <- (data.frame(c=1:.C,
                                cell=factor(colnames(expr),
                                            levels=levels(cell.meta$cell)))
                    %>% left_join(cell.meta)
                    %>% left_join(cell.expr))
        stopifnot(! is.na(cell.map %>% dplyr::select(cell, capture, obstime)))
        #
        # Calculate the time points at which to make predictions
        if (opts$periodic) {
            test.input <- seq(0, opts$period, length.out=num.test)
        } else {
            test.input <- seq(
                time.range[1] - 2 * opts$sigma.tau,
                time.range[2] + 2 * opts$sigma.tau,
                length.out=num.test)
        }
        #
        # Gather all the data into one list
        stan.data <- c(
            # Hyper-parameters
            hyper,
            list(
                # Dimensions
                C=.C,
                G=.G,
                H=hold.out,
                M=num.inducing,
                # Data
                time=cell.map$obstime,
                expr=expr,
                # Inducing pseudotimes
                u=calc.inducing.pseudotimes(dl, num.inducing, period),
                # Held out parameters
                heldout_psi=filter(gene.map, g > .G)$psi.hat,
                heldout_omega=filter(gene.map, g > .G)$omega.hat,
                heldout_phi=filter(gene.map, g > .G)$phi.hat,
                # Generated quantities
                numtest=opts$num.test,
                testinput=test.input,
                # Periodic?
                periodic=opts$periodic,
                period=opts$period
            )
        )
        # If we're not estimating phi, add it to the data
        if (! opts$estimate.phi) {
            stan.data$phi <- gene.map$x.mean
        }
    })
}

#' Compile the model and cache the DSO to avoid unnecessary recompilation.
#'
#' @param dl de.lorean object
#'
#' @export
#'
compile.model <- function(dl) {
  stan.model.file <- system.file(file.path('Stan',
                                            sprintf('%s.stan',
                                                    dl$opts$model.name)),
                                  package='DeLorean',
                                  mustWork=TRUE)
  data.dir <- system.file('extdata', package='DeLorean')
  compiled.model.file <- paste(data.dir,
                               sprintf("%s.rds", dl$opts$model.name),
                               sep='/')
  within(dl, {
    if (file.exists(compiled.model.file)
        &&
        file.info(compiled.model.file)$mtime
            > file.info(stan.model.file)$mtime)
    {
      message("Loading pre-compiled model from ", compiled.model.file)
      compiled <- readRDS(compiled.model.file)
    } else {
      message("Compiling model")
      compiled <- rstan::stan(file=stan.model.file, chains=0,
                              data=stan.data)
      message("Saving compiled model to ", compiled.model.file)
      saveRDS(compiled, compiled.model.file)
    }
    # Try one iteration to check everything is OK
    # message("Trying iteration")
    fit <- rstan::stan(
      fit=compiled,
      data=stan.data,
      init=make.chain.init.fn(dl),
      warmup=1,
      iter=1,
      chains=1)
  })
}


# Define a function to initialise the chains
#
# @param dl de.lorean object
#
make.chain.init.fn <- function(dl) {
  function() {
    with(dl$stan.data, {
      # message("Creating initialisation")
      init <- list(
        S=dl$cell.map$S.hat,
        tau=rnorm(C, mean=time, sd=sigma_tau),
        psi=rlnorm(G, meanlog=mu_psi, sdlog=sigma_psi),
        omega=rlnorm(G, meanlog=mu_omega, sdlog=sigma_omega)
      )
      # If estimating phi, include it.
      if (dl$opts$estimate.phi) {
        init$phi <- rnorm(G, mean=mu_phi, sd=sigma_phi)
      }
      init
    })
  }
}


even.tau.spread <- function(dl) {
    with(dl$stan.data,
         seq(min(time) - sigma_tau,
             max(time) + sigma_tau,
             length=C))
}

# From http://stackoverflow.com/questions/11094822/numbers-in-geometric-progression
geom.series <- function(base, max) {
  base^(0:floor(log(max, base)))
}

#' Use seriation package to find good orderings
#'
#' @param dl de.lorean object
#' @param .methods The seriation methods to apply
#' @param scaled Whether to use the scaled and/or unscaled expression data
#' @param dim.red Dimension reduction methods to apply
#' @param dims Number of dimensions to reduce to
#' @param num.cores Number of cores to use in parallel
#' @param num.tau.to.keep How many initialisations to keep
#'
#' @export
#'
seriation.find.orderings <- function(
  dl,
  # .methods = c("ARSA", "TSP", "R2E", "HC", "GW", "OLO"),
  .methods = c("TSP", "R2E", "HC", "GW", "OLO"),
  scaled = c('scaled', 'unscaled'),
  dim.red = c('none', 'pca', 'kfa', 'ica', 'mds'),
  # dim.red = c('mds'),
  dims = geom.series(base=2, max=min(8, nrow(dl$expr)-1)),
  num.cores = default.num.cores(),
  num.tau.to.keep = default.num.cores())
{
  #
  # Calculate all combinations of parameters
  combinations <- expand.grid(method=.methods,
                              scaled=scaled,
                              dim.red=dim.red,
                              dims=dims,
                              stringsAsFactors=FALSE) %>%
    # Only keep one combination with no dimensionality reduction
    dplyr::filter(dim.red != 'none' | 1 == dims) %>%
    # Cannot do KFA with 1 dimension
    dplyr::filter(dim.red != 'kfa' | 1 != dims)
  get.expr <- function(scaled) switch(scaled,
    scaled = scale(t(dl$expr)),
    unscaled = t(dl$expr),
    stop('scaled must be "scaled" or "unscaled"')
  )
  get.expr.mem <- memoise::memoise(get.expr)
  get.red <- function(scaled, dim.red, dims=5) {
    expr <- get.expr(scaled)
    switch(dim.red,
      none = expr,
      pca = prcomp(expr)$x[,1:dims],
      ica = fastICA::fastICA(expr, n.comp=dims)$S,
      kfa = t(kernlab::kfa(t(expr), features=dims)@xmatrix),
      mds = cmdscale(dist(expr), k=dims),
      stop('dim.red must be "none", "ica", "kfa", "mds" or "pca"'))
  }
  get.red.mem <- memoise::memoise(get.red)
  get.dist <- function(scaled, dim.red, dims=5) {
    dist(get.red.mem(scaled, dim.red, dims))
  }
  get.dist.mem <- memoise::memoise(get.dist)
  result <- parallel::mclapply(
    1:nrow(combinations),
    function(i) with(combinations[i,], within(list(), {
      method.name <- stringr::str_c(method,
                                    ':', scaled,
                                    ':', dim.red,
                                    ':', dims)
      elapsed <- system.time(
        per <- seriation::seriate(get.dist.mem(scaled, dim.red, dims),
                                   method=method))
      ser.order <- seriation::get_order(per)
      ll <- dl$ordering.ll(ser.order)
      message(method.name, '; time=', elapsed[3], 's; LL=', ll)
    })),
    mc.cores=num.cores)
  result
}


#' Run a find good ordering method and append results to existing orderings
#'
#' @param dl de.lorean object
#' @param method Function that runs the method
#' @param ... Any other arguments for the method
#'
#' @export
#'
find.good.ordering <- function(dl, method, ...)
{
  dl <- create.ordering.ll.fn(dl)
  dl$order.inits <- c(
    dl$order.inits,
    lapply(
      method(dl, ...),
      function(i) {
        i$ser.order <- rev.order.if.better(dl, i$ser.order)
        i
      }))
  dl
}


#' Plot likelihoods of orderings against elapsed times taken
#' to generate them
#'
#' @param dl The DeLorean object
#'
#' @export
#'
orderings.plot <- function(dl) with(dl, {
  results.df <- data.frame(
    method=sapply(order.inits, function(r) r$method.name),
    elapsed=sapply(order.inits, function(r) r$elapsed[3]),
    ll=sapply(order.inits, function(r) r$ll))
  ggplot2::ggplot(results.df, aes(x=elapsed, y=ll, label=method)) +
    ggplot2::geom_text()
})

#' Use Magda's code to find good orderings
#'
#' @param dl de.lorean object
#' @param number_paths Number of paths for each starting point
#'
#' @export
#'
magda.find.orderings <- function(
    dl,
    number_paths = 5)
{
  #
  # Determine which cell indexes have either the highest or lowest
  # observation times, we will use these as starting points
  max.obs <- max(dl$cell.map$obstime)
  min.obs <- min(dl$cell.map$obstime)
  starting_points <-
    (1:dl$.C)[which(dl$cell.map$obstime %in% c(max.obs, min.obs))]
  #
  # Call Magda's function to generate paths from these starting points
  elapsed <- system.time(
    magda.paths <- CombfuncPaths(
      dl$expr,
      starting_points=starting_points,
      number_paths=number_paths))
  #
  # Convert Magda's paths into our format
  lapply(
    1:ncol(magda.paths),
    function(c) within(list(), {
      method.name <- stringr::str_c('Magda:', c)
      elapsed <- elapsed / ncol(magda.paths)
      ser.order <- magda.paths[,c]
      ll <- dl$ordering.ll(ser.order)
    }))
}

# Reverse ordering if it is better correlated with observed times
rev.order.if.better <- function(dl, ser.order) {
  rev.order <- rev(ser.order)
  if (cor(ser.order, dl$cell.map$obstime) < cor(rev.order, dl$cell.map$obstime)) {
    ser.order <- rev.order
  }
  ser.order
}

# Deduplicate orderings
#
deduplicate.orderings <- function(dl) {
  orderings <- sapply(dl$order.inits, function(i) i$ser.order)
  dupd <- duplicated(t(orderings))
  dl$order.inits <- dl$order.inits[!dupd]
  dl
}

#' Convert best orderings into initialisations
#'
#' @param dl The DeLorean object
#' @param num.to.keep The number to keep (defaults to default.num.cores())
#'
#' @export
#'
pseudotimes.from.orderings <- function(
  dl,
  num.to.keep=default.num.cores())
within(deduplicate.orderings(dl), {
  order.orderings <- order(sapply(order.inits, function(i) i$ll),
                           decreasing=TRUE)
  #
  # Make sure we don't try to keep too many (due to duplicates, etc...)
  actually.keep <- min(length(order.orderings), num.to.keep)
  if (actually.keep < num.to.keep) {
    warning("Don't have enough ", num.to.keep, " pseudotimes to keep,",
            " only have ", actually.keep)
  }
  best.orderings <- order.inits[order.orderings[1:actually.keep]]
  tau.inits <- lapply(
    best.orderings,
    function(O) {
      message('Using ordering ', O$method.name, '; LL=', O$ll)
      # Create an initialisation using the ordering
      init <- init.chain.sample.tau(dl)
      init$tau <- even.tau.spread(dl)[O$ser.order]
      init
    })
})


#' Find best order of the samples assuming some smooth GP prior on the
#' expression profiles over this ordering.
#'
#' @param dl de.lorean object
#' @param psi Temporal variation
#' @param omega Noise
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#' @param num.tau.to.try How many initialisations to try
#' @param num.tau.to.keep How many initialisations to keep
#' @param method Method to use "maximise" or "metropolis"
#' @param ... Extra arguments to method
#'
#' @export
#'
find.smooth.tau <- function(
    dl,
    psi = exp(dl$hyper$mu_psi),
    omega = exp(dl$hyper$mu_omega),
    num.cores = default.num.cores(),
    num.tau.to.try = num.cores,
    num.tau.to.keep = num.cores,
    method = "metropolis",
    ...
) {
  dl <- create.ordering.ll.fn(dl)
  log.likelihood <- dl$ordering.ll
  dl$tau.inits <- with(dl, {
    # Maximise the sum of the log marginal likelihoods
    ordering.search <- function(seed) {
      set.seed(seed)
      # Choose a starting point by random projection
      expr.centre <- t(scale(t(expr), center=T, scale=F))
      init.ordering <- order(rnorm(nrow(expr.centre)) %*% expr.centre)
      # init.ordering <- sample(stan.data$C)
      metropolis.fn <- function(ordering, log.likelihood, ...) {
        mh.run <- ordering.metropolis.hastings(
          ordering,
          log.likelihood,
          proposal.fn=ordering.random.block.move,
          ...)
        best.sample <- which.max(mh.run$log.likelihoods)
        #ordering.maximise(mh.run$chain[best.sample,], log.likelihood)
        mh.run$chain[best.sample,]
      }
      method.fn <- switch(method,
                          "maximise"=ordering.maximise,
                          "metropolis"=metropolis.fn,
                          NA)
      ordering <- method.fn(init.ordering, log.likelihood, ...)
      stopifnot(! is.null(ordering))
      # Reverse the ordering if it makes it correlate better with
      # the capture times
      capture.order <- order(stan.data$time)
      if (cor(capture.order, ordering) <
          cor(capture.order, rev(ordering)))
      {
        ordering <- rev(ordering)
      }
      ordering
    }
    # Choose seeds
    seeds <- sample.int(.Machine$integer.max, num.tau.to.try)
    # Run in parallel or not?
    if (num.cores > 1) {
      orderings <- parallel::mclapply(seeds,
                                      mc.cores=num.cores,
                                      ordering.search)
    } else {
      orderings <- lapply(seeds, ordering.search)
    }
    # Order the taus by the best orderings
    lls <- sapply(orderings, function(o) -log.likelihood(o))
    best.order <- order(lls)
    # Make the complete chain initialisation with the tau.
    lapply(orderings[best.order[1:num.tau.to.keep]],
           function(ordering) {
             init <- init.chain.sample.tau(dl)
             init$tau <- even.tau.spread(dl)[ordering.invert(ordering)]
             init
           })
  })
  dl
}


#' Test ordering Metropolis-Hastings sampler.
#'
#' @param dl de.lorean object
#' @param psi Temporal variation
#' @param omega Noise
#' @param num.cores Number of cores to run on.
#'          Defaults to getOption("DL.num.cores", max(parallel::detectCores()-1, 1))
#' @param iterations Number of iterations
#' @param thin Thin the samples
#'
test.mh <- function(
    dl,
    psi = mean(dl$gene.map$psi.hat),
    omega = mean(dl$gene.map$omega.hat),
    num.cores = getOption("DL.num.cores", max(parallel::detectCores() - 1, 1)),
    iterations = 1000,
    thin = 15
) {
    dl <- create.ordering.ll.fn(dl)
    log.likelihood <- dl$ordering.ll
    dl$tau.inits <- with(dl, {
        # Maximise the sum of the log marginal likelihoods
        ordering.search <- function(seed) {
            set.seed(seed)
            # Choose a random starting point
            init.ordering <- sample(stan.data$C)
            mh.run <- ordering.metropolis.hastings(
                init.ordering,
                log.likelihood,
                proposal.fn=ordering.random.block.move,
                iterations=iterations,
                thin=thin)
        }
        # Choose seeds
        seeds <- sample.int(.Machine$integer.max, num.cores)
        # Run in parallel or not?
        if (num.cores > 1) {
            orderings <- parallel::mclapply(seeds,
                                  mc.cores=num.cores,
                                  ordering.search)
        } else {
            orderings <- lapply(seeds, ordering.search)
        }
        orderings
    })
}

# The covariance function for the DeLorean object.
cov.fn.for <- function(dl) {
    cov.matern.32
}

# Choose an initialisation by sampling tau from the prior.
#
# @param dl de.lorean object
#
init.chain.sample.tau <- function(dl) {
    with(dl$stan.data, {
        init <- list(
            alpha=dl$cell.map$alpha.hat,
            beta=rep(0, G),
            S=dl$cell.map$S.hat,
            tau=rnorm(C, time, sd=sigma_tau),
            phi=dl$gene.map$phi.hat[1:G],
            psi=dl$gene.map$psi.hat[1:G],
            omega=dl$gene.map$omega.hat[1:G]
        )
        init$tauoffsets <- init$tau - time
        # If not estimating phi, don't include it.
        if (! dl$opts$estimate.phi) {
            init$phi <- NULL
        }
        init
    })
}


#' Find best tau to initialise chains with by sampling tau from the prior
#' and using empirical Bayes parameter estimates for the other parameters.
#'
#' @param dl de.lorean object
#' @param num.tau.candidates How many candidates to examine. Defaults to 6000.
#' @param num.tau.to.keep How many candidates to keep. Defaults to num.cores.
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#'
#' @export
#'
find.best.tau <- function(
  dl,
  num.tau.candidates = 6000,
  num.tau.to.keep = num.cores,
  num.cores = default.num.cores())
{
  within(dl, {
    # Define a function that calculates log probability for
    # random seeded tau
    try.tau.init <- function(i) {
      set.seed(i)
      pars <- init.chain.sample.tau(dl)
      lp <- rstan::log_prob(fit, rstan::unconstrain_pars(fit, pars))
      list(lp=lp, tau=pars$tau)
    }
    # Choose tau several times and calculate log probability
    if (num.cores > 1) {
      tau.inits <- parallel::mclapply(1:num.tau.candidates,
                                      mc.cores=num.cores,
                                      try.tau.init)
    } else {
      tau.inits <- lapply(1:num.tau.candidates, try.tau.init)
    }
    # qplot(sapply(tau.inits, function(init) init$lp))
    # Which tau gave highest log probability?
    tau.inits.order <- order(sapply(tau.inits, function(init) -init$lp))
    # Just keep so many best tau inits
    tau.inits <- tau.inits[tau.inits.order[1:num.tau.to.keep]]
    rm(tau.inits.order, try.tau.init)
  })
}


#' Perform all the steps necessary to fit the model.
#' - prepare the data
#' - compile the model
#' - find suitable initialisations
#' - fit the model using the specified method (sampling or variational Bayes)
#' - process the posterior.
#'
#' @param dl de.lorean object
#' @param method Fitting method:
#'   \itemize{
#'     \item 'sample': Use a Stan sampler.
#'       See \code{\link{fit.model.sample}}.
#'     \item 'vb': Use Stan ADVI variational Bayes algorithm.
#'       See \code{\link{fit.model.vb}}.
#'   }
#' @param ... Extra arguments for fitting method
#'
#' @export
#'
fit.dl <- function(
    dl,
    method = 'sample',
    ...)
{
  dl <- prepare.for.stan(dl)
  dl <- compile.model(dl)
  dl <- find.good.ordering(dl, seriation.find.orderings)
  dl <- pseudotimes.from.orderings(dl)
  dl <- fit.model(dl, method=method, ...)
  dl <- process.posterior(dl)
  dl <- analyse.noise.levels(dl)
}


#' Fit the model using specified method (sampling or variational Bayes).
#'
#' @param dl de.lorean object
#' @param method Fitting method:
#'   \itemize{
#'     \item 'sample': Use a Stan sampler.
#'       See \code{\link{fit.model.sample}}.
#'     \item 'vb': Use Stan ADVI variational Bayes algorithm.
#'       See \code{\link{fit.model.vb}}.
#'   }
#' @param ... Extra arguments for method
#'
#' @export
#'
fit.model <- function(
    dl,
    method = 'sample',
    ...)
{
    switch(
        method,
        "sample" = fit.model.sample(dl, ...),
        "vb" = fit.model.vb(dl, ...),
        stop('Unknown method'))
}


#' Fit the model using Stan sampler
#'
#' @param dl de.lorean object
#' @param num.cores Number of cores to run on.
#'   Defaults to getOption("DL.num.cores", max(parallel::detectCores()-1, 1))
#' @param chains Number of chains to run on each core
#' @param thin How many samples to generate before retaining one
#' @param ... Extra arguments for rstan::stan() sampling call
#'
#' @export
#'
fit.model.sample <- function(
    dl,
    num.cores = getOption("DL.num.cores", max(parallel::detectCores() - 1, 1)),
    chains = 1,
    thin = 50,
    ...)
{
  init.chain.good.tau <- make.init.fn(dl)
  # Run the chains in parallel
  sflist <- parallel::mclapply(
    1:num.cores,
    mc.cores=num.cores,
    function(i)
      rstan::stan(
        fit=dl$fit,
        data=dl$stan.data,
        thin=thin,
        init=init.chain.good.tau,
        seed=i,
        chains=chains,
        chain_id=i,
        refresh=-1,
        ...))
  dl$fit <- rstan::sflist2stanfit(sflist)
  dl$compiled <- NULL  # Delete large unneeded object
  return(dl)
}


#' Returns a function that constructs parameter settings with good tau.
#'
#' @param dl de.lorean object
#'
make.init.fn <- function(dl) {
  function(chain_id) {
    stopifnot(chain_id <= length(dl$tau.inits))
    #
    # Create random parameters
    pars <- make.chain.init.fn(dl)()
    #
    # Replace tau with good tau
    pars$tau <- dl$tau.inits[[chain_id]]$tau
    pars
  }
}


#' Average across a parameters samples.
#'
#' @param s An array of any dimension in which the first dimensions
#'   indexes the samples
#'
avg.par.samples <- function(s) {
  n.dims <- length(dim(s))
  if (n.dims > 1) {
    apply(s, 2:n.dims, mean)
  } else {
    mean(s)
  }
}

#' Get posterior mean of samples
#'
#' @param extract A named list of samples
#'
#' @export
#'
get.posterior.mean <- function(extract) lapply(extract, avg.par.samples)


#' Fit the model using Stan variational Bayes
#'
#' @param dl de.lorean object
#' @param num.cores Number of cores to run on. Defaults to default.num.cores()
#' @param num.inits Number initialisations to try. Defaults to num.cores
#' @param init.idx Which initialisation to use if only using one
#' @param ... Extra arguments for rstan::vb()
#'
#' @export
#'
fit.model.vb <- function(
    dl,
    num.cores = default.num.cores(),
    num.inits = num.cores,
    init.idx = 1,
    ...)
{
  init.chain.good.tau <- make.init.fn(dl)
  if (num.cores > 1) {
    #
    # Run variational Bayes in parallel
    sflist <- parallel::mclapply(
      1:num.inits,
      mc.cores=num.cores,
      # mc.cores=1,
      function(i) within(list(), {
        fit <- rstan::vb(
          rstan::get_stanmodel(dl$fit),
          data=dl$stan.data,
          seed=i,
          init=init.chain.good.tau(i),
          ...)
        pars <- get.posterior.mean(rstan::extract(fit))
        upars <- rstan::unconstrain_pars(fit, pars)
        lp <- rstan::log_prob(fit, upars)}))
    #
    # Only keep results that worked
    sflist <- Filter(function(x) ! is.null(x$lp), sflist)
    n.worked <- length(sflist)
    if (0 == n.worked) {
      stop('No VB fits worked.')
    } else if (n.worked < num.inits) {
      warning('Only ', n.worked, '/', num.inits, ' VB fits succeeded.')
    }
    #
    # Get the log likelihoods as a vector
    dl$vb.lls <- sapply(sflist, function(sf) sf$lp)
    #
    # Calculate which run had best lp for posterior mean parameters
    best.idx <- which.max(dl$vb.lls)
    #
    # Save the estimated tau for analysis
    dl$vb.tau <- sapply(sflist, function(sf) sf$pars$tau)
    #
    # Use those results
    dl$fit <- sflist[[best.idx]]$fit
  } else {
    # Run single variational Bayes
    dl$fit <- rstan::vb(
      rstan::get_stanmodel(dl$fit),
      data=dl$stan.data,
      seed=init.idx,
      init=init.chain.good.tau(init.idx),
      ...)
  }
  return(dl)
}


#' Analyse the samples and gather the convergence statistics. Note this
#' only makes sense if a sampling method was used to fit the model as
#' opposed to variational Bayes.
#'
#' @param dl de.lorean object
#'
#' @export
#'
examine.convergence <- function(dl) {
    within(dl, {
        pars <- c("tau", "psi", "S", "omega")
        if (opts$estimate.phi) {
            pars <- c(pars, "phi")
        }
        summ <- rstan::monitor(fit,
                        print=FALSE,
                        pars=pars)
        ignore.names <- stringr::str_detect(rownames(summ),
                                   "^(predictedvar|predictedmean)")
        rhat.sorted <- sort(summ[! ignore.names, "Rhat"])
        rhat.df <- data.frame(
            rhat=rhat.sorted,
            param=names(rhat.sorted),
            parameter=stringr::str_match(names(rhat.sorted), "^[[:alpha:]]+"))
        rm(summ)
        rm(pars)
    })
}


# The dimensions of the model parameters
#
# @param dl de.lorean object
#
model.parameter.dimensions <- function(dl) {
    sample.dims <- list(
        lp__=c(),
        S=c("c"),
        tau=c("c"),
        phi=c("g"),
        psi=c("g"),
        omega=c("g"),
        predictedmean=c("g", "t"),
        predictedvar=c("g", "t"),
        logmarglike=c("g")
    )
    if (! dl$opts$estimate.phi) {
        sample.dims$phi <- NULL
        sample.dims$S <- NULL
    }
    sample.dims
}


# Sample melter
#
# @param dl de.lorean object
#
sample.melter <- function(dl, include.iter=TRUE) {
    function(sample.list, sample.dims) {
        melt.var <- function(param) {
            # message(param)
            if (include.iter) {
                varnames <- c("iter", sample.dims[[param]])
            } else {
                varnames <- sample.dims[[param]]
            }
            melt(sample.list[[param]], varnames, value.name=param)
        }
        sapply(names(sample.dims), melt.var)
    }
}


# Join extra data to tau samples.
#
# @param dl de.lorean object
#
join.tau.samples <- function(dl, tau.samples) {
    with(dl,
         tau.samples
             %>% left_join(cell.map)
             %>% mutate(tau.offset=tau-obstime))
}


#' Process the posterior, that is extract and reformat the samples from
#' Stan. We also determine which sample has the highest likelihood, this
#' is labelled as the 'best' sample.
#'
#' @param dl de.lorean object
#'
#' @export
#'
process.posterior <- function(dl) {
    within(dl, {
        # Define a function to melt samples into a long format
        samples.l <- sample.melter(dl)(rstan::extract(dl$fit, permuted=TRUE),
                                       model.parameter.dimensions(dl))
        best.sample <- which.max(samples.l$lp__$lp__)
        if (TRUE %in% samples.l$logmarglike$is.held.out) {
            mean.held.out.marg.ll <- mean(
                (samples.l$logmarglike
                %>% left_join(gene.map)
                %>% filter(is.held.out))$logmarglike)
            message('Mean held out marginal log likelihood per cell: ',
                    mean.held.out.marg.ll / stan.data$C)
        }
        # Include meta data in tau samples
        samples.l$tau <- join.tau.samples(dl, samples.l$tau)
    })
}


# Optimise the log posterior starting at a particular sample or from some
# other set of parameters.
#
# @param dl de.lorean object
# @param sample.iter Sample to optimise (defaults to best sample).
#
optimise.sample <- function(
    dl,
    parameters=sample.parameters(dl, sample.iter=sample.iter),
    sample.iter=dl$best.sample,
    ...)
{
    with(dl, {
        optimised <- rstan::optimizing(dl$fit@stanmodel,
                                data=dl$stan.data,
                                init=parameters,
                                as_vector=FALSE,
                                ...)
        dims <- model.parameter.dimensions(dl)
        # Don't melt lp__ value
        dims$lp__ <- NULL
        samples.l <- sample.melter(dl, include.iter=FALSE)(optimised$par, dims)
        # Include meta data in tau samples
        samples.l$tau <- join.tau.samples(dl, samples.l$tau)
        # Include lp__ value
        samples.l$lp__ <- data.frame(lp__=optimised$value)
        samples.l
    })
}


# Bind a sample.
#
# @param dl de.lorean object
# @param samples Samples to bind to existing samples.
# @param sample.iter Iteration (defaults to -1).
#
bind.sample <- function(dl, samples, sample.iter=-1)
{
    within(
        dl,
        for (param in names(samples.l)) {
            samples.l[[param]] <- rbind(samples[[param]]
                                            %>% mutate(iter=sample.iter),
                                        samples.l[[param]])
        })
}


#' Optimise the best sample and update the best.sample index.
#'
#' @param dl de.lorean object
#' @param sample.to.opt Sample to optimise
#' @param new.best.sample Update to best sample index
#'
#' @export
#'
optimise.best.sample <- function(
    dl,
    sample.to.opt=dl$best.sample,
    new.best.sample=-1)
{
    dl$sample.optimised <- sample.to.opt
    dl <- bind.sample(dl, optimise.sample(dl, sample.iter=sample.to.opt))
    dl$sample.old.best <- dl$best.sample
    dl$best.sample <- new.best.sample
    dl
}


#' Analyse noise levels and assess which genes have the greatest
#' ratio of temporal variance to noise. This are labelled as the
#' 'gene.high.psi' genes.
#'
#' @param dl de.lorean object
#' @param num.high.psi How many genes with high variance to examine
#'
#' @export
#'
analyse.noise.levels <- function(dl, num.high.psi=25) {
    within(dl, {
        noise.levels <- (
            with(samples.l, left_join(psi, omega))
            %>% left_join(gene.map))
        # Summarise by gene
        gene.noise.levels <- (
            noise.levels
            %>% group_by(g)
            %>% dplyr::summarise(omega=mean(omega), psi=mean(psi))
            %>% left_join(gene.map)
            %>% arrange(-psi/omega))
        genes.high.psi <- head(gene.noise.levels$gene, num.high.psi)
    })
}


# Get the sampled parameter for the gene
#
# @param dl de.lorean object
# @param gene.idx Gene index
# @param param Parameter
# @param sample.iter Iteration to use (defaults to best.sample)
#
sampled.gene.param <- function(dl,
                               gene.idx,
                               param,
                               sample.iter=dl$best.sample) {
    with(dl, {
        filter(samples.l[[param]], gene.idx == g, sample.iter == iter)[[param]]
    })
}


#' Make predictions
#'
#' @param dl de.lorean object
#'
#' @export
#'
make.predictions <- function(dl) {
    within(dl, {
        predictions <- with(samples.l,
                            predictedmean
                            %>% left_join(predictedvar)
                            # %>% left_join(S)
                            %>% mutate(tau=test.input[t]))
        if (opts$estimate.phi) {
            predictions <- predictions %>% left_join(samples.l$phi)
        }
    })
}


#' Fit held out genes
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#' @param sample.iter The sample to use to fit with
#'
fit.held.out <- function(
    dl,
    expr.held.out,
    sample.iter=dl$best.sample)
{
    with(dl, {
        if (opts$adjust.cell.sizes) {
            cell.posterior <- samples.l$S %>% filter(sample.iter == iter)
            expr.held.out <- t(t(expr.held.out) + cell.posterior$S)
        }
        #' Calculate covariance over pseudotimes and capture times
        calc.K <- functional::Curry(cov.calc.gene,
                        dl,
                        include.test=F,
                        psi=exp(stan.data$mu_psi),
                        omega=exp(stan.data$mu_omega))
        tau <- tau.for.sample(dl, sample.iter=sample.iter)
        obstime <- cell.map$obstime
        K.tau <- calc.K(tau=tau)
        K.capture <- calc.K(tau=obstime)
        #' Evaluate the held out gene under the GP model using pseudotimes
        #' and a model without.
        #'
        calc.gp.marginals <- function(expr) {
            c(
                gp.log.marg.like(expr, K.tau),
                gp.log.marg.like(expr, K.capture))
        }
        fit.model <- function(expr, model=loess) {
            list(tau=model(expr~s(tau)))
        }
        apply(expr.held.out, 1, calc.gp.marginals)
        # list(
            # gp.marginals=sapply(held.out.genes, calc.gp.marginals),
            # loess=lapply(held.out.genes, functional::Curry(fit.model, model=gam)))
            # gam=lapply(held.out.genes, functional::Curry(fit.model, model=gam)))
    })
}


# Parameter values for sample
#
# @param dl de.lorean object
# @param sample.iter The sample we want the parameters for.
#
sample.parameters <- function(dl,
                              sample.iter=dl$best.sample,
                              param.names=names(dl$samples.l))
{
    parameters <- lapply(param.names,
                         function(param)
                             filter(dl$samples.l[[param]],
                                    iter == sample.iter)[[param]])
    names(parameters) <- param.names
    parameters
}


#' Test fit for log normal and gamma
#'
#' @param vars Data to fit
#'
#' @export
#'
test.fit <- function(vars) {
    fit.gamma <- MASS::fitdistr(vars, 'gamma')
    fit.lognormal <- MASS::fitdistr(vars, 'lognormal')
    gp <- (
        ggplot(data.frame(V=vars), aes(x=V))
        + geom_density()
        + stat_function(fun=functional::Curry(dgamma,
                                shape=fit.gamma$estimate['shape'],
                                rate=fit.gamma$estimate['rate']),
                        linetype='dashed')
        + stat_function(fun=functional::Curry(dlnorm,
                                meanlog=fit.lognormal$estimate['meanlog'],
                                sdlog=fit.lognormal$estimate['sdlog']),
                        linetype='dotted')
    )
    list(gamma=fit.gamma, lognormal=fit.lognormal, gp=gp)
}


# The samples
#
# @param dl de.lorean object
#
sample.iters <- function(dl) dl$samples.l$lp__$iter
