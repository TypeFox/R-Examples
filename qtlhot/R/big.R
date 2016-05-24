######################################################################
# big.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: big.phase0, big.phase1, big.phase2, big.phase3, do.big.phase2
######################################################################

big.phase0 <- function(dirpath, cross, trait.file, trait.matrix,
                       droptrait.names = NULL,
                       keeptrait.names = NULL,
                       lod.thrs,
                       sex = "Sex", trait.index, batch.effect = NULL, size.set = 250,
                       offset = 0, subset.sex = NULL,
                       verbose = TRUE)
{
  ## Set up cross object and trait.data objects.

  ## Cross object with sex, trait.index and batch.effect as only phenotypes.
  traits <- match(c(sex,trait.index, batch.effect), names(cross$pheno), nomatch = 0)
  cross$pheno <- cross$pheno[, traits, drop = FALSE]
  
  ## Subset to individuals with batch if used.
  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  ## Subset to individuals with sex if desired.
  if(!is.null(subset.sex))
    cross <- subset(cross, ind = (cross$pheno[[sex]] %in% subset.sex))
  
  ## Save cross object and other pertinent information.
  if(verbose)
    cat("Saving cross object in cross.RData\n")
  save(cross, lod.thrs, file = "cross.RData", compress = TRUE)

  ## Get trait values for selected traits.
  load(file.path(dirpath, trait.file))
  trait.names <- dimnames(get(trait.matrix))[[2]]
  all.traits <- seq(trait.names)
  if(!is.null(droptrait.names))
    all.traits[match(droptrait.names, trait.names, nomatch = 0)] <- 0
  if(!is.null(keeptrait.names))
    all.traits[-match(keeptrait.names, trait.names, nomatch = 0)] <- 0
  all.traits <- all.traits[all.traits > 0]
  n.traits <- length(all.traits)

  num.sets <- ceiling(n.traits / size.set)
  trait.nums <- seq(size.set)

  ## Save trait names, including those dropped, and size of sets.
  tmp <- paste("Trait", 0, "RData", sep = ".")
  if(verbose)
    cat("Saving trait metadata in", tmp, "\n")
  save(trait.file, trait.matrix, droptrait.names, size.set,
       trait.names, all.traits, num.sets,
       file = tmp, compress = TRUE)

  ## Cycle through all the phenotypes in sets of size size.set. Keeps R object smaller.
  for(pheno.set in seq(num.sets)) {
    traits <- all.traits[trait.nums[trait.nums <= n.traits]]
    
    trait.data <- get(trait.matrix)[, trait.names[traits], drop = FALSE]

    tmp <- paste("Trait", offset + pheno.set, "RData", sep = ".")
    if(verbose)
      cat("Saving", length(traits), "trait data as set", offset + pheno.set, "in", tmp, "\n")
    save(trait.data, offset, pheno.set, file = tmp, compress = TRUE)

    ## Get next size.set traits.
    trait.nums <- trait.nums + size.set
  }
}
#############################################################################################
big.phase1 <- function(dirpath = ".", cross.index = 0, params.file,
                       cross, lod.thrs, n.quant = 2000, n.perm = 1, n.split = 1,
                       batch.effect = NULL,
                       drop.lod = 1.5, lod.min = min(lod.thrs),
                       window = 5, seed = 0, big = TRUE, ...,
                       ## Following are loaded with Trait0.RData created in big.phase0. 
                       trait.index, trait.names, all.traits, size.set, seeds, lod.sums,
                       ##
                       verbose = FALSE)
{
  ## Want this to depend on parallel.phase1 as much as possible.
  ## Just include new stuff here.
  
  ## Phase 1 for big dataset.
  ## Sets up initial permutations.
  
  ## Get any parameters in file. These will overwrite passed arguments.
  eval(parse(file.path(dirpath, params.file)))
  
  ## Calculate genotype probabilities.
  cross <- calc.genoprob(cross, step=0.5, error.prob = 0.002,
                         map.function = "c-f", stepwidth = "max")
  
  ## Subset to individuals with batch if used.
  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  ## Load trait.
  load(file.path(dirpath, "Trait.0.RData"))

  ## Random numbers.
  if(n.perm > 1) {
    if(seed[1] > 0)
      set.seed(seed[1])

    n.ind <- nind(cross)
    seeds <- sample(c(98765:987654321), n.perm, replace = FALSE)
  }
  else
    seeds <- 0

  ## Big assumes n.split is n.perm and does one perm each run.
  n.split <- max(1, n.perm)
  
  save(cross, n.perm, seeds, n.quant, drop.lod, lod.thrs, lod.min, batch.effect,
       trait.index, trait.names, all.traits, size.set, ## From Trait.0.RData
       window, cross.index, n.split, big,
       file = "Phase1.RData", compress = TRUE)
  invisible()
}
#############################################################################################
big.phase2 <- function(dirpath = ".", index, ...,
                       ## Following are loaded with Phase1.RData created in big.phase1.
                       batch.effect = NULL, all.traits, trait.index, lod.min,
                       drop.lod, n.quant, lod.thrs, cross.index, seeds,
                       ##
                       addcovar = NULL, intcovar = NULL,
                       remove.files = TRUE, verbose = FALSE)
{
  ## Phase 2.
  ## Loop on sets of phenotypes, adding only what is needed.
  ## Loop on permutations internal to scanone.permutations.

  ## This loads cross object and many other things from phase1.
  load(file.path(dirpath, "Phase1.RData"))

  if(verbose)
    cat("compute covariates\n")
  if(!is.null(batch.effect))
    covars <- sexbatch.covar(cross, batch.effect, verbose = TRUE)
  else
    covars <- list(addcovar = addcovar, intcovar = intcovar)

  n.traits <- length(all.traits)
  perms <- NULL

  if(index <= 1) {
    ## Original data.
    do.big.phase2(dirpath, cross, covars, perms, 1, trait.index,
                  lod.min, drop.lod, remove.files, n.quant, lod.thrs, window, n.traits, cross.index,
                  ..., verbose)
  }
  else {
    ## Random permutation. Use preset seed if provided.
    seed <- seeds[index]
    if(seed > 0)
      set.seed(seed[1])
    if(verbose)
      cat("sample permutation", seed[1], "\n")
    n.ind <- nind(cross)
    perms <- sample(seq(n.ind), n.ind, replace = FALSE)
    cross$pheno <- cross$pheno[perms,]

    do.big.phase2(dirpath, cross, covars, perms, index, trait.index,
                  lod.min, drop.lod, remove.files, n.quant, lod.thrs, window, n.traits, cross.index,
                  ..., verbose)
  }
}
do.big.phase2 <- function(dirpath, cross, covars, perms, index, trait.index,
                          lod.min, drop.lod, remove.files, n.quant, lod.thrs, window, n.traits, cross.index,
                          ..., verbose,
                          ## Following supplied by Trait.*.RData created in big.phase0.
                          trait.data, offset, pheno.set)
{
  ## Cycle through all the phenotypes in sets of size size.set. Keeps R object smaller.
  ## Assume large trait matrix has been broken into Trait.i.RData, each with trait.data.

  filenames <- list.files(dirpath, "Trait.[1-9][0-9]*.RData")
  if(!length(filenames))
    parallel.error(5, 5, index)

  for(pheno.index in seq(length(filenames))) {
    if(verbose)
      cat(pheno.index, "\n")

    if(exists("trait.data"))
      rm("trait.data")
    ## Uses trait.data from Trait*.RData file.
    out <- with(file.path(dirpath, filenames[pheno.index]), {
      list(perm.cross = add.phenos(cross, trait.data, index = trait.index),
           pheno.col = find.pheno(perm.cross, dimnames(trait.data)[[2]]))
    })

    per.scan <- scanone(out$perm.cross, pheno.col = out$pheno.col, method = "hk", 
                        addcovar = covars$addcovar, intcovar = covars$intcovar, ...)

    per.scan.hl <- highlod(per.scan, lod.thr = lod.min, drop.lod = drop.lod,
                                 restrict.lod = TRUE)$highlod

    save(per.scan.hl, perms,
         file=file.path(dirpath, paste("per.scan",pheno.index, index,"RData",sep=".")))
  }

  chr.pos <- per.scan[,1:2]

  filenames <- list.files(dirpath, paste("per.scan", "[0-9][0-9]*", index, "RData",
                                         sep = "."))

  scan.hl <- cat.scanone(dirpath, filenames, chr.pos)

  if(remove.files)
    file.remove(dirpath, filenames)

  lod.sums <- as.list(max(scan.hl, lod.thrs, window), c("max.N", "max.N.window"))
  lod.sums$max.lod.quant <- quantile(scan.hl, n.quant = n.quant)
  
  save(scan.hl, lod.sums, cross.index, index,
       file = paste("perm", ".", cross.index, "_", index, ".RData", sep = ""))
}
#############################################################################################
big.phase3 <- function(dirpath = ".", index, cross.index, ...,
                       ## Following are loaded with Phase1.RData created in big.phase1.
                       n.perm, lod.thrs, n.quant, lod.sums,
                       ##
                       outfile = phase3name, verbose = FALSE)
{
  ## Phase 3. Merge back together.

  load(file.path(dirpath, "Phase1.RData"))

  filenames <- list.files(dirpath, paste("perm.", cross.index, "_[0-9][0-9]*.RData", sep = ""))
  n.file <- length(filenames)
  if(n.file != n.perm)
    warning(paste("Number of perm files", n.file, "does not match number of permutations", n.perm))
  ## First filename perm.*_1.RData is raw data. Drop it.
  tmp <- paste("perm.", cross.index, "_1.RData", sep = "")
  tmp <- match(tmp, filenames)
  if(is.na(tmp)) stop("perm run 1 with raw data not present")
  filenames <- filenames[-tmp]
  n.file <- length(filenames)

  if(n.file == 0)
    return(NULL)
  
  n.thrs <- length(lod.thrs)
  
  max.lod.quant <- matrix(NA, n.file, n.quant)
  max.N <- max.N.window <- matrix(NA, n.file, n.thrs)

  for(i.perm in seq(n.file)) {
    if(verbose)
      cat(i.perm, "\n")
    if(exists("lod.sums"))
      rm("lod.sums")
    out <- with(file.path(dirpath, filenames[i.perm]), {
      n.quant <- length(lod.sums$max.lod.quant)
      max.lod.quant <- lod.sums$max.lod.quant
      ## legacy adjustment
      if(length(lod.sums$max.N) == 2) {
        max.N <- lod.sums$max.N$max.N
        max.N.window <- lod.sums$max.N$max.N.win
      }
      else {
        max.N <- lod.sums$max.N
        max.N.window <- lod.sums$max.N.window
      }
      list(max.lod.quant = max.lod.quant, max.N = max.N, max.N.window = max.N.window)
    })
    max.lod.quant[i.perm, seq(length(out$max.lod.quant))] <- out$max.lod.quant
    max.N[i.perm,] <- out$max.N[i.perm,]
    max.N.window[i.perm,] <- out$max.N.window[i.perm,]
  }
  phase3name <- paste("Phase3", ifelse(cross.index > 0, cross.index, ""), ".RData", sep = "")

  save(max.lod.quant, max.N, max.N.window, file = outfile, compress = TRUE)
}
