niche.equivalency.test <- function(p, env, n = 99, app, dir){
	
  # checks and definitions
  # ----------------------
  layer.names <- names(env)
  species <- sort(unique(levels(p[, 1])[p[, 1]]))
  
  # append covariate data to presence points
  # ----------------------------------------
  attributes <- extract(x = env, y = SpatialPoints(p[, 2:3]))
  p <- data.frame(p, attributes)
  NAs <- which(is.na(p), arr.ind = TRUE)
  NAs <- unique(NAs[, 1])
  if (length(NAs) > 0){
    p <- p[-NAs, ]
    warning(length(NAs), " presence points with missing environmental data removed")
  }
		
	# create permutated presence records pp
	# -------------------------------------
  permutate <- function(n, presence){
    presence[, 1] <- presence[sample(nrow(presence)), 1]
    presence[, 1] <- paste(presence[, 1], n, sep = "_")
    return(presence)
  }
	pp <- lapply(1:n, permutate, presence = p)
  pp <- do.call(rbind, pp)
  
  ## sample background points from env
  ## ---------------------------------
  bg <- sampleRandom(env, size = 9999, na.rm = TRUE, sp = TRUE)
  bg <- data.frame("background", coordinates(bg), slot(bg, "data"))
  
	# save input files:
	# -----------------
  if (missing(dir)) {
    DIR <- "R.phyloclim.temp"
  }
  else {
    DIR <- dir
  }
  if (file.exists(DIR))
    unlink(DIR, recursive = TRUE)
  dir.create(DIR)
  dir.create(ODIR <- paste(DIR, "out/", sep = "/"))
  dir.create(PDIR <- paste(DIR, "proj/", sep = "/"))
  
  write.table(bg, paste(DIR, "background.csv", sep = "/"), 
              row.names = FALSE, col.names = TRUE, sep = ",")
  write.table(rbind(p, pp), paste(DIR, "samples.csv", sep = "/"),
              row.names = FALSE, col.names = TRUE, sep = ",")
  fn <- paste(PDIR, layer.names, ".asc", sep = "")
  env <- unstack(env)
  for (i in seq_along(fn)){
    writeRaster(x = env[[i]], filename = fn[i], format = "ascii", 
                overwrite = TRUE, NAflag = -9999)
  }
  
	# call MAXENT:
	# ------------
  togglelayertype <- ifelse(length(grep("cat_", layer.names)) > 0, "-t cat_", "")
	CALL <- paste("java -jar", app , 		
    "-e ", paste(DIR, "background.csv", sep = "/"),
		"-s ", paste(DIR, "samples.csv", sep = "/"),
		"-j ", PDIR, 
		"-o ", ODIR,
    togglelayertype,
		"-r removeduplicates nopictures autorun")
	system(CALL, wait = TRUE)
	
	# calculate D and I for original parameters and null distribution 
	# ---------------------------------------------------------------
	fns <- paste(ODIR, species, "_proj.asc", sep = "")
  x <- read.asciigrid(fns[1])
  y <- read.asciigrid(fns[2])
  di <- di.enm(x = x, y = y)
  
  di.random <- sapply(X = 1:n, FUN = di.enm, x = fns[1], y = fns[2])
  di.random <- t(di.random)
	
	# assess significance:
	# --------------------
  # The observed values of I(pX,pY) and D(pX,pY) are compared to the percentiles 
  # of these null distributions in a one-tailed test to evaluate the hypothesis 
  # that the niche models for X and Y are not statistically significantly different.
  m <- colMeans(di.random)
  s <- apply(di.random, 2, sd)
  p.D <- pnorm(di[1], m[1], s[1], lower.tail = TRUE)
  p.I <- pnorm(di[2], m[2], s[2], lower.tail = TRUE)
  
	# remove MAXENT output:
	# ---------------------
  if (DIR == "R.phyloclim.temp") unlink(DIR, recursive = TRUE)
	
	# create output object:
	# ---------------------
	out <- list(
		method = "niche identity test",
		species = species,
		null = "niche models are not identical/equivalent",
    statistic = di,
    p.value = c(p.D, p.I),
		null.distribution = di.random
	)
  class(out) <- "ntest"
  out
}