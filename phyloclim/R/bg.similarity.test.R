bg.similarity.test <- function(p, env, n = 99, conf.level = .95, app, dir){
  
  # checks and definitions
  # ----------------------
  names(p) <- c("species", "long", "lat")
  layer.names <- names(env)
  species <- sort(unique(levels(p[, 1])[p[, 1]]))
	
  ## sample background points from env
  ## ---------------------------------
  bg <- sampleRandom(env, size = 9999, na.rm = TRUE, sp = TRUE)
  bg <- data.frame("background", coordinates(bg), slot(bg, "data"))
  
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
  
  ## sample n(spec1) and n(spec2) points from background
  ## ---------------------------------------------------
  nb.occ <- table(p[, 1])[species] # corresponds to 'o' (p.2872)
  spec.vect <- sort(as.character(levels(p[, 1])[p[, 1]]))
  random.presence <- function(i, x, nbo, p, name){
    name <- paste(name, i, sep = "_")
    s <- rbind(sampleRandom(x, size = nbo[1], na.rm = TRUE, sp = TRUE),
               sampleRandom(x, size = nbo[2], na.rm = TRUE, sp = TRUE))
    s <- data.frame(name, coordinates(s), slot(s, "data"))
    colnames(s)[1:3] <- c("species", "long", "lat")
    return(s)
  }
  rp <- lapply(1:n, FUN = random.presence, x = env, nbo = nb.occ, name = spec.vect)
  rp <- do.call(rbind, rp)
	
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
  write.table(rbind(p, rp), paste(DIR, "samples.csv", sep = "/"),
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
	
  # calculate D and I for actual models
  # -----------------------------------
  fns <- paste(ODIR, species, "_proj.asc", sep = "")
  x <- read.asciigrid(fns[1])
  y <- read.asciigrid(fns[2])
  di <- di.enm(x = x, y = y)
  
  # calculate D and I for null distributions 
  # ----------------------------------------
  di.x.randomY <- sapply(X = 1:n, FUN = di.enm, x = x, y = fns[2])
  di.x.randomY <- t(di.x.randomY)
  di.y.randomX <- sapply(X = 1:n, FUN = di.enm, x = fns[1], y = y)
  di.y.randomX <- t(di.y.randomX)
  
  # CIs for null distributions
  # ------------------------------
  conf.limits <- c((1 - conf.level) / 2, 1- (1 - conf.level) / 2)
  ci.x.randomY <- apply(di.x.randomY, 2, quantile, probs = conf.limits)
  ci.y.randomX <- apply(di.y.randomX, 2, quantile, probs = conf.limits)
  
#   h0.x.randomY <- ci.x.randomY[1, ] < di & di < ci.x.randomY[2, ]
#   h0.y.randomX <- ci.y.randomX[1, ] < di & di < ci.y.randomX[2, ]
  
  # The null hypothesis that measured niche overlap between species is explained
  # by regional similarities or differences in available habitat (and not by
  # niche conservatism), is rejected if the actual similarity between two
  # species falls outside of the 95% confidence limits of the null distribution.

  # remove MAXENT output:
  # ---------------------
  if (DIR == "R.phyloclim.temp") unlink(DIR, recursive = TRUE)
	
	# create output object:
	# ---------------------
  out <- list(
    method = "background similarity test",
    species = species,
    null = paste("niche models are either more similar \n", 
                 paste(rep(" ", 24), collapse = ""), 
                 "or more different than expected by chance", sep = ""),
    statistic = di,
    ci.x.randomY = ci.x.randomY,
    ci.y.randomX = ci.y.randomX,
    nd.x.randomY = di.x.randomY,
    nd.y.randomX = di.y.randomX
    )
  class(out) <- "ntest"
  out
}