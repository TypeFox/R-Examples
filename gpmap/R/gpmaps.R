## This file contains functions for creating, printing, summarizing and 
## plotting gpmap objects.

generate_gpmap <- function(y, locinames = NULL, allelenames = NULL, mapnames = NULL) {
  # Create a gpmap object 

  y <- as.matrix(y)

  #find number of loci and check number of genotypic values
  ngenotypes <- dim(y)[1]  # number of genotypic values
  nloci <- round(log(ngenotypes, 3))  # 3^nloci genotypes for biallelic loci
  if (3^nloci!=ngenotypes) { 
    stop(paste('Wrong number of genotypic values ', ngenotypes, ', should be 3^nloci',sep=""))
  }

  # names of gpmap
  npheno <- dim(y)[2]
  if (is.null(mapnames) & !is.null(colnames(y))) {
    mapnames <- colnames(y)
  } else if (is.null(mapnames)) { 
    mapnames <- paste("GPmap_", 1:npheno, sep = "") 
  }

  #create gpmap object
  gp <- NULL
  genotype <- enumerate_genotypes(nloci, locinames, allelenames)
  gp$nloci <- as.integer(nloci)
  gp$locinames <- colnames(genotype)
  gp$nmaps <- npheno
  gp$mapnames <- mapnames
  gp$genotype <- genotype
  gp$values <- y
  colnames(gp$values) <- mapnames
  class(gp) <- "gpmap"
  return(gp)
}

plot.gpmap <- function(x, show=1, decomposed=FALSE, ...) {
  ##TODO: pick out 
  
  gpmap <- x
  if(x$nmaps>1){
    if(show==1){print('Plotting gpmap number 1, change with argument "show".')}
    gpmap$nmaps <- 1
    gpmap$mapnames <- gpmap$mapnames[show]
    gpmap$values <- gpmap$values[,show]
    dim(gpmap$values) <- c(3^gpmap$nloci,1)
    if(decomposed){
      gpmap$values.mono <- gpmap$values.mono[,show]
      dim(gpmap$values.mono) <- c(3^gpmap$nloci,1)
      }
    }
  
  if(gpmap$nloci>3) { 
      print('No plot method for gpmap with more than 3 loci')
      return()
      }
  
  if (!decomposed) { 
    return(eval(parse(text=paste('plot',gpmap$nloci,'_orig(gpmap)',sep=""))))
  } else {
    return(eval(parse(text=paste('plot',gpmap$nloci,'_dec(gpmap)',sep=""))))
  }
}
  
plot1_orig <- function(gpmap) {
  # Plot gpmap with 1 locus
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp  <- cbind(gpmap$genotype, gpmap$values[, i])
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i])
    gp <- data.frame(gp)
    p <- ggplot(data = gp, aes_string(x = gpmap$locinames[1],  
                               y = gpmap$mapnames[i],
                               group = gpmap$locinames[1]))
    p <- p + geom_line(size = 2) + labs(x = gpmap$locinames[1],
                                        y = "Genotypic value")
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Original GPmap)",sep="\n"))
    p
  }
  return(figs)
}

plot2_orig <- function(gpmap) {
  # Plot gpmap object with 2 loci
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp <- cbind(gpmap$genotype, gpmap$values[, i]) 
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i])
    p <- ggplot(data = gp, aes_string(x     = gpmap$locinames[1], 
                                      y     = gpmap$mapnames[i], 
                                      group = gpmap$locinames[2], 
                                     colour = gpmap$locinames[2]))
    p <- p + geom_line(size = 2) + labs(x = gpmap$locinames[1], 
                                      y = "Genotypic value", 
                                      colour = gpmap$locinames[2]) 
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Original GPmap)",sep="\n"))
    p
  }
  return(figs)
}

plot3_orig <- function(gpmap) {
  # Plot gpmap object witho 3 loci
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp <- cbind(gpmap$genotype, gpmap$values[, i]) 
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i])
    p <- ggplot(data = gp, aes_string(x     = gpmap$locinames[1], 
                                      y     = gpmap$mapnames[i], 
                                      group = gpmap$locinames[2],  
                                     colour = gpmap$locinames[2])) 
    p <- p + geom_line(size = 2) + 
        facet_wrap(as.formula(paste('~',gpmap$locinames[3],sep="")), nrow = 1)
    p <- p + labs(y = 'Genotypic value', x = gpmap$locinames[1])
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Original GPmap)",sep="\n"))
    p
  }
  return(figs)
}

plot1_dec <- function(gpmap) {
  # Plot decomposed gpmap object with 1 locus
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp <- data.frame("values" = c(gpmap$values[, i], gpmap$values.mono[, i]),
                "map" = c(rep("GP map", dim(gpmap$values)[1]),
                          rep("Monotone map", dim(gpmap$values.mono)[1])))
    gp <- cbind(rbind(gpmap$genotype, gpmap$genotype), gp)
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i],"map")
    
    p <- ggplot(data = gp, aes_string(x = gpmap$locinames,  
                                      y = gpmap$mapnames[i],
                                      group = gpmap$locinames)) 
    p <- p + geom_point(size = 2)
    p <- p + facet_wrap(~map, nrow = 1, scales = "free")
    p <- p + labs(x = gpmap$locinames, y = 'Genotypic value')
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Decomposed GPmap)",sep="\n"))
    p
  }
  return(figs)
}

plot2_dec <- function(gpmap) {
  # Plot decomposed gpmap object with 2 loci
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp <- data.frame("values" = c(gpmap$values[, i], gpmap$values.mono[, i]),
                "map" = c(rep("GP map", dim(gpmap$values)[1]),
                          rep("Monotone map", dim(gpmap$values.mono)[1])))
    gp <- cbind(rbind(gpmap$genotype, gpmap$genotype), gp)
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i],"map")
    
    p <- ggplot(data = gp, aes_string(x = gpmap$locinames[1], 
                                      y = gpmap$mapnames[i],  
                                      group = gpmap$locinames[2], 
                                      colour = gpmap$locinames[2]))
    p <- p + geom_line(size = 2)
    p <- p + facet_wrap(~map, nrow = 2, scales = "free")
    p <- p + labs(x = gpmap$locinames[1], y = 'Genotypic value')
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Decomposed GPmap)",sep="\n"))
    p
  }
  return(figs)
}

plot3_dec <- function(gpmap) {
  # Plot decomposed gpmap object with 3 loci.
    i <- NULL
  figs <- foreach (i = 1:gpmap$nmaps) %do% {
    gp <- data.frame("values" = c(gpmap$values[, i], gpmap$values.mono[, i]),
                     "map" = c(rep("GP map", dim(gpmap$values)[1]),
                               rep("Monotone map", dim(gpmap$values.mono)[1])))
    gp <- cbind(rbind(gpmap$genotype, gpmap$genotype), gp)
    colnames(gp) <- c(gpmap$locinames, gpmap$mapnames[i], "map")
    p <- ggplot(data = gp, aes_string(x = gpmap$locinames[i],  
                                      y = gpmap$mapnames[i], 
                                      group = gpmap$locinames[2], 
                                      colour = gpmap$locinames[2]))
    p <- p + geom_line(size = 2)
    p <- p + facet_grid(as.formula(paste("map", "~", gpmap$locinames[3], 
                                         sep = "")), space = "free", scales = "free")
    p <- p + labs(x = gpmap$locinames[1], y = 'Genotypic value') 
    p <- p + labs(title = paste(gpmap$mapnames[i], "(Decomposed GPmap)",sep="\n"))
    p
    p
  }
  return(figs)
}

## print method for objects of class 'gpmap'

print.gpmap <- function(x, ...) {

  # map name(s)
  if (x$nmaps==1) {
    cat(sprintf("%d GP map: \n", x$nmaps))
    cat("\t", x$mapnames, "\n\n")
  } else {
    cat(sprintf("%d GP maps: \n", x$nmaps))
    cat("\t", x$mapnames, "\n\n")
  }

  # names of loci
  if (x$nloci==1) {
    cat('1 locus:')
  } else {
    cat(sprintf("%d loci: \n", x$nloci))
  }  
  cat("\t", x$locinames, "\n\n")

  #genotypic values
  cat("Summary of genotypic values:\n")
  print(summary(x$values))
  cat("\n\n")
  
  #degree of monotonicity if available
  if("degree.monotonicity" %in% names(x)) {
    cat("Degree of monotonicity (m):\n")
    print(x$degree.monotonicity)
    cat("\n")
  }

  #R^2_mono if available
  if("monoR2" %in% names(x)) {
    cat("R-square from isotone regression:\n")
    print(x$monoR2)
    cat("\n")
  }
}

#treat ´i´ as global during checking (to avoid NOTE)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("i"))