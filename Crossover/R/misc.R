# Alternative to getDesign without quotes. (Not exported yet.)
design <- function(d, character.only = FALSE) {
  if (!character.only) {
    d <- as.character(substitute(d))
  }
  return(get(d, envir=Crossover.env))
}

dput2 <- function(x) {
  paste(capture.output(dput(x)), collapse = " ")
}

# Substitute Treatments with other treatments (needed for combine.PBIB):
subs.treatment <- function(design, new.treatments) {
  x <- matrix(0, dim(design)[1], dim(design)[2])
  for (t in 1:max(design)) {
   x[design==t] <- new.treatments[t] 
  }    
  return(x)
}

# pbib <- clatworthy.r94
# coDesign <- get("williams4t", Crossover:::Crossover.env)

combine.PBIB <- function() {
  
  balancedDesigns1 <- c("andersonPreece", "anderson", "archdeacon", "atkinson3t", "atkinson4t", 
                        "atkinson5t", "balaam3t", "balaam4t", "balaam5t", "balaam6t", 
                        "berenblut3t", "berenblut4t", "berenblut5t", "davisHall7tb", 
                        "federerAtkinson3ta", "federerAtkinson3tb", "federerAtkinson4ta", 
                        "federerAtkinson4tb", "federerAtkinson5ta", "fletcher9", "iqbalJones1", 
                        "iqbalJones2", "iqbalJones3", "iqbalJones4", "iqbalJones5", "iqbalJones7", 
                        "iqbalJones11", "iqbalJones17", "iqbalJones23", "lewisFletcherMatthews3", 
                        "orthogonalLatinSquare3t", "orthogonalLatinSquare4t", "orthogonalLatinSquare5t", 
                        "orthogonalLatinSquare7t", "pattersonLucasExtraPeriod30", "pattersonLucasExtraPeriod31", 
                        "pattersonLucasExtraPeriod32", "pattersonLucasExtraPeriod33", 
                        "pattersonLucasExtraPeriod34", "pattersonLucasExtraPeriod35", 
                        "pattersonLucasExtraPeriod36", "pattersonLucasExtraPeriod37", 
                        "pattersonLucasExtraPeriod38", "pattersonLucasExtraPeriod39", 
                        "pattersonLucasExtraPeriod40", "pattersonLucasExtraPeriod41", 
                        "pattersonLucasExtraPeriod42", "pattersonLucasExtraPeriod43", 
                        "pattersonLucasExtraPeriod44", "pattersonLucasExtraPeriod45", 
                        "pattersonLucasExtraPeriod46", "pattersonLucasExtraPeriod47", 
                        "pattersonLucasExtraPeriod48", "pattersonLucasExtraPeriod49", 
                        "pattersonLucasExtraPeriod86", "pattersonLucasPltT1", "pattersonLucasPltT3", 
                        "pattersonLucasPltT4", "pattersonLucasPltT5", "pattersonLucasPltT7", 
                        "pattersonLucasPltT8", "pattersonLucasPltT9", "pattersonLucasPltT10", 
                        "pattersonLucasPltT12", "pattersonLucasPltT13", "pattersonLucasPltT15", 
                        "pattersonLucasPltT16", "pattersonLucasPltT17", "pattersonLucasPltT18", 
                        "pattersonLucasPltT19", "pattersonLucasPltT20", "pattersonLucasPltT21", 
                        "pattersonLucasPltT22", "pattersonLucasPltT23", "prescott1", 
                        "prescott2", "quenouille3t2", "quenouille4t1", "quenouille4t2", 
                        "quenouille4t3", "switchback3t", "switchback4t", "switchback5t", 
                        "switchback6t", "switchback7t", "williams3t", "williams4t", "williams5t", 
                        "williams6t", "williams7t", "williams8t", "williams9t")
  balancedDesigns2 <- c("bateJones5t", "bateJones8t", "iqbalJones31", "quenouille3t1")
  balancedDesigns <- c(balancedDesigns1, balancedDesigns2)
  
  # List of PBIB(2):
  # path <- system.file("data", package="Crossover")
  x <- load("pkg/Crossover/data/clatworthy1.rda")
  designs <- c()
  # Available crossover designs:
  summary <- buildSummaryTable()
  #summary <- summary[grep("williams", summary$dataset), ]
  summary <- summary[summary$dataset %in% balancedDesigns, ]
  for (pbib.name in x) {
    pbib <- get(pbib.name)
    t.per.block <- dim(pbib)[1] # verified with book
    pbib.title <-  attr(pbib, "title")
    cat("  Creating designs for ", pbib.name, "\n")
    coDesigns <- summary[ summary$t==t.per.block, ]
    if (dim(coDesigns)[1]>0) {
      for (row in 1:(dim(coDesigns)[1])) {
        #coDesign <- get(coDesigns$dataset[row], envir=Crossover:::Crossover.env)
        coDesign <- get(coDesigns$dataset[row], envir=Crossover.env)
        coTitle <- coDesigns[row, "title"]
        coReference <- coDesigns[row, "reference"]
        title <- paste(paste("PB2.", tail(strsplit(pbib.title, " ")[[1]], 1), "-",  coDesigns$dataset[row], sep=""), "- Design from combining PBIB(2)", pbib.title, "with", coTitle)
        cat("Creating design ", title, "\n")
        design <- c()
        for (col in 1:(dim(pbib)[2])) {
          design <- cbind(design, subs.treatment(coDesign, pbib[, col]))
        }
        #print(design)
        design.name <- tolower(paste("pb2.", coDesigns$dataset[row], ".", tail(strsplit(pbib.title, " ")[[1]], 1), sep=""))
        attr(design, "reference") <- paste(coReference, "Clatworthy, Willard H. Tables of two-associate-class partially balanced designs. US Government Printing Office, 1973.", sep="\n")
        attr(design, "signature") <- paste("p=", dim(design)[1],", n=", dim(design)[2],", t=", max(design), sep="")
        attr(design, "title") <- title
        #designI <- as.matrix(as.integer(design), nrow=dim(design)[1])
        #if (!isTRUE(all.equal(max(abs(design - designI)), 0))) stop("Design and DesignI differ!")
        assign(design.name, design)
        designs <- c(designs, design.name)
      }
    } else {
      cat(" == No crossover designs for ",t.per.block," treatments.\n")
    }    
  }
  save(list=designs, file="/home/kornel/pbib2combine.rda")
  
  # Generating alias
  # designs <- load(paste(system.file("data", package="Crossover"), "pbib2combine.rda", sep="/"), envir=Crossover:::Crossover.env)
  designs <- load(paste(system.file("data", package="Crossover"), "pbib2combine.rda", sep="/"), envir=Crossover.env)
  for (design in designs) {
    cat(paste("\\alias{",design,"}\n", sep=""))
  }  
}