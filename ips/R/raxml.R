## PACKAGE: ips
## CALLED BY: USER
## AUTHOR: Christoph Heibl (at gmx.net)
## LAST UPDATE: 2014-07-30

raxml <- function(DNAbin, m = "GTRCAT", f, N, p, b, x, k,
                  partitions, outgroup, backbone = NULL, 
                  file = "fromR", exec, threads){
				
	# number of threads (PTHREADS only)
	# ---------------------------------
	if ( !missing(threads) )
    exec <- paste(exec, "-T", threads)
  
	# clear previous runs
	# -------------------
	unlink(list.files(pattern = "RAxML_")) 
  
  # substitution model
  # ------------------
	m <- match.arg(m, c("GTRCAT", "GTRCATX", 
	                    "GTRCATI", "GTRCATIX",
	                    "ASC_GTRCAT", "ASC_GTRCATX", 
	                    "GTRGAMMA", "GTRGAMMAX", 
	                    "GTRGAMMAI", "GTRGAMMAIX",
	                    "ASC_GTRGAMMA", "ASC_GTRGAMMAX"))
  m <- paste("-m", m)
  
  ## number of searches/replicates
  ## -----------------------------
  if ( is.character(N) ){
    N <- match.arg(N, c("autoFC", "autoMR", "autoMRE", "autoMRE_IGN"))
  }
  N <- paste("-N", N)
  
  ## random seeds
  ## ------------
	rs <- function(rseed, type = "p"){
	  if ( missing(rseed) ) rseed <- sample(1:999999, 1)
	  paste("-", type, " ",  rseed, sep = "")
	}
	p <- rs(p); x <- rs(x, type = "x")
  
  ## input file names
  ## ----------------
	rin <- c(s = paste("-s ", file, ".phy", sep = ""),
	         n = paste("-n", file),
	         napt = paste("-n ", file, ".APT", sep = ""))
 
	write.phy(DNAbin, paste(file, "phy", sep = ".")) ## write sequence file

  
  ## rout: raxml output file names
  ## -----------------------------
  output.types <- c("info", "bestTree", "bootstrap", "bipartitions")
  rout <- paste("RAxML_", output.types, ".", file, sep = "")
  names(rout) <- output.types
  
  ## algorithms
  ## ----------
  if ( missing(f) ) f <- "d"
  f <- match.arg(f, c("d", "a"))
  alg <- paste("-f", f, p) # add parsimony seed to algorithm
  if ( f == "a" ) alg <- paste(alg, x)
  if ( !missing(b) ) alg <- paste(alg, "-b", b)
  if ( missing(N) ) stop("the number of runs must be given (N)")
  
  ## outgroup
  ## --------
  if ( missing(outgroup) ){
    o <- ""
  } else {
    o <- outgroup %in% rownames(DNAbin)
    if ( !all(o) ){
      o <- paste(paste("\n  -", outgroup[!o]), collapse = "")
      stop(paste("outgroup names not in 'DNAbin':",   o))
    }
    o <- paste(outgroup, collapse = ",")
    o <- paste("-o", o)
  }

	# write partition file
	## ------------------
	if ( !missing(partitions) ) {
    if ( is.character(partitions) ){
      q <- partitions
    } else {
      q <- paste(partitions$type, ", ", 
                 partitions$locus, " = ", 
                 partitions$begin, "-", 
                 partitions$end, sep = "")
    }
		write(q, "partitionsFromR")
    multipleModelFileName <- " -q partitionsFromR "
	} else multipleModelFileName <- ""


	if ( !is.null(backbone) ){
	  write.tree(backbone, "backbone.tre")
	  g <- " -g backbone.tre"
	} else {
	  g <- " "
	}
  
  ## save branch lengths of bootstrap replicates
  ## -------------------------------------------
  if ( missing(k) ) k <- FALSE 
  k <- ifelse(k, "-k", "")
  
  
	## prepare and execute call
  ## ------------------------
	CALL <- paste(exec, alg, m, o, k,
	              multipleModelFileName, N, g, 
	              rin["s"], rin["n"])
  
	if ( length(grep("MPI", exec) > 0) ) system(paste("mpirun", CALL))
  print(CALL)
	system(CALL)
  
	res <- scan(rout["info"], quiet = TRUE, what = "char", sep = "\n")
	if ( length(grep("exiting", res)) > 0 )
	  stop("\n", paste(res, collapse = "\n"))
	
	## read results
	## ------------
	bestTree <- bipartitions <- bootstrap <- NULL
  info <- scan(rout["info"], what = "c", sep = "\n", quiet = TRUE)
	if ( f %in% c("a", "d") & missing(b) )
	  bestTree <- read.tree(rout["bestTree"])
	if ( f %in% c("a") )
	  bipartitions <- read.tree(rout["bipartitions"])
	if ( !missing(b) | !missing(x) )
	  bootstrap <- read.tree(rout["bootstrap"])
	
	obj <- list(info = info,
	            bestTree = bestTree,
	            bipartitions = bipartitions,
	            bootstrap = bootstrap)
	obj[sapply(obj, is.null)] <- NULL
	obj
}

# ##########################################################
# # 	test initial rearrangement
# ##########################################################
# 
# if ( optimize ){
#   n <- 0:4
#   call.opt <- c(paste("./raxmlHPC ", rs(), " -y -s ", file , " -m GTRCAT -n ST", n, sep = ""),
#                 paste("./raxmlHPC -f d -i 10 -m GTRCAT -s ", file ,
#                       " -t RAxML_parsimonyTree.ST", n, " -n FI", n,  sep = ""),
#                 paste("./raxmlHPC -f d -m GTRCAT -s ", file , 
#                       " -t RAxML_parsimonyTree.ST", n, " -n AI", n, sep = ""))
#   lapply(call.opt, system)
#   fixed <- paste("FI", n, sep = "")
#   auto <- paste("AI", n, sep = "")
#   res <- paste("RAxML_info", c(fixed, auto), sep = ".")
#   res <- lapply(res, scan, what = "c", sep = "", quiet = TRUE)
#   res <- as.numeric(sapply(res, function(x) x[grep("Program", x) - 1]))
#   names(res) <- c(fixed, auto)
#   FIXED <- mean(res[fixed])
#   AUTO <- mean(res[auto])
#   # set i
#   if ( AUTO > FIXED ) {
#     i <- x[grep("setting", x) + 1]
#     i <- as.numeric(gsub(",", "", i))
#   }
#   
#   
#   ##########################################################
#   # 	number of categories
#   ##########################################################
#   x <- paste("./raxmlHPC -f d -i ", i, " -m GTRCAT -s ", file ," -t RAxML_parsimonyTree.ST", sep = "")
#   system(paste(x, "0 -c 10 -n C10_0", sep = ""))
#   system(paste(x, "1 -c 10 -n C10_1", sep = ""))
#   system(paste(x, "2 -c 10 -n C10_2", sep = ""))
#   system(paste(x, "3 -c 10 -n C10_3", sep = ""))
#   system(paste(x, "4 -c 10 -n C10_4", sep = ""))
#   system(paste(x, "0 -c 40 -n C40_0", sep = ""))
#   system(paste(x, "1 -c 40 -n C40_1", sep = ""))
#   system(paste(x, "2 -c 40 -n C40_2", sep = ""))
#   system(paste(x, "3 -c 40 -n C40_3", sep = ""))
#   system(paste(x, "4 -c 40 -n C40_4", sep = ""))
#   system(paste(x, "0 -c 55 -n C55_0", sep = ""))
#   system(paste(x, "1 -c 55 -n C55_1", sep = ""))
#   system(paste(x, "2 -c 55 -n C55_2", sep = ""))
#   system(paste(x, "3 -c 55 -n C55_3", sep = ""))
#   system(paste(x, "4 -c 55 -n C55_4", sep = ""))
#   x <- scan("RAxML_info.C10_0", what="c", sep = "", 
#             quiet = TRUE)
#   C10_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_1", what="c", sep = "", 
#             quiet = TRUE)
#   C10_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_2", what="c", sep = "", 
#             quiet = TRUE)
#   C10_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_3", what="c", sep = "", 
#             quiet = TRUE)
#   C10_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_4", what="c", sep = "", 
#             quiet = TRUE)
#   C10_4 <- x[grep("Final", x)-1]
#   C10 <- mean(as.numeric(c(C10_0, C10_1, C10_2, C10_3, 			C10_4)))
#   
#   x <- scan("RAxML_info.C40_0", what="c", sep = "", quiet = TRUE)
#   C40_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_1", what="c", sep = "", quiet = TRUE)
#   C40_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_2", what="c", sep = "", quiet = TRUE)
#   C40_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_3", what="c", sep = "", quiet = TRUE)
#   C40_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_4", what="c", sep = "", quiet = TRUE)
#   C40_4 <- x[grep("Final", x)-1]
#   C40 <- mean(as.numeric(c(C40_0, C40_1, C40_2, C40_3, C40_4)))
#   
#   x <- scan("RAxML_info.C55_0", what = "c", sep = "", 		quiet = TRUE)
#   C55_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_1", what = "c", sep = "", 		quiet = TRUE)
#   C55_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_2", what = "c", sep = "", 		quiet = TRUE)
#   C55_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_3", what = "c", sep = "", 		quiet = TRUE)
#   C55_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_4", what = "c", sep = "", 		quiet = TRUE)
#   C55_4 <- x[grep("Final", x)-1]
#   C55 <- mean(as.numeric(c(C55_0, C55_1, C55_2, C55_3, C55_4)))
#   
#   C <- if (FIXED > AUTO) c(C10, FIXED, C40, C55) 		else c(C10, AUTO, C40, C55)
#   catnum <- c(10, 25, 40, 55)
#   DF <- cbind(catnum, C)
#   DF <- DF[order(DF[,2], decreasing = TRUE),]
#   numcat <- DF[1,1]
#   
# } # end of OPTIMIZE
# else {
#   i <- 10
#   numcat = 25
# }
