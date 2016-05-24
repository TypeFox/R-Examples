#' Solve a Linear Program (Maximization)
#'
#' \code{latteMax} uses LattE's maximize function to find the maximum of a linear objective function over the integers satisfying linearity constraints.  This makes use of the digging algorithm; see the LattE manual at \url{http://www.math.ucdavis.edu/~latte} for details.
#' 
#' @param objective a linear polynomial to pass to \code{\link{mp}}, see examples
#' @param constraints a collection of linear polynomial (in)equalities that define the feasibility region, the integers in the polytope
#' @param method method LP or cones
#' @param dir directory to place the files in, without an ending /
#' @param opts options; see the LattE manual at \url{http://www.math.ucdavis.edu/~latte}
#' @param quiet show latte output
#' @return the count.  if the count is a number has less than 10 digits, an integer is returned.  if the number has 10 or more digits, an integer in a character string is returned. you may want to use the gmp package's as.bigz to parse it.
#' @export latteMax
#' @examples
#' \dontrun{
#' 
#' latteMax("-2 x + 3 y", c("x + y <= 10", "x >= 0", "y >= 0"))
#'
#' 
#' df <- expand.grid(x = 0:10, y = 0:10)
#' df <- subset(df, x + y <= 10)
#' df$val <- apply(df, 1, function(v) -2*v[1] + 3*v[2])
#' df[which.max(df$val),]
#' 
#' library(ggplot2)
#' qplot(x, y, data = df, size = val)
#' 
#' 
#'
#' }
#' 
latteMax <- function(objective, constraints, 
  method = c("lp","cones"), dir = tempdir(), 
  opts = "", quiet = TRUE
){
	
  warning("this function is experimental.", call. = FALSE)	
	
  method <- match.arg(method)	

  ## parse objective
  if(is.character(objective)) objective <- mp(objective)
  stopifnot(is.linear(objective))
  

  ## parse constraints into the poly <= 0 format
  if(is.character(constraints) && length(constraints) > 1){
  	
  	parsedCons <- as.list(rep(NA, length(constraints)))
  	
  	geqNdcs <- which(str_detect(constraints, " >= "))
  	leqNdcs <- which(str_detect(constraints, " <= "))  	
  	eeqNdcs <- which(str_detect(constraints, " == "))  	  	
  	eqNdcs <- which(str_detect(constraints, " = "))  	  	  	

    if(length(geqNdcs) > 0){
      tmp <- strsplit(constraints[geqNdcs], " >= ")
      parsedCons[geqNdcs] <- 
        lapply(tmp, function(v) mp(v[2]) - mp(v[1]))
    } 
    
    if(length(leqNdcs) > 0){
      tmp <- strsplit(constraints[leqNdcs], " <= ")
      parsedCons[leqNdcs] <- 
        lapply(tmp, function(v) mp(v[1]) - mp(v[2]))
    }     
    
    if(length(eeqNdcs) > 0){
      tmp <- strsplit(constraints[eeqNdcs], " == ")
      parsedCons[eeqNdcs] <- 
        lapply(tmp, function(v) mp(v[1]) - mp(v[2]))
    }          	
    
    if(length(eqNdcs) > 0){
      tmp <- strsplit(constraints[eqNdcs], " = ")
      parsedCons[eqNdcs] <- 
        lapply(tmp, function(v) mp(v[1]) - mp(v[2]))
    } 
    
    linearityNdcs <- sort(c(eeqNdcs, eqNdcs))
      	
    constraints <- parsedCons
    class(constraints) <- "mpolyList"
    
    if(!all(is.linear(constraints))){
      stop("all polynomials must be linear.", call. = FALSE)
    }    
  }
  
  
  ## mpolyListToMat is in file count.r
  matFull <- mpolyListToMat(c(list(objective),constraints))
  
  
  
  ## convert constraints to latte hrep code
  mat <- cbind(
    -matFull[-1,"coef",drop=FALSE], 
    -matFull[-1,-ncol(matFull)]
  )
  consCode <- paste(nrow(mat), ncol(mat))
  consCode <- paste0(consCode, "\n")
  consCode <- paste0(consCode, 
    paste(apply(unname(mat), 1, paste, collapse = " "), 
    collapse = "\n")
  )   
  
  if(length(linearityNdcs) > 0){
    consCode <- paste0(consCode, "\n")
    consCode <- paste0(consCode,
      paste0("linearity ", length(linearityNdcs), " ", 
        paste(linearityNdcs, collapse = " "))
    )
  }  
  
  
  
  ## convert objective to latte hrep code
  mat <- cbind(
    -matFull[1,"coef",drop=FALSE], 
    -matFull[1,-ncol(matFull),drop=FALSE]
  )
  objCode <- paste(nrow(mat), ncol(mat))
  objCode <- paste0(objCode, "\n")
  objCode <- paste0(objCode, 
    paste(apply(unname(mat), 1, paste, collapse = " "), 
    collapse = "\n")
  )     
  
      
	
  ## make dir to put latte files in (within the tempdir) timestamped
  timeStamp <- as.character(Sys.time())
  timeStamp <- chartr("-", "_", timeStamp)
  timeStamp <- chartr(" ", "_", timeStamp)
  timeStamp <- chartr(":", "_", timeStamp)
  dir2 <- file.path2(dir, timeStamp)
  suppressWarnings(dir.create(dir2))
	


  ## write code file
  writeLines(consCode, con = file.path2(dir2, "maxCode"))
  writeLines(objCode, con = file.path2(dir2, "maxCode.cost"))  



  ## switch to temporary directory
  oldWd <- getwd()
  setwd(dir2)
  
  
  ## run count
  if(.Platform$OS.type == "unix"){    	
  	# bizarrely, latte-maximize returns its output as stderr
    system(
      paste(
        file.path2(getOption("lattePath"), "latte-maximize"), 
        opts, 
        file.path2(dir2, "maxCode 2> out.txt")
      ),
      intern = FALSE, ignore.stderr = FALSE
    )     
    outPrint <- readLines(file.path2(dir2, "out.txt"))
  } else { # windows    	
    matFile <- file.path2(dir2, "maxCode 2> out.txt")
    matFile <- chartr("\\", "/", matFile)
    matFile <- str_c("/cygdrive/c", str_sub(matFile, 3))
    system(
      paste(
        paste0("cmd.exe /c env.exe"),
        file.path(getOption("lattePath"), "latte-maximize"), 
        opts, 
        matFile
      ),
      intern = FALSE, ignore.stderr = FALSE
    )
    outPrint <- readLines(file.path2(dir2, "out.txt"))
  }
  

  ## print count output when quiet = FALSE
  if(!quiet) cat(outPrint, sep = "\n")

  
  ## migrate back to original working directory
  setwd(oldWd)
  

  ## parse output
  lookFor <- ifelse(method == "cones", "An optimal", 
    "A vertex which we found via LP")
  par <- outPrint[which(str_detect(outPrint, lookFor))]
  par <- strsplit(par, ": ")[[1]][2]
  if(method == "cones") par <- str_sub(par, 2, -3)
  if(method == "lp") par <- str_sub(par, 2, -2)
  par <- as.integer(strsplit(par, " ")[[1]])
  names(par) <- colnames(matFull)[1:(ncol(matFull)-1)]
  
  lookFor <- ifelse(method == "cones", "The optimal value is", 
    "The LP optimal value is")  
  val <- outPrint[which(str_detect(outPrint, lookFor))]
  val <- strsplit(val, ": ")[[1]][2]  
  if(method == "cones") val <- str_sub(val, 1, -2)
  val <- as.integer(val)
  
  ## out
  list(par = par, value = val)
}




