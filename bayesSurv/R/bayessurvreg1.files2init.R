#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       bayessurvreg1.files2init.R          ####
####                                                 ####
#### FUNCTIONS:  bayessurvreg1.files2init            ####
#########################################################

### ======================================
### bayessurvreg1.files2init
### ======================================
bayessurvreg1.files2init <- function(dir = getwd(), 
                                     row, 
                                     kmax)
{
  files <- dir(dir)    ## character vector with available files (except D and mixture files)
  fnames <- paste(c("iteration", "beta", "b", "Y", "r", "otherp", "u"), ".sim", sep = "")  
  compnames <- c("iter", "beta", "b", "y", "r", "otherp", "u")

  if (missing(row)){
    if (sum(!is.na(match(files, "iteration.sim")))){
      iter <- read.table(file = paste(dir, "/", "iteration.sim", sep = ""), header = TRUE, skip = 0)
      skip <- nrow(iter) - 1
    }
    else
      skip <- 0
  }
  else
    skip <- row - 1
  
  init <- list()

    ## Read everything except mixture and D matrix
    ##  (some files (r, Y, u, b) can be without header)
    ## * if there is no header, read the first line as initial values
  for (i in 1:length(fnames)){
    init[[compnames[i]]] <- NULL
    skip.now <- skip
    if (sum(!is.na(match(files, fnames[i])))){
      help <- read.table(file = paste(dir, "/", fnames[i], sep = ""), header = FALSE, nrows = 1, as.is = TRUE)
      if (is.character(help[1, 1]))
        header <- TRUE
      else{
        header <- FALSE
        skip.now <- 0
      }        
      help <- scan(file = paste(dir, "/", fnames[i], sep = ""), skip = skip.now + 1*header, nlines = 1)
      if (!length(help)) stop("Incorrect 'row' or 'header' parameter.")
      init[[compnames[i]]] <- as.numeric(help)
    }
  }

    ## Read D matrix
  init$D <- NULL
  if (sum(!is.na(match(files, "D.sim")))){
    header <- TRUE
    help <- scan(file = paste(dir, "/", "D.sim", sep = ""), skip = skip + 1*header, nlines = 1)
    if (!length(help)) stop("Incorrect 'row' or 'header' parameter.")
    init$D <- as.numeric(help[-1])     ## Remove determinant
  }
  
    ## Read mixture
  init$mixture <- NULL
  if (sum(!is.na(match(files, "mixmoment.sim"))) + sum(!is.na(match(files, "mweight.sim"))) +
      sum(!is.na(match(files, "mmean.sim"))) + sum(!is.na(match(files, "mvariance.sim"))) == 4){

    header <- TRUE
    
    if (header){
      mix <- read.table(file = paste(dir, "/mweight.sim", sep = ""), nrows = 1)
      kmax <- length(mix)
      mix <- read.table(paste(dir, "/mmean.sim", sep = ""), nrows = 1)
      kmax2 <- length(mix)
      if (kmax != kmax2) stop("Different kmax indicated by files mweight.sim and mmean.sim.")
      mix <- read.table(paste(dir, "/mvariance.sim", sep = ""), nrows = 1)     
      kmax2 <- length(mix)
      if (kmax != kmax2) stop("Different kmax indicated by files mweight.sim and mvariance.sim.")         
    }
    if (!header & missing(kmax)){
      stop("kmax must be given")
    }      
    
    kk <- scan(file = paste(dir, "/", "mixmoment.sim", sep = ""), skip = skip + 1*header, nlines = 1)[1]
    mweight <- scan(file = paste(dir, "/", "mweight.sim", sep = ""), skip = skip + 1*header, nlines = 1)
    mmean <- scan(file = paste(dir, "/", "mmean.sim", sep = ""), skip = skip + 1*header, nlines = 1)    
    mvariance <- scan(file = paste(dir, "/", "mvariance.sim", sep = ""), skip = skip + 1*header, nlines = 1)
    k.now <- length(mweight)
    if (k.now == 0) stop("Invalid mixture weights read.")
    k.now2 <- length(mmean)
    if (k.now != k.now2) stop("Different k indicated by files mweight.sim and mmean.sim.")
    k.now2 <- length(mvariance)
    if (k.now != k.now2) stop("Different k indicated by files mweight.sim and mmean.sim.")    

    init$mixture <- c(kk, mweight, rep(0, kmax-k.now), mmean, rep(0, kmax-k.now), mvariance, rep(0, kmax-k.now))
  }      

  return(init)
}

