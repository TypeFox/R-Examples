#require("mva")

number.of.tips.sim <- function(x){
  ## extracts number of tips from a simulation file
  scan(file = x, skip = 3, nlines = 1, what = list("", 0), strip.white = TRUE,
       sep = "=", quiet = TRUE)[[2]][1]
}

number.of.tips.inp <- function(x){
  ## extracts number of tips from an inp file
  scan(file = x, skip = 3, nlines = 1,
       what = list("", "", 0, "", "", "", "", "", ""), strip.white = TRUE,
       quiet = TRUE)[[3]][1]
}

number.of.tips.pdi <- function(x){
  ## extracts number of tips from a pdi file
  scan(file = x, skip = 0, nlines = 1, quiet = TRUE)
}

tips.names.inp <- function(x){
  ## extracts names of tips from an inp file
  scan(x, list("", 0, "", 0), skip = 8, nmax = number.of.tips.inp(x),
       quiet = TRUE)[[1]]
}

tips.names.pdi <- function(x){
  ## extracts names of tips from an pdi file
  scan(x, list("", 0, 0), skip = 1 + 2 * number.of.tips.pdi(x),
       quiet = TRUE)[[1]]
}

number.of.simulations <- function(x){
    ## number of simulations
    scan(file = x, skip = 17 + number.of.tips.sim(x), nlines = 1,
         what = list("", "", "", 0, "", ""), quiet = TRUE)[[4]][1]
}

num.sim.tips <- function(x){
          ##returns the number of tips and the number of simulations
    tips <- scan(file = x, skip = 3, nlines = 1, what = list("", 0),
                 strip.white = TRUE, sep = "=", quiet = TRUE)[[2]][1]
    c(tips, scan(file = x, skip = 17 + tips, nlines = 1,
                 what = list("", "", "", 0, "", ""), quiet = TRUE)[[4]][1])
}


scan.simulation.file <- function(x, max.num = 0){
  ## extract, from each simulation file,
    ##conly the two columns with simulated data
  num.tips <- number.of.tips.sim(x)
  matrix(unlist(scan(x, list("", 0, 0, "", 0, 0, 0), skip = 21 + num.tips,
                     nlines = max.num * num.tips, quiet = TRUE)[2:3]),
         ncol = 2)
}



scan.inp.file <- function(x){
  ## extract, form an inp file, the two columns with the data
  matrix(unlist(scan(x, list("", 0, "", 0), skip = 8,
                     nmax = number.of.tips.inp(x),
                     quiet = TRUE)[c(2, 4)]), ncol = 2)
}


scan.pdi.file <- function(x){
  ## extract, form a pdi file, the two columns with the data
  matrix(unlist(scan(x, list("", 0, 0), skip = 1 + 2 * number.of.tips.pdi(x),
                     quiet = TRUE)[c(2, 3)]), ncol = 2)
}




read.pdi.data <- function(input.pdi.files, variable.names = NULL){
  ## this function takes one or more pdi files and extracts
  ##  the tips names and character values
  ## you can provide a column of names for the variables
  
  ## error checks
  ## (they only apply when more than one inp file, but to make things
  ##  simpler I use the same code for all cases).
    numbers.tips.check <- as.vector(sapply(input.pdi.files,
                                           number.of.tips.pdi))
    if(length(unique(numbers.tips.check))>1)
        stop("Pdi files have different numbers of tips")
 
    tips.order.check <- sapply(input.pdi.files, tips.names.pdi)
    if(!all(apply(tips.order.check, 2,
                  function(x){tips.order.check[, 1] == x})))
        stop("Species names (and/or their order) are different")    


  ## Put the files together
  tmp <- as.data.frame(matrix(unlist(lapply(input.pdi.files, scan.pdi.file)),
                              ncol = 2 * length(input.pdi.files)))
  if(!is.null(variable.names))
    if(length(variable.names) != length(names(tmp)))
        warning("The vector of variable names is not of the \nsame length as the number of variables. \nVector of variable names not added") else names(tmp) <- variable.names
  tmp <- cbind(tips.order.check[, 1], tmp)
  names(tmp)[1] <- "Tips"
  class(tmp) <- c("pdi.file", "data.frame")
  tmp
}




read.inp.data <- function(input.inp.files, variable.names = NULL){
  ## this function takes one or more inp files and extracts the tips
  ##  names and character values
  ## you can provide a column of names for the variables
  
  ## error checks
  ## (they only apply when more than one inp file, but to make
  ##  things simpler I use the same
  ## code for all cases).
    ##    numbers.tips.check <- t(sapply(input.inp.files, number.of.tips.inp)) #up to version 1.0.1

    numbers.tips.check <- as.vector(sapply(input.inp.files, number.of.tips.inp))
    if(length(unique(numbers.tips.check))>1)
        stop("Inp files have different numbers of tips")
 
    tips.order.check <- sapply(input.inp.files, tips.names.inp)
    if(!all(apply(tips.order.check, 2,
                  function(x){tips.order.check[, 1] == x})))
        stop("Species names (and/or their order) are different")    


  # Put the files together
  tmp <- as.data.frame(matrix(unlist(lapply(input.inp.files, scan.inp.file)),
                              ncol = 2 * length(input.inp.files)))
  if(!is.null(variable.names))
    if(length(variable.names) != length(names(tmp)))
        warning("The vector of variable names is not of the \nsame length as the number of variables. \nVector of variable names not added") else names(tmp) <- variable.names
  tmp <- cbind(tips.order.check[, 1], tmp)
  names(tmp)[1] <- "Tips"
  class(tmp) <- c("inp.file", "data.frame")
  tmp
}


read.phylip.data <- function(input.phylip.file, variable.names = NULL) {
    ## takes a data input file such as used by PHYLIP
    num.species.and.variables <- scan(file = input.phylip.file, n = 2,
                                      quiet = TRUE)
    dataB <- scan(file = input.phylip.file, what = "", skip = 1, quiet = TRUE)
    dataB <- matrix(dataB, ncol = num.species.and.variables[2] + 1, byrow = TRUE)
    species <- dataB[, 1]
    dataB <- dataB[, -1]
    dataB <- as.matrix(dataB)
    dataB <- as.data.frame(matrix(as.numeric(dataB),
                                  nrow = num.species.and.variables[1], byrow = FALSE))
    if(!is.null(variable.names))
    if(length(variable.names) != length(names(dataB)))
        warning("The vector of variable names is not of the \nsame length as the number of variables. \nVector of variable names not added") else names(dataB) <- variable.names
    dataB <- cbind(species, dataB)
    names(dataB)[1] <- "Tips"
    class(dataB) <- c("phylip.file", "data.frame")
    dataB
}



read.sim.data <- function(sim.files, inp.files = NULL, pdi.files = NULL,
                          phylip.file = NULL, variable.names = NULL,
                          other.variables = NULL, max.num = 0){
  # can enter sim and inp file without using explicit arguments; otherwise, for pdi, 
  # have to make it explicit that is pdi.file.

    
  # inp: one inp file; can either be after reading it (use function read.inp.data)
  # or before having read it.
  # If the file has already been read, then it is just a regular
  # R data frame, and the name is NOT quoted.
  # IF the file has still not been read, its name is quoted, and you
  # need to specify the appropriate path.
  # IMPORTANT: the inp file, as output from read.inp.data, is just a data frame with
  # one column of tips names and columns (multiples of 2) of variables.
  # If the inp file is given it is added to the output file;

  # In most cases you will want to have it done this way
  # variable.names: names of variables
  # other.variables: a vector of other variables to be added, that were not originally on the
  # sim files
  # max.num: if given, the number of simulations to use 


  ## To improve: deal with both inp and pid, in any possible order.
  ## recognize that a certain file is pdi.

rep.matrix.data.frame <- function(x, times){
#  browser()
  if(is.matrix(x)) {
    dimnames(x) <- NULL
    structure(apply(x, 2, rep, times))
  }

  else {as.data.frame(lapply(x, rep, times))}
##                                        # i.e., it is data.frame
##    tmp <- data.frame(matrix(0, nrow = nrow(x) * times, ncol = ncol(x)))
##    for (i in 1:ncol(x)) {tmp[[i]] <- rep(as.matrix(x[[i]]), times)}                     
##    names(tmp) <- NULL
##    tmp}
}
  original.data.files <- NULL  
  number.sim.files <- length(sim.files)
  numbers.check <- t(sapply(sim.files, num.sim.tips))
  if(is.null(inp.files) & is.null(pdi.files) & is.null(phylip.file))
      start.seq <- 1 else start.seq <- 0
##used to determine how many copies of the sim.counter, tips, etc,
   ##vectors to make when producing the output
  #some error checks
  if(!is.null(inp.files) & !is.null(pdi.files) & !is.null(phylip.file))
      stop("Only reads either inp or pdi or phylip, not both")
  if(length(unique(numbers.check[, 1]))>1)
      stop("Simulation files have different numbers of tips")
  if(length(unique(numbers.check[, 2]))>1)
      stop("Simulation files have different numbers of simulations")
  number.of.tips <- numbers.check[1, 1]
  number.of.simulations <- numbers.check[1, 2]
  ### I could do without these, or make a single object,
     ##  but they are not a big deal, and greatly enhance readability of code.

  # get species ---tips--- names and simulation counter
  tips.in.simulation <- unlist(scan(sim.files[1], list("", 0, 0, "", 0, 0, 0),
                                    skip = 21 + number.of.tips, quiet = TRUE,
                                    nlines = number.of.tips)[1])

  # error check for inp & pdi files
  if( is.null(inp.files) & is.null(pdi.files) & is.null(phylip.file))
      name.of.columns <- paste("V", seq(1:(2 * number.sim.files)), sep = "")
  else { if(!is.null(inp.files)){
# first determine if inp.files have already been read or not,  and process accordingly
    if(!is.data.frame(inp.files)) inp.files <- read.inp.data(inp.files)
    else inp.files <- inp.files
    if(number.of.tips != length(inp.files[, 1]))
        stop("Number of tips of inp and simualtion file(s) differ")
    if(!all(inp.files[, 1] == tips.in.simulation))
        stop("Species names or their order in inp and sim are different")
#    sim.counter <- rep(0, length(inp.files[, 1]))
    inp.files <- cbind(sim.counter = rep(0, length(inp.files[, 1])), inp.files)
                         ## added the simulation counter, which takes value 0
    name.of.columns <- names(inp.files)[-c(1, 2)] # will use later

    # to allow to combine with the simulation files later
    # and to minimize memory use and using rbind and cbind on a data frame, 
    # I make inp.files a matrix
    inp.files <- as.matrix(inp.files[, -c(1, 2)])
    dimnames(inp.files) <- NULL
    original.data.files <- inp.files
    rm(inp.files)
  }
  if(!is.null(pdi.files)){
# first determine if pdi.files have already been read or not, and process accordingly
    if(!is.data.frame(pdi.files))
        pdi.files <- read.pdi.data(pdi.files)
    else pdi.files <- pdi.files
    if(number.of.tips != length(pdi.files[, 1]))
        stop("Number of tips of pdi and simualtion file(s) differ")
    if(!all(pdi.files[, 1] == tips.in.simulation))
        stop("Species names or their order in pdi and sim are different")
#    sim.counter <- rep(0, length(inp.files[, 1]))
    pdi.files <- cbind(sim.counter = rep(0, length(pdi.files[, 1])), pdi.files)
         ## added the simulation counter, which takes value 0
    name.of.columns <- names(pdi.files)[-c(1, 2)] # will use later

    # to allow to combine with the simulation files later
    # and to minimize memory use and using rbind and cbind on a data frame, 
    # I make inp.files a matrix
    pdi.files <- as.matrix(pdi.files[, -c(1, 2)])
    dimnames(pdi.files) <- NULL
    original.data.files <- pdi.files
    rm(pdi.files)
  } 
  if(!is.null(phylip.file)){
# first determine if phylip.file has already been read or not, and process accordingly
    if(!is.data.frame(phylip.file))
        phylip.file <- read.phylip.data(phylip.file)
    else phylip.file <- phylip.file
    if(number.of.tips != length(phylip.file[, 1]))
        stop("Number of tips of phylip and simualtion file(s) differ")
    if(!all(phylip.file[, 1] == tips.in.simulation))
        stop("Species names or their order in phylip and sim are different")
#    sim.counter <- rep(0, length(inp.files[, 1]))
    phylip.file <- cbind(sim.counter = rep(0, length(phylip.file[, 1])),
                         phylip.file)
        ## added the simulation counter, which takes value 0
    name.of.columns <- names(phylip.file)[-c(1, 2)] # will use later

    # to allow to combine with the simulation files later
    # and to minimize memory use and using rbind and cbind on a data frame, 
    # I make inp.files a matrix
    phylip.file <- as.matrix(phylip.file[, -c(1, 2)])
    dimnames(phylip.file) <- NULL
    original.data.files <- phylip.file
    rm(phylip.file)
  } }

##   names of columns:
##   if there is not, then name of columns is whatever was found above
##   (null if no inp file, inp cols otherwise); thus, no need to deal
##   with this case
##     if there is such a variable, check it is of the right size;
##      if it is, assign it to name of columns
  
##     if there is a variable.names vector, make this the name of
##       colums for the inp file (fi it exists)
##     and the simulation file; otherwise, use the col. names
##    from the the inp file, if it exists.

if(!is.null(variable.names))
    if(length(variable.names) != (number.sim.files * 2))
        warning("The vector of variable names is not the same length as the number of variables. \nVector of variable names not added")
    else name.of.columns <- variable.names
  # now prepare the final name.of.columns
  name.of.columns <- c("sim.counter", "Tips", name.of.columns)
     
  if(max.num) number.of.simulations <- max.num
  
##   Create the return data frame by combining the above
##   what we do: rbind the inp to the simulation output; then cbind
##       the counter and tips; then we give the
##             columns their names.
##       This is a little bit hard to read, but avoids creating
##            intermediate objects.
##       The simulation files are read using lapply of the scan.simulation.file
##           on all sim files, and turning that into a matrix.
##     If there is the other.variables column, it is added here
##       too (after checking right length)
                               
    if (is.null(other.variables)) 
        structure(data.frame(rep(seq(from = start.seq,
                                     to = number.of.simulations), 
            rep(number.of.tips, number.of.simulations + 1 - start.seq)), 
            rep(tips.in.simulation, number.of.simulations + 1 - start.seq), 
            rbind(original.data.files, matrix(unlist(lapply(sim.files, 
                scan.simulation.file, max.num = number.of.simulations)),
                                              ncol = 2 * number.sim.files))),
                  names = name.of.columns, 
            class = c("simul.phylog", "data.frame"))
    else {  # there are other variables
        if (is.null(ncol(other.variables))) {
            name.of.vector <- deparse(substitute(other.variables))
            dim(other.variables) <- c(length(other.variables), 
                1)  #always a matrix
            dimnames(other.variables) <- list(NULL, name.of.vector)
            ##but try to get name to something reasonable
        }
        if (length(other.variables[, 1]) == number.of.tips) {
            if (is.data.frame(other.variables)) 
                name.of.columns <- c(name.of.columns, names(other.variables))
            else name.of.columns <- c(name.of.columns,
                                      dimnames(other.variables)[[2]])
            structure(data.frame(rep(seq(from = start.seq,
                                         to = number.of.simulations), 
                rep(number.of.tips, number.of.simulations + 1-start.seq)), 
                rep(tips.in.simulation, number.of.simulations + 1-start.seq),
                                 rbind(original.data.files,
                                       matrix(unlist(lapply(sim.files, 
                                                            scan.simulation.file,
                                                            max.num = number.of.simulations)), 
                                              ncol = 2 * number.sim.files)),
                                 rep.matrix.data.frame(other.variables, 
                                                       number.of.simulations + 1-start.seq)),
                      names = name.of.columns, 
                      class = c("simul.phylog", "data.frame"))
        }
        else {
            warning("The length of the 'other variables' is different from the number of tips; the 'other variables' have not been added")
            structure(data.frame(rep(seq(from = start.seq,
                                         to = number.of.simulations), 
                                     rep(number.of.tips,
                                         number.of.simulations + 1-start.seq)), 
                                 rep(tips.in.simulation,
                                     number.of.simulations + 1-start.seq),
                                 rbind(original.data.files,
                                       matrix
                                       (unlist(lapply(sim.files, 
                                                      scan.simulation.file,
                                                      max.num =
                                                      number.of.simulations +
                                                      1-start.seq)), 
                                        ncol = 2 * number.sim.files))),
                      names = name.of.columns, 
                      class = c("simul.phylog", "data.frame"))
        }
    }
}


# summary.simul.phylog <- function(x){}

# plot.simul.phylog <- function(x){}



summary.phylog.lm <- function(object, ...){

## note: we use quantile, so there is linear interpolation to find quantiles
##  if(class(object) != "phylog.lm")
## stop("Object mas be of class phylog.lm,
    ##such as obtained from the lm.phylog function")
  num.simul <- dim(object$MarginalTests)[1]-1
  F.value.sim <- as.matrix(object$MarginalTests[-1, ][-1])
  F.value.original <- object$MarginalTests[1, ][-1]
  ## make that into a matrix by repeating the value as many times as needed
  ## takes up space but seems more efficient than other solutions
  F.value.original <- matrix(rep(as.matrix(F.value.original), num.simul),
                             nrow = num.simul, byrow = TRUE)
##  sim.larger <- F.value.sim>F.value.original
##  count.larger <- apply(sim.larger, 2, function(object)
        ##{length(object[object == TRUE])})
##  count.larger <- apply(F.value.sim>F.value.original, 2, sum)
##  p.value <- (count.larger + 1)/(num.simul + 1)

## this is what p-value does: count the number of cases where the simulated F is
## larger than the one from inp, add one, and divide by num.simulations + 1.
  p.value <- (apply(F.value.sim>F.value.original, 2, sum) + 1)/(num.simul + 1) 
  
  quant.F.value <- apply(F.value.sim, 2, quantile,
                         probs = c(.5, .9, .95, .99, .999))
  structure(list(call = object$call, original.model = object$Fits[1, ][-1],
                 original.Fvalue = object$MarginalTests[1, ][-1], 
                 p.value = p.value, quant.Fvalue = quant.F.value,
                 num.simul = num.simul), class = "summary.phylog.lm")
}
           
print.summary.phylog.lm <- function(x, ...){
  cat("\nCall: \n")
  print(x$call)
  cat("\nModel from inp data:\n")
  print(unlist(x$original.model))
  cat("\nF-test's from inp data:\n")
  print(unlist(x$original.Fvalue))
  cat("\nP-value:\n")
  print(x$p.value)
  cat("\nQuantiles of F-value distribution:\n")
  print(x$quant.Fvalue)
  cat("\nNumber of simulations used in analyses: ",
      format(x$num.simul, ...), "\n")
  invisible(x)
}
    

plot.phylog.lm <- function(x, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  num.plot.fit <- length(names(x$Fits))-1
  num.plot.F <- length(names(x$MarginalTest))-1
  num.simul <- dim(x$MarginalTests)[1]-1

  figs.to.plot <- num.plot.fit + num.plot.F
  nb.fig <- prod(par("mfcol"))
  par(ask = interactive() && nb.fig < figs.to.plot)

  par(las = 1)
  for(i in 1:num.plot.fit) {
    hist(x$Fits[, i + 1], xlab = names(x$Fits)[i + 1], main = "Histogram of model terms")
    abline(v = x$Fits[1, i + 1])}
  for(j in 1:num.plot.F) {
    p.value <-
        (sum(x$MarginalTest[-1, j + 1] > x$MarginalTest[1, j + 1]) + 1) /
            (num.simul + 1)
    hist(x$MarginalTest[, j + 1],
         xlab = paste(names(x$MarginalTest)[j + 1], " (p = ",
         format.pval(p.value), ")", sep = ""),
         main = "Histogram of marginal F-values")
    abline(v = x$MarginalTest[1, j + 1])}
}
  


lm.phylog <- function(formula, data, max.num = 0, weights = NULL,
                      exclude.tips = NULL, lapply.size = 100){
    ## using looping over lapply and keeping my.drop outside  
    
    my.drop <- function (object) 
    {
## Returns ONLY the marginal F's.
## This function is based on drop1 (a function in the R base distribution),
## by B.D. Ripley;
## all I have done here is eliminate things I didn't need.

        ## From R-1.9.0 deviance.lm no longer available. Kludge here:
        deviance.lm <- function(object, ...)
            sum(weighted.residuals(object)^2, na.rm=TRUE)
        
        x <- model.matrix(object)
        iswt <- !is.null(wt <- object$weights)
        n <- nrow(x)
        asgn <- attr(x, "assign")
        tl <- attr(object$terms, "term.labels")
        scope <- drop.scope(object)
        ndrop <- match(scope, tl)
        ns <- length(scope)
        rdf <- object$df.resid
        chisq <- deviance.lm(object)
        dfs <- numeric(ns)
        RSS <- numeric(ns)
        y <- object$residuals + predict(object)
        rank <- object$rank
        for (i in 1:ns) {
            ii <- seq(along = asgn)[asgn == ndrop[i]]
            jj <- setdiff(seq(ncol(x)), ii)        
            z <- if (iswt) 
                lm.wfit(x[, jj, drop = FALSE], y, wt)
            else lm.fit(x[, jj, drop = FALSE], y)
            dfs[i] <- z$rank
            RSS[i] <- deviance.lm(z)
        }
        ##   scope <- c("<none>", scope)
        dfs <- c(object$rank, dfs)
        RSS <- c(chisq, RSS)
        dfs <- (dfs[1] - dfs)[-1]
        rdf <- object$df.resid
        dev <- RSS[-1]-RSS[1]
        rms <- RSS[1]/rdf
        Fs <- (dev/dfs)/rms
        Fs
    }
       
    
    if(min(data$sim.counter) == 1)
        stop("You need the input file to include original data from an inp or similar file; these are the data whose coefficients you are trying to test")

    ## I do not use match.call, and evaluating it, to avoid creating a new mf,
    ##  which could be huge.
    ## weights, to make things comparable to other R functions,
    ## can either be a column
    ## in the data set data or something else in the parent environment.
    ## First look for it in "data"; if it is not there, go to parent environment
    ## Note, though, that if it is a vector in the parent environment it MUST have
    ## the same length as any column in data, before data is reduced with arguments
    ## max.num and exclude.tips.
    ## Things seem to work when weights is called as a function of a
    ## variable in data, 
    ## but the data frame needs to be given explicitly;
    ## ej: weights = sqrt(dataxx$variableyy)
    
    ## Note that there have been several changes, w.r.t.
    ## initial versions, because of
    ## scoping differences between R-1.1.1 and R-1.2.0; I have used a solution
    ## provided by P. Dalgaard (email of 12-Oct-2000)
    
    
    if(!is.null(substitute(weights))){
        if(match(deparse(substitute(weights)), names(data), nomatch = 0))
            data$w <- eval(parse(text = paste(deparse(substitute(data)),
                                 "$", substitute(weights), sep = "")))
        else {if(length(weights) != length(data$sim.counter))
                  stop("Length of weights not equal to number of rows")
              data$w <- weights}
        if(any(data$w <= 0)) stop("Weights have to be positive") }
    else data$w <- NULL
    
    if(!is.null(exclude.tips)) data <-
        data[match(data$Tips, exclude.tips, nomatch = 0) == 0, ] 
    if(max.num) {
        max.num <- min(max(data$sim.counter), max.num)
        data <- data[data$sim.counter<max.num + 1, ]
    }
    else max.num <- max(data$sim.counter)
    
    
    names.vars <- drop.scope(formula)
    ##  terms.in.model <- names(lm(formula = formula, data = data, weights = w, subset = sim.counter == 0)[[1]])
    ## there should be a simpler way
    
    ## terms.in.model <- names(lm(formula = formula, data = data,
    ##                            subset = sim.counter == 0)[[1]])

    ## I am getting a "note": no visible binding for global variable ‘sim.counter’
    ## let's see if this solves it
    terms.in.model <- names(lm(formula = formula,
                               data = data[data$sim.counter == 0,])[[1]])


    ##note that there is no weights argument here; problems getting it to run
    ## in both R-1.1.1 and R-1.2.0; after all, the names are the same with or
    ##without weights. Still, there should be a simpler way!!
    
    
    
    
    if(terms.in.model[1] == "(Intercept)") terms.in.model[1] <- "Intercept"
    
    loop.counter <- (max.num + 1)%/%lapply.size
    rest.of.data <- (max.num + 1)%%lapply.size
    
    tmp <- matrix(nrow = max.num + 1,
                  ncol = length(names.vars) + length(terms.in.model) + 1)
    i <- 0
    if (loop.counter) { #only enter in the loop if needed
        for(i in 1:loop.counter){
            datai <-
                data[data$sim.counter <= ((i * lapply.size)-1) &
                     data$sim.counter >= ((i-1) * lapply.size), ]
            
            ## obtain output of interest ---fm[[1]], my.drop(fm)---
            ##by applying that function
            ## over the subset of data within loop; as result is list,
            ##unlist and turn into a matrix
            
            tmp[(((i-1) * lapply.size) + 1):(i * lapply.size), ] <-
                matrix(unlist(lapply(split(datai, datai$sim.counter), 
                                     function(datos, formula){
                                         environment(formula) <- environment()    
                                         fm <- lm(formula = formula,
                                                  data = datos, weights = datos$w);
                                         c(fm$coefficients, summary(fm)$fstatistic[[1]],
                                           my.drop(fm))
                                     }, 
                                     formula = formula)), 
                       nrow = lapply.size, byrow = TRUE)
        }
    }
    
    if (rest.of.data){
        datai <- data[data$sim.counter >= (loop.counter * lapply.size), ]
        tmp[(((i * lapply.size) + 1):(max.num + 1)), ] <-
            matrix(unlist(lapply(split(datai, datai$sim.counter), 
                                 function(datos, formula){
                                     environment(formula) <- environment()  
                                     fm <- lm(formula = formula,
                                              data = datos, weights = datos$w);
                                     c(fm$coefficients, summary(fm)$fstatistic[[1]],
                                       my.drop(fm))
                                 }, 
                                 formula = formula)), 
                   nrow = rest.of.data, byrow = TRUE)
    }
    
    structure(list(call = match.call(), 
                   Fits = data.frame(structure(cbind(seq(from = 0, to = max.num), 
                   tmp[, c(1:length(terms.in.model))]),
                   dimnames = list(NULL, c("sim.counter", terms.in.model)))), 
                   MarginalTests =
                   data.frame(structure
                              (cbind(seq(from = 0, to = max.num), 
                                     tmp[, -c(1:length(terms.in.model))]),
                               dimnames = list(NULL, c("sim.counter",
                               "Overall F", names.vars))))), 
              class = c("phylog.lm", "list"))
    
}




prcomp.phylog <- function(x, max.num = 0, exclude.tips = NULL,
                          lapply.size = 100, center = TRUE, scale = TRUE, ...){
##  require(mva, quietly = FALSE, warn.conflicts = TRUE) no longer needed
## pass as x the columns you want analyzed AND, as first column, sim.counter
    ## and as second columns Tips
## using looping over lapply and keeping my.drop outside  
  if(min(x$sim.counter) == 1)
      stop("You need the input file to include original x from an inp or similar file; these are the data whose coefficients you are trying to test")

  if(!is.null(exclude.tips)) x <-
      x[match(x$Tips, exclude.tips, nomatch = 0) == 0, ] 
  if(max.num) {max.num <- min(max(x$sim.counter), max.num)
               x <- x[x$sim.counter<max.num + 1, ] }
  else max.num <- max(x$sim.counter)

  x <- x[, -2]
  
  loop.counter <- (max.num + 1)%/%lapply.size
  rest.of.data <- (max.num + 1)%%lapply.size
  
  tmp <- matrix(nrow = max.num + 1, ncol = dim(x)[2]-1)
  
  
  i <- 0
  if (loop.counter) { #only enter in the loop if needed
      for(i in 1:loop.counter){
          xi <-
              x[x$sim.counter <= ((i * lapply.size)-1) &
                   x$sim.counter >= ((i-1) * lapply.size), ]  
          tmp[(((i-1) * lapply.size) + 1):(i * lapply.size), ] <- 
              matrix(unlist(lapply(split(xi, xi$sim.counter),
                                   function(datos, center, scale){
                                       (prcomp(datos[, -1], center = center,
                                               scale = scale)[[1]])^2
                                   },
                                   center = center, scale = scale)),
                     nrow = lapply.size, byrow = TRUE)
      }
  }

  if (rest.of.data){
      xi <- x[x$sim.counter >= (loop.counter * lapply.size), ]
      tmp[(((i * lapply.size) + 1):(max.num + 1)), ] <-
          matrix(unlist(lapply(split(xi, xi$sim.counter),
                               function(datos, center, scale){
                                   (prcomp(datos[, -1], center = center,
                                           scale = scale)[[1]])^2
                               },
                               center = center, scale = scale)),
                 nrow = rest.of.data, byrow = TRUE)
  }


  dimnames(tmp) <- list(NULL, paste("lambda", 1:(dim(x)[2]-1), sep = ""))  
  
  structure(list(call = match.call(), 
                 Eigenvalues =
                 data.frame("sim.counter" = seq(from = 0, to = max.num), tmp)), 
            class = c("phylog.prcomp", "list"))
  
}

summary.phylog.prcomp <- function(object, ...){

## note: we use quantile, so there is linear interpolation to find quantiles
##if(class(object) != "phylog.lm")
##stop("Object mas be of class phylog.lm, such as obtained from the lm.phylog function")
  num.simul <- dim(object$Eigenvalues)[1]-1
  eigenvalue.sim <- as.matrix(object$Eigenvalues[-1, ][-1])
  eigenvalue.original <- object$Eigenvalues[1, ][-1]
  # make that into a matrix by repeating the value as many times as needed
  # takes up space but seems more efficient than other solutions
  eigenvalue.original <-
      matrix(rep(as.matrix(eigenvalue.original), num.simul),
             nrow = num.simul, byrow = TRUE)
  

## this is what p-value does: count the number of cases where the simulated eigenvalue is
## larger than the one from inp, add one, and divide by num.simulations + 1.
  p.value <- (apply(eigenvalue.sim>eigenvalue.original, 2, sum) + 1)/(num.simul + 1)

## take into account that, to look whetehr the second eigenvalue
## is "significant" the first must have been to.
## thus, you get support agains the null hypothesis (that, say, lambda 2
## is really there) if BOTH the first and the second eigenvalues
## from the inp data are larger than most of the simulated ones.
## Thus, the p-value is computed as the number of simulations
## in which the second simulated eigenvalue is larger than the observed
## or the first simulated is larger than the observed, or both.
## In other words, the only cases that contribute to provide
## support against the null are simulations in which
## BOTH the first and the second eigenvalues are smaller
## than the observed ones.

## In this code, I first look at wheter (cell-wise, as above)
## the simulated is larger than the observed.
## Then, I go row by row and compute the cumulative sum
## so if I have FALSE TRUE that is a 0, 1, if I have FALSE FASLE that is 0, 0, 
## if I have TRUE, FALSE that is 1, 1, etc.
## But I don't care about the value of the sum, only whether or not it
## is larger than 0.
## Once in each cell I have a TRUE or FALSE again, 
## I add over columns to get the p-values.
  p.value2 <- (apply(t(apply(eigenvalue.sim>eigenvalue.original, 1,
                             function(x){as.logical(cumsum(x))})),
                     2, sum) + 1)/(num.simul + 1)
  names(p.value2) <- names(p.value)
  mean.eigenvalue <- apply(eigenvalue.sim, 2, mean)
  quant.eigenvalue <- apply(eigenvalue.sim, 2, quantile,
                            probs = c(.5, .9, .95, .99, .999))
  structure(list(call = object$call,
                 original.eigenvalues = object$Eigenvalues[1, ][-1], 
                 p.value.component = p.value, p.value.multiple = p.value2,
                 mean.eigenvalue = mean.eigenvalue,
                 quant.eigenvalue = quant.eigenvalue,
                 num.simul = num.simul),
            class = "summary.phylog.prcomp")
}


print.summary.phylog.prcomp <- function(x, ...){
  cat("\nCall: \n")
  print(x$call)
  cat("\nEigenvalues from inp data:\n")
  print(unlist(x$original.eigenvalue))
  cat("\n'Component-wise' P-value for eigenvalues:\n")
  print(x$p.value.component)
  cat("\n'Multiple eigenvalue' P-value:\n")
  print(x$p.value.multiple)
  cat("\nMean of simulated eigenvalues: \n")
  print(x$mean.eigenvalue)
  cat("\nQuantiles of eigenvalues' distributions:\n")
  print(x$quant.eigenvalue)
  cat("\nNumber of simulations used in analyses: ",
      format(x$num.simul, ...), "\n")
  invisible(x)}


plot.phylog.prcomp <- function(x, ...){
  oldpar <- par(no.readonly = TRUE, ...)
  on.exit(par(oldpar))
  num.plots <- dim(x$Eigenvalues)[2]-1
  num.simul <- dim(x$Eigenvalues)[1]-1

  nb.fig <- prod(par("mfcol"))
  par(ask = interactive() && nb.fig < num.plots)

  par(las = 1)
  mean.eigenvalue <- apply(as.matrix(x$Eigenvalues[-1, ][-1]), 2, mean)
  plot(c(1, num.plots), c(0, max(x$Eigenvalues[, -1])),
       type = "n", ylab = "Value of eigenvalue",
       xlab = "Eigenvalue", main = "Scree plot as in Horn (1965)")
  lines(c(1:num.plots), x$Eigenvalues[1, -1], lty = 1)
  lines(c(1:num.plots), mean.eigenvalue, lty = 2)
  legend(num.plots/2, max(x$Eigenvalues[, -1]), lty = c(1, 2),
         c("Observed", "Simulated"), bty = "n", xjust = 0)

  for(i in 1:num.plots) {
    p.value <- (sum(x$Eigenvalues[-1, i + 1] > x$Eigenvalues[1, i + 1]) + 1) /
        (num.simul + 1)
    hist(x$Eigenvalues[, i + 1],
         xlab = paste(names(x$Eigenvalues)[i + 1], " (p = ",
         format.pval(p.value), ")", sep = ""),
         main = "Histogram of eigenvalues")
    
##hist(x$Eigenvalues[, i + 1], xlab = names(x$Eigenvalues)[i + 1], main = "Histogram of eigenvalues")
    abline(v = x$Eigenvalues[1, i + 1])}
}


cancor.phylog <- function(data1, data2, max.num = 0, exclude.tips = NULL,
                          lapply.size = 100,
                          xcenter = TRUE, ycenter = TRUE){

## pass as data the columns you want analyzed AND, as first column, sim.counter
## and as second columns Tips
## using looping over lapply and keeping my.drop outside  
  if(min(data1$sim.counter) == 1)
      stop("You need the input file to include original data from an inp or similar file; these are the data whose coefficients you are trying to test")

  if(!is.null(exclude.tips)) {
      data1 <- data1[match(data1$Tips, exclude.tips, nomatch = 0) == 0, ]
      data2 <- data2[match(data2$Tips, exclude.tips, nomatch = 0) == 0, ]
  }
  if(max.num) {
      max.num <- min(max(data1$sim.counter), max.num)
      data1 <- data1[data1$sim.counter<max.num + 1, ]
      data2 <- data2[data2$sim.counter<max.num + 1, ]
  }
  else max.num <- max(data1$sim.counter)
  
  if(!all(data1$Tips == data2$Tips))
      stop("Different Tips in the two data sets")
  if(max(data1$sim.counter) != max(data2$sim.counter))
      stop("Different number of simulations in the two data sets")

  vars.x1 <- dim(data1)[2]-2
  vars.x2 <- dim(data2)[2]-2
  
  data <- cbind(data1[, -2], data2[, -c(1, 2)])
  number.of.tips <- length(data1$Tips[data1$sim.counter == 0])
  rm(data1, data2)

  loop.counter <- (max.num + 1)%/%lapply.size
  rest.of.data <- (max.num + 1)%%lapply.size

  tmp <- matrix(nrow = max.num + 1, ncol = min(vars.x1, vars.x2))
  
  i <- 0
  if (loop.counter) { #only enter in the loop if needed
      for(i in 1:loop.counter){
          datai <-
              data[data$sim.counter <= ((i * lapply.size)-1) &
                   data$sim.counter >= ((i-1) * lapply.size), ]  
          tmp[(((i-1) * lapply.size) + 1):(i * lapply.size), ] <- 
              matrix(unlist(lapply(split(datai, datai$sim.counter),
                                   function(datos, xcenter, ycenter){
                                       cancor(datos[, 2:(vars.x1 + 1)],
                                              datos[, -(1:vars.x1 + 1)],
                                              xcenter = xcenter,
                                              ycenter = ycenter)[[1]]
                                   },
                                   xcenter = xcenter,
                                   ycenter = ycenter)),
                     nrow = lapply.size, byrow = TRUE)
      }
  }
  
  if (rest.of.data){
      datai <-
          data[data$sim.counter >= (loop.counter * lapply.size), ]
      tmp[(((i * lapply.size) + 1):(max.num + 1)), ] <-
          matrix(unlist(lapply(split(datai, datai$sim.counter),
                               function(datos, xcenter, ycenter){
                                   cancor(datos[, 2:(vars.x1 + 1)],
                                          datos[, -(1:vars.x1 + 1)],
                                          xcenter = xcenter,
                                          ycenter = ycenter)[[1]]
                               },
                               xcenter = xcenter,
                               ycenter = ycenter)),
                 nrow = rest.of.data, byrow = TRUE)
  }
  
  
  dimnames(tmp) <- list(NULL, paste("corr", 1:(dim(tmp)[2]), sep = ""))  
  
### for likelihood-ratio based test
 ## the LR.stat is from expression 14.23 in Krzanowski (p. 447)
  n.prime <- number.of.tips - 0.5 * (vars.x1 + vars.x2 + 3)
  ## for compacteness don't create an intermediate object to hold the LR statistic.
                                 
  structure(list(call = match.call(), 
                 CanonicalCorrelations =
                 data.frame("sim.counter" = seq(from = 0, to = max.num), tmp), 
                 LR.statistic =
                 data.frame("sim.counter" = seq(from = 0, to = max.num), 
                            lambda =
                            apply(tmp, 1,
                                  function(y) -n.prime * sum(log(1 - (y^2)))))), 
            class = c("phylog.cancor", "list"))
}

summary.phylog.cancor <- function(object, ...){
    
## note: we use quantile, so there is linear interpolation to find quantiles

  num.simul <- dim(object$CanonicalCorrelations)[1]-1
  cancor.sim <- as.matrix(object$CanonicalCorrelations[-1, ][-1])
  cancor.original <- object$CanonicalCorrelations[1, ][-1]
  overall.sim <- object$LR.statistic[-1, 2]
##browser()
  overall.original <- object$LR.statistic[1, 2]

  ## make that into a matrix by repeating the value as many times as needed
  ## takes up space but seems more efficient than other solutions
  cancor.original <-
      matrix(rep(as.matrix(cancor.original), num.simul), nrow = num.simul,
             byrow = TRUE)

## this is what p-value does: count the number of cases where the
  ##simulated eigenvalue is
## larger than the one from inp, add one, and divide by num.simulations + 1.

  p.value <- (apply(cancor.sim >= cancor.original, 2, sum) + 1)/(num.simul + 1)

  ## the following is the same idea, but just for the LR overall statistic of
  # of hypothesis that all canonical corr. are zero
  overall.p.value <- (sum(overall.sim >= overall.original) + 1)/(num.simul + 1)
  names(overall.p.value) <- "P-value"

## see code for phylog.prcomp

  p.value2 <-
      (apply(t(apply(cancor.sim>cancor.original, 1,
                     function(x){as.logical(cumsum(x))})), 2, sum) + 1)/
                         (num.simul + 1)
  names(p.value2) <- names(p.value)
  
  names(overall.original) <- "lambda"
  
  quant.cancor <- apply(cancor.sim, 2, quantile,
                        probs = c(.5, .9, .95, .99, .999))
  structure(list(call = object$call,
                 original.LR.statistic = overall.original, 
                 original.canonicalcorrelations =
                 object$CanonicalCorrelations[1, ][-1], 
                 p.value.overall.test = overall.p.value, 
                 p.value.corwise = p.value, p.value.multiple = p.value2,
                 quant.canonicalcorrelations = quant.cancor,
                 num.simul = num.simul), class = "summary.phylog.cancor")
}


print.summary.phylog.cancor <- function(x, ...){
    cat("\nCall: \n")
    print(x$call)
    cat("\nCanonical correlations from original data:\n")
    print(unlist(x$original.canonicalcorrelations))
    cat("\nLR statistic from original data to test that all canonical correlations are zero:\n")
    print(x$original.LR.statistic)
    cat("\nTest that all canonical correlations are zero:\n")
    print(x$p.value.overall.test)
    cat("\n'Correlation-wise' P-value for canonical correlations:\n")
    print(x$p.value.corwise)
    cat("\n'Multiple' P-value for canonical correlations:\n")
    print(x$p.value.multiple)  
    cat("\nQuantiles of canonical correlations' distributions:\n")
    print(x$quant.canonicalcorrelations)
    cat("\nNumber of simulations used in analyses: ", format(x$num.simul, ...), "\n")
    invisible(x)
}



plot.phylog.cancor <- function(x, ...){
  oldpar <- par(no.readonly = TRUE, ...)
  on.exit(par(oldpar))
  num.plots <- dim(x$CanonicalCorrelations)[2]-1
  num.simul <- dim(x$CanonicalCorrelations)[1]-1

  nb.fig <- prod(par("mfcol"))
  par(ask = interactive() && nb.fig < num.plots)

  par(las = 1)
  for(i in 1:num.plots) {
      p.value <- (sum(x$CanonicalCorrelations[-1, i + 1] >
                      x$CanonicalCorrelations[1, i + 1]) + 1) /
                          (num.simul + 1)
      hist(x$CanonicalCorrelations[, i + 1],
           xlab = paste(names(x$CanonicalCorrelations)[i + 1],
           " (p = ", format.pval(p.value), ")", sep = ""),
           main = "Histogram of canonical correlations")
    
## hist(x$CanonicalCorrelations[, i + 1],
      ##xlab = names(x$CanonicalCorrelations)[i + 1], main = "Histogram of canonical correlations")
      abline(v = x$CanonicalCorrelations[1, i + 1])}
}
  

read.phylog.matrix <- function(x){read.table(x, skip = 1, header = TRUE)}

matrix.D <- function(x){
  ## To obtain the matrix D in p. 361 of Garland & Ives, 2000.
  ## you don't need to call it directly. This function is called by other functions.
  
  ## In the terminology of Garland & Ives, 2000 we want D %*% C %*% D' = I
  ## We can follow Schott, 1997, "Matrix analysis for statistics", pp. 139 and 141.
  ## Let D = inverse of T.
  ## Then T is any square root matrix of C, s.t., T %*% T' = C.
  ## T can be obtained with an eigenvalue-eigenvector (spectral) decomposition
  ## where: T = P %*% sqrt(lambda) %*% t(P),
    ## with P eigenvector of C and lmbda the diagonal matrix of eigenvalues.

  eigx <- eigen(x)
  
  T <- eigx$vectors  %*%  diag(sqrt(eigx$values))  %*%  t(eigx$vectors)
  matrixD <- solve(T)
  dimnames(matrixD) <- dimnames(x)
  matrixD
}



phylog.gls.fit <- function(x, y, cov.matrix, intercept = TRUE, exclude.tips = NULL){
  # to fit a linear model using gls as in Garland & Ives, 2000.
  # the data have to be continuous variables. If you have a
  # factor, convert it to the appropriate dummy variables (you
  # can use function "contrasts").
  
  mD <- matrix.D(cov.matrix)
  y <- mD %*% y
  if(is.vector(x)) {
    tmp.name <- deparse(substitute(x))
    x <- as.data.frame(x)
    names(x) <- tmp.name
    rm(tmp.name)}
    else x <- as.data.frame(x)
  dimnames(x)[[1]] <- dimnames(y)[[1]] # to allow to select tips
  if(intercept){x1 <- mD %*% cbind(rep(1, length(y)), as.matrix(x))
                if(!is.null(exclude.tips)) {
                    x1 <- x1[match(dimnames(x1)[[1]], exclude.tips,
                                   nomatch = 0) == 0, ]
                  y <- y[match(dimnames(y)[[1]], exclude.tips,
                               nomatch = 0) == 0, ]
                }
                fit1 <- lm(y ~ x1-1)
                names(fit1$coefficients) <- c("(Intercept)", names(x))
               }
  else {x1 <- mD %*% as.matrix(x)
        if(!is.null(exclude.tips)) {
                  x1 <- x1[match(dimnames(x1)[[1]], exclude.tips,
                                 nomatch = 0) == 0, ]
                  y <- y[match(dimnames(y)[[1]], exclude.tips,
                               nomatch = 0) == 0, ]
                }
        fit1 <- lm(y ~ x1-1)
        names(fit1$coefficients) <- names(x)  
    } 
  fit1
}


cor.origin <- function(x, y){lm(y ~ x-1)$coefficients * sqrt(x %*% x)/sqrt(y %*% y)}
