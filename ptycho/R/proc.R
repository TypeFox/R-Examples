indicVarCols <- function(obj) {
  cols <- ptycho.cols(obj)
  grep("^indic.var", cols)
}

indicGrpCols <- function(obj) {
  cols <- ptycho.cols(obj)
  grep("^indic.grp", cols)
}

ptycho.cols <- function(z) {
  if (inherits(z, "mcmc.list")) {
    colnames(z[[1]])
  } else if (is.numeric(z)) {
    if (is.null(dim(z))) names(z) else colnames(z)
  } else {
    stop("ptycho.cols: invalid input")
  }
}

# Don't print attributes by default.  Probably won't be used very often because
# there will usually be so many columns in a ptycho object that one doesn't want
# to print all of them ... But it is doubly annoying to see all the attributes
# when I forget to type foo[,1:5].
print.ptycho <- function(x, ...) { print(x[,], ...) }

parsePtychoFilename <- function(fn) {
  # parse out rpl and col
  # value returned by strsplit() will be list of length 1 with entry
  # c("",rpl,col)
  ns <- as.numeric(strsplit(basename(fn), "[^[:digit:]]+")[[1]][-1])
}

Evar <- function(eta2, tau, n, A, B, V, s) {
  z <- tau * sqrt(n)
  z <- z * (B+V-s-1)/(A+s)
  z <- 1 + z * (1-eta2)^(n/2)
  z <- 1/z
}

getLastIters <- function(obj) {
  nmax <- max(obj[,"iter"])
  obj <- obj[ obj[,"iter"] == nmax, ]
}

### Compute means across chains.  Also, variance of tau.
meanTau <- function(obj) {
  as.numeric(mean(getLastIters(obj)[,"tau"]))
}
varTau <- function(obj) {
  z <- getLastIters(obj)[,c("tau","tau2")]
  z <- colMeans(z)
  as.numeric(z["tau2"] - z["tau"]^2)
}
meanIndicators <- function(obj) {
  colMeans(getLastIters(obj)[ , c(indicGrpCols(obj),indicVarCols(obj)) ])
}
meanVarIndicators <- function(obj) {
  colMeans(getLastIters(obj)[ , indicVarCols(obj) ])
}
meanGrpIndicators <- function(obj) {
  colMeans(getLastIters(obj)[ , indicGrpCols(obj) ])
}

### Compute the difference between chains.
drng <- function(x) { diff(range(x)) }
checkConvergence <- function(obj, doLastIterOnly=TRUE) {
  obj <- obj[ , -c(1,ncol(obj)) ]
  if (doLastIterOnly) obj <- getLastIters(obj)
  # A data frame cannot have [ or ] or , in its column names.  They all get
  # replaced with periods.  The original names are easier to work with.
  z <- ldply(strsplit(colnames(obj)[-1], "[],[]"),
             function(b) {
               if (b[1] == "tau") {
                 c(b,"1","")
               } else {
                 b[1] <- sub("indic.", "", b[1])
                 # Grouping variable for Across Traits
                 if (length(b) == 2) b[3] <- ""
                 b
               }
             })
  # Now back to the original task.
  obj <- ddply(data.frame(obj), c("iter"), colwise(drng))
  niter <- nrow(obj)
  obj <- melt(obj, id.vars=c("iter"), value.name="range")
  if (length(unique(obj$variable[seq_len(niter)])) > 1) {
    stop("checkConvergence: melt() behaves badly")
  }
  obj <- data.frame(iter=obj$iter,
                    type=rep(z$V1, each=niter),
                    index=rep(as.numeric(z$V2), each=niter),
                    y=rep(z$V3, each=niter),
                    range=obj$range)
  obj
}
