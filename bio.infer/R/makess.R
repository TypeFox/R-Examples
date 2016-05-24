"makess" <-
function(bcnt,tname = "OTU", plothist=FALSE,prints=FALSE,nview=0) {

  names0 <- names(bcnt)
  site <- names0[1]
  abund <- names0[3]

  if (is.na(match(tname, names0))) {
    stop(paste("No field name",tname,"in the benthic count file"))
  }

  if (is.factor(bcnt[, site])) {
    sitenames <- levels(bcnt[, site])[bcnt[, site]]
  }
  else {
    if (is.numeric(bcnt[, site])) {
      sitenames <- as.character(bcnt[, site])
    }
    else {
      sitenames <- bcnt[, site]
    }
  }

  sitenames.u <- unique(sitenames)
  numsites <- length(sitenames.u)

  if (is.factor(bcnt[, tname])) {
    taxa.all <- levels(bcnt[, tname])[bcnt[, tname]]
  }
  else {
    taxa.all <- bcnt[, tname]
  }
  taxa.u <- sort(unique(taxa.all))
  ntaxa <- length(taxa.u)

  abund.all <- bcnt[, abund]
  ss1 <- matrix(0, ncol = ntaxa, nrow = numsites)

  for (j in 1:ntaxa) {
    flush.console()
    incvec <- taxa.all == taxa.u[j]
    incvec[is.na(incvec)] <- FALSE
    sitenames.g <- sitenames[incvec]
    abund.g <- abund.all[incvec]
    abund.s <- tapply(abund.g, sitenames.g, sum)
    sitenames.s <- names(abund.s)
    for (k in 1:length(sitenames.s)) {
      ss1[match(sitenames.s[k], sitenames.u), j] <- abund.s[k]
    }
  }
  ss.df <- data.frame(sitenames.u, ss1)
  names(ss.df) <- c(site, taxa.u)

  totabund.t <- apply(ss1,1, sum)
  df1 <- data.frame(sitenames.u, totabund.t)
  totabund.r <- tapply(bcnt[, abund], bcnt[, site], sum, na.rm = TRUE)
  df2 <- data.frame(names(totabund.r), totabund.r)
  names(df2) <- c("sitenames.u", "totabund.r")

  df3 <- merge(df1, df2, by = "sitenames.u")
  prop.t <- df3$totabund.t/df3$totabund.r
  if( plothist ){
#    windows(width = 4, height = 4, pointsize = 10, restoreConsole = TRUE)
    dev.new()
    hist(prop.t, xlab = "Proportion of abundance",
         main = "OTU abundance/total abundance" )
  } else {
    dfhist <- hist( prop.t, plot=FALSE )
    dfhist$xname <- "Proportion of abundance"
    attr( ss.df, "histogram" ) <- dfhist
  }
  if( prints ){
    cat("Summary statistics for OTU abundance/total abundance\n")
    print(summary(prop.t))
    cat("Taxa lists at sites in which OTU abundance is low can be informative.\n")
  } else {
    attr( ss.df, "summary" ) <- summary(prop.t)
  }

  if (nview > 0) {
    df3 <- df3[order(prop.t), ]

    for (i in 1:nview) {
      if (is.factor(df3$sitenames.u)) {
        siteid <- levels(df3$sitenames.u)[df3$sitenames.u][i]
      }
      else {
        siteid <- df3$sitenames.u[i]
      }
      cat(siteid, "\n")
      incvec <- bcnt[, site] == siteid
      print(bcnt[incvec,])
    }
  }
  return(ss.df)
}
