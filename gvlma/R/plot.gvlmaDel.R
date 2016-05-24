"plot.gvlmaDel" <-
function(x, which = 1:2, TukeyStyle = TRUE,
         ask = prod(par("mfcol")) <  max(c(10,5)[which]) && dev.interactive(),
         pointlabels, ...
         )
{
  gvlmaDelobj <- x
  if (missing(pointlabels)) pointlabels <- rownames(gvlmaDelobj)
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  timeseq <- attr(gvlmaDelobj, "timeseq")
  if (show[1]) # plot deleted statistics and p-values versus time sequence
    {
      for (w in seq(1, 10, by = 2))
        {
          statname <- names(gvlmaDelobj)[w]
          statnm <-
            switch(statname,
                   "DeltaGlobalStat" =
                   expression(paste("Deleted ", {G[4]}^2, " statistic (% change)")),
                   "DeltaStat1" =
                   expression(paste("Deleted ", {S[1]}^2, " statistic (% change)")),
                   "DeltaStat2" =
                   expression(paste("Deleted ", {S[2]}^2, " statistic (% change)")),
                   "DeltaStat3" =
                   expression(paste("Deleted ", {S[3]}^2, " statistic (% change)")),
                   "DeltaStat4" =
                   expression(paste("Deleted ", {S[4]}^2, " statistic (% change)"))
                   )
          pvalname <- names(gvlmaDelobj[w+1])
          pvalnm <- switch(pvalname,
                           "GStatpvalue" =
                           expression(paste("Deleted ", {G[4]}^2, " p-value")),
                           "Stat1pvalue" =
                           expression(paste("Deleted ", {S[1]}^2, " p-value")),                           
                           "Stat2pvalue" =
                           expression(paste("Deleted ", {S[2]}^2, " p-value")),                           
                           "Stat3pvalue" =
                           expression(paste("Deleted ", {S[3]}^2, " p-value")),                           
                           "Stat4pvalue" =
                           expression(paste("Deleted ", {S[4]}^2, " p-value")),                           
                           )
          plot(timeseq, gvlmaDelobj[,w], xlab = "Time sequence",
               ylab = statnm)
          plot(timeseq, gvlmaDelobj[, w+1], xlab = "Time sequence",
               ylab = pvalnm, ylim = c(0,1))
        }
    }
  if (show[2])
    {
      for (w in seq(1, 10, by = 2))
        {
          nm <- names(gvlmaDelobj)[w]
          nmm <- switch(nm,
                        "DeltaGlobalStat" = "G",
                        "DeltaStat1" = "S1",
                        "DeltaStat2" = "S2",
                        "DeltaStat3" = "S3",
                        "DeltaStat4" = "S4")
          display.delstats(gvlmaDelobj[,w], gvlmaDelobj[,w+1],
                           nsd = 3,
                           TukeyStyle, statname = nmm, pointlabels
                           )
        }
    } 
}

