# Generate a PSF chart

# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2014-09-06 18:51:58 gjw>
#
# Implement evaluate functionality.
#
# Copyright (c) 2009-2013 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

psfchart <- function(predicted,
                     actual,
                     bins=100,          # Number of bins to use for the plot.
                     threshold=0.5,     # The decision threshold.
                     splits=NULL,	# E.g., c(0.2, 0.8)
                     split.lables=c("Low", "Medium", "High"),
                     tic.size=0.2,      # proportional gap between axis ticks.
                     gg=TRUE,
                     verbose=FALSE)
{
  dosplits <- ! is.null(splits)

  if (is.factor(actual)) actual <- as.numeric(actual)-1
    
  doAggregate <- function()
  {
          
    # Bin the scores into "bins" bins and store the bins into
    # variable "bin". If there are more bins specified than
    # scores, then simply rank the scores. In the end we have for
    # each score a "rank", which is an integer in 1:bins or
    # 1:length(actual).
      
    if (length(actual) >= bins)
    {
      bin <- as.numeric(binning(predicted, bins, method="quantile",
                                ordered=FALSE, labels=FALSE))
    }
    else
    {
      bin <- reshape::rescaler(predicted, "rank")
    }

    # Check whether the pr and target agree.

    agree <- as.numeric(as.numeric(predicted > threshold) == actual)

    # Get the size of each bin (should all be the same +/-
    # 1). Also get a count of the positives in each bin, assumming
    # a 0/1 value for actual (so sum will work) and the accuracy
    # of each bin.

    agg <- aggregate(actual, list(bin), length)
    names(agg) <- c("bin", "size")
    agg$pos <- aggregate(actual, list(bin), sum)[[2]]
    agg$acc <- aggregate(agree, list(bin), sum)[[2]]
    agg$max <- aggregate(predicted, list(bin), max)[[2]]
    agg$tdiff <- agg$max - threshold
        
    # Rescale bins to be between 0 and 1 so AUC is more sensible.
  
    agg$rbin <- agg$bin/bins
  
    # Calulate the proportion accuracy
  
    agg$pacc <- agg$acc/agg$size

    return(agg)
  }
  
  # Determine the model's decision based on the score and threshold
  # and save that as pr.

  pr <- as.numeric(predicted > threshold)

  tp <- round(100 * sum(pr==1 & actual==1)/length(actual))
  fp <- round(100 * sum(pr==1 & actual==0)/length(actual))
  tn <- round(100 * sum(pr==0 & actual==0)/length(actual))
  fn <- round(100 * sum(pr==0 & actual==1)/length(actual))
    
  if (verbose)
    cat("\nData Summary:",
        sprintf("    Obs: %s\n",
                format(length(actual), big.mark=",")),
        sprintf("    Targets: %s; Rate: %0.2f%%\n",
                format(sum(actual==1), big.mark=","),
                100*sum(actual==1)/length(actual)),
        sprintf("    Model TP: %9s FN: %9s\n",
                format(sum(pr==1 & actual==1), big.mark=","),
                format(sum(pr==0 & actual==1), big.mark=",")),
        sprintf("          FP: %9s TN: %9s\n",
                format(sum(pr==1 & actual==0), big.mark=","),
                format(sum(pr==0 & actual==0), big.mark=",")))
  
  agg <- doAggregate()

    if (gg)
    {
        if (dosplits)
            classes <- data.frame(x=c(splits[1]/2,
                                      (splits[1]+splits[2])/2,
                                      (1+splits[2])/2),
                                  lbl=split.labels)
        
        quads <- data.frame(x=c(0, 1, 0, 1), hj=c(0, 1, 0, 1),
                            y=c(0, 0, 1, 1), vj=c(0, 0, 1, 1),
                            lbl=c(sprintf("True Negatives (%s%%)", tn),
                                sprintf("True Positives (%s%%)", tp),
                                sprintf("False Negatives (%s%%)", fn),
                                sprintf("False Positives (%s%%)", fp)))

        xthresh <- agg$rbin[which(abs(agg$tdiff) == min(abs(agg$tdiff)))][1]

        tics <- seq(0, 1, tic.size)
        ord <- order(predicted)
        scores <- data.frame(x=tics,
                             score=round(predicted[ord][c(1,
                                 round(tics*length(ord)))], 2))

        p <- ggplot2::ggplot(agg, ggplot2::aes(rbin, pacc))
        p <- p + ggplot2::geom_line()
        p <- p + ggplot2::ggtitle("Proportional Score Function (PSF) Curve")
        p <- p + ggplot2::scale_y_continuous("% Accuracy", limits=c(0,1),
                                    labels=100*tics, breaks=tics)
        p <- p + ggplot2::scale_x_continuous(paste("Proportion of Cases",
                                          "\nSorted by Increasing Risk Scores"),
                                    breaks=tics)
        p <- p + ggplot2::geom_text(data=quads,
                           ggplot2::aes(x=x, y=y, label=lbl, hjust=hj, vjust=vj, size=5))
        p <- p + ggplot2::geom_text(x=xthresh, ggplot2::aes(y=0, size=5),
                           label=sprintf("Threshold (%s)", threshold),
                           vjust=2, hjust=1.1)
        p <- p + ggplot2::geom_vline(xintercept=xthresh)
        p <- p + ggplot2::geom_text(data=scores, ggplot2::aes(x=x, y=1, label=score, size=5),
                           vjust=-0.5)
        p <- p + ggplot2::theme(legend.position="none")
        if (dosplits)
        {
            p <- p + ggplot2::geom_text(data=classes, ggplot2::aes(x=x, y=0.2, label=lbl))
            p <- p + ggplot2::geom_vline(xintercept=low, linetype="twodash", color="grey")
            p <- p + ggplot2::geom_vline(xintercept=high, linetype="twodash", color="grey")
        }
        return(p)
    }
    else
    {
        plot(agg$rbin, agg$pacc, type="l", xlim=c(0,1), ylim=c(0,1),
             xlab="Proportion of Cases\nSorted by Risk Score", ylab="% Accuracy")
        title(main="PSF\n")
  
        abline(v=0.25, lty=3)
        abline(v=0.75, lty=3)
        text(0.08, 0.08, "Low")
        text(0.5, 0.08, "Medium")
        text(0.92, 0.08, "High")
        
        # Add annotations for the sinlge plot.
        
        abline(v=agg$rbin[which(abs(agg$tdiff) == min(abs(agg$tdiff)))], lty=1)
        xthresh <- agg$rbin[which(abs(agg$tdiff) == min(abs(agg$tdiff)))]
        text(xthresh, 0.15, "Threshold", pos=2)
        text(xthresh, 0.1, threshold, pos=2)
        
        # TODO NEED TO PROGRAMMATICALLY DETERMINE THE LABELS FROM MIN SCORE TO MAX SCORE
        ord <- order(predicted)
        scores <- predicted[ord]
        axis(3, at=seq(0, 1, 0.2), padj=1.5, lwd.ticks=0,
             labels=round(scores[c(1, round(seq(0, 1, 0.2)*length(scores)))], 2))
        
        text(0.92, 1, "False Positives")
        text(0.085, 1, "False Negatives")
        text(0.92, 0, "True Positives")
        text(0.08, 0, "True Negatives")
        
        opar <- par(xpd=TRUE)
        text(0.5, 1.08,"Scores")
        par(opar)
    }
}
