verifRankHist <-
function(forecasts, observations) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
#rank <- apply((forecasts > observations),1,sum)+1
 rank <- apply(cbind(observations,forecasts),1,function(x)
               rank(x,ties="random")[1] )
 k <- ncol(forecasts)
 hist(rank, breaks = 0:(k+1), prob = TRUE, xaxt = "n",
     xlab = "", ylab = "", main = "Verification Rank Histogram")
 axis(1, at = seq(.5, to = k+.5, by = 1), labels = 1:(k+1))
 abline(h=1/(k+1), lty=2)
 invisible(rank)
}

