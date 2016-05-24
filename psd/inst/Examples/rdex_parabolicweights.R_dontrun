\dontrun{#REX
library(psd)
library(grDevices)
library(RColorBrewer)

##
## Show parabolic weighting factors as a function of maximum tapers
##

# maximum number of tapers
maxx <- 1e3
# sequence in logspace
xseq <- seq(from=1,to=2.8,by=0.2)

# plot palette
pal <- "Spectral"
npal <- switch(pal, RdYlBu=11, Spectral=11, Blues=9)
pal.col <- RColorBrewer::brewer.pal(npal, pal)
cols <- rev(grDevices::colorRampPalette(pal.col)(maxx))

to_df <- function(W){
  # convert parabolic results to data.frame
  with(W, data.frame(taper_seq=as.vector(taper_seq), taper_weights=as.vector(taper_weights)))
}

## a roundabout way of bootstrapping y-axis limits:
#  upper
WgtsU <- parabolic_weights_fast(5)
DfU <- to_df(WgtsU)
#  lower
WgtsL <- parabolic_weights_fast(maxx)
DfL <- to_df(WgtsL)

ylims <- range(pretty(dB(c(DfL$taper_weights, DfU$taper_weights)))) + c(-2,5)

# function for plotting text
TFUN <- function(Df.){
  tx <- max(Df.$taper_seq)
  ty <- mean(Df.$taper_weights)
  text(log10(tx)+0.1, dB(ty), sprintf("%i", tx), col=cols[tx])
}

# function for weighting factors and plotting
WFUN <- function(x){
  message(x)
  Wgts <- parabolic_weights_fast(x)
  Df <- to_df(Wgts)
  lcol <- cols[x]
  lines(dB(taper_weights) ~ log10(taper_seq), Df, type="s", lwd=2, col=lcol)
  TFUN(Df)
}

## Plot parabolic weighting, in dB, colored by maximum num tapers
plot(dB(taper_weights) ~ log10(taper_seq), DfU, type="s", 
     xlim=c(0, log10(maxx)+0.2), 
     ylim=ylims, yaxs="i",
     col=cols[5], lwd=2,  
     main="Multitaper weighting factors by maximum tapers applied",
     xlab="log10 taper sequence", 
     ylab="dB")
TFUN(DfU)
invisible(lapply(round(10**xseq), FUN=WFUN))
WFUN(maxx)

##
}#REX
