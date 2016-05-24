## ----results="hide",echo=FALSE-------------------------------------------
cairo <- function(name, width, height, ...) grDevices::cairo_pdf(file = paste0(name, ".pdf"), width = width, height = height)
png <- function(name, width, height, ...) grDevices::png(file = paste0(name, ".png"), width = width*300, height = height*300)

## ----fig1plot,fig.width=8,fig.height=3-----------------------------------
library(tagcloud)
data(gambia)
tags <- strmultline(gambia$Term)[1:40]
weights <- -log(gambia$Pvalue)[1:40]
or <- gambia$OddsRatio[1:40]
colors <- smoothPalette(or, max=4)
tagcloud(tags, weights=weights, col=colors)

## ----fig2plot,fig.width=8,fig.height=12----------------------------------
par( mfrow=c( 3, 2 ) )
tagcloud(tags, weights=weights, col=colors, algorithm="oval")
tagcloud(tags, weights=weights, col=colors, algorithm="fill")
tagcloud(tags, weights=weights, col=colors, algorithm="snake")
tagcloud(tags, weights=weights, col=colors, algorithm="random")
tags2 <- gambia$Term[1:20]
cols2 <- colors[1:20]
wei2 <- weights[1:20]
tagcloud(tags2, weights=wei2, col=cols2, algorithm="list")
tagcloud(tags2, weights=wei2, col=cols2, algorithm="clist")

## ----fig3plot,fig.width=8,fig.height=4-----------------------------------
par(mfrow=c(1, 2))
tagcloud(tags, weights=weights, col=colors, fvert=0.3)
tagcloud(tags, weights=weights, col=colors, fvert=0.7)

## ----fig4plot,fig.width=8,fig.height=4-----------------------------------
par(mfrow=c(1, 2))
tagcloud(tags, weights=weights, col=colors, order="size")
tagcloud(tags, weights=weights, col=colors, order="random")

## ----fig5plot,eval=FALSE-------------------------------------------------
#  library(extrafont)
#  library(RColorBrewer)
#  fnames <- sample(fonts(), 40)
#  fweights <- rgamma(40, 1)
#  fcolors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( 40 )
#  tagcloud( fnames, weights=fweights, col=fcolors, family=fnames )

## ----fig6plot------------------------------------------------------------
library(RColorBrewer)
colors <- smoothPalette(weights, pal= brewer.pal( 11, "Spectral" ) )
tagcloud(tags, weights=weights, col=colors, order="size")

## ----fig7plot------------------------------------------------------------
palf <- colorRampPalette( c( "blue", "grey", "red" ) )
colors <- smoothPalette(weights, palfunc= palf )
tagcloud(tags, weights=weights, col=colors, order="size")

