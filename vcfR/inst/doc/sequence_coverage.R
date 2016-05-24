## ------------------------------------------------------------------------
library(vcfR)
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

## ------------------------------------------------------------------------
vcf

## ------------------------------------------------------------------------
head(vcf)

## ------------------------------------------------------------------------
strwrap(grep('DP', vcf@meta, value = TRUE))

## ------------------------------------------------------------------------
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

## ---- fig.align='center', fig.width=7, fig.height=7----------------------
par(mar=c(8,4,1,1))
#boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", log='y', las=2)
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)
abline(h=seq(0,1e4, by=100), col="#C0C0C088")
par(mar=c(5,4,4,2))

## ---- fig.align='center', fig.width=7------------------------------------

if( require(reshape2) & require(ggplot2) ){
  dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
  dpf <- dpf[ dpf$Depth > 0,]
  p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                         scale = "count", trim=TRUE)
  p <- p + theme_bw()
  p <- p + theme(axis.title.x = element_blank(), 
                 axis.text.x = element_text(angle = 60, hjust = 1))
#  p <- p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black")
  p <- p + scale_y_continuous(trans=scales::log2_trans(), breaks=c(1, 10, 100, 800))
  p
} else {
  message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed.  Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
}


## ------------------------------------------------------------------------
sums <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[1,])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[2,])
dp[dp2 > 0] <- NA
dp[dp < 4] <- NA

## ---- fig.align='center', fig.width=7, fig.height=7----------------------
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth")
abline(h=seq(0,200, by=20), col="#C0C0C088")
par(mar=c(5,4,4,2))

