### R code from vignette source 'NanoStringNorm_Introduction.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=40, signif=3, digits=3)
set.seed(0xdada)

## To create bitmap versions of plots with many dots, circumventing
##   Sweave's fig=TRUE mechanism...
##   (pdfs are too large)
openBitmap = function(nm, rows=1, cols=1) {
  png(paste("NSN-", nm, ".png", sep=""), 
       width=600*cols, height=700*rows, pointsize=14)
  par(mfrow=c(rows, cols), cex=2)
}


###################################################
### code chunk number 2: load.package
###################################################
require("NanoStringNorm")
data("NanoString")


###################################################
### code chunk number 3: eg.basic.usage (eval = FALSE)
###################################################
## require('NanoStringNorm');
## data("NanoString");


###################################################
### code chunk number 4: eg.basic.usage (eval = FALSE)
###################################################
## # example 1
## NanoString.mRNA.raw = NanoStringNorm(NanoString.mRNA);
## pdf('my_first_nsn_plot.pdf');
## Plot.NanoStringNorm(NanoString.mRNA.raw, plot.type = 'all');
## dev.off();


###################################################
### code chunk number 5: eg.read.xls
###################################################
# directly import the nCounter output
path.to.xls.file <- system.file("extdata", "RCC_files", "RCCCollector1_rat_tcdd.xls", 
                                package = "NanoStringNorm");
NanoString.mRNA <- read.xls.RCC(x = path.to.xls.file, sheet = 1);
# only keep the counts and not the header
NanoString.mRNA <- NanoString.mRNA$x;


###################################################
### code chunk number 6: eg.standard.1
###################################################
# specify housekeeping genes in annotation 
data(NanoString);
NanoString.mRNA[NanoString.mRNA$Name %in%  
	c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';


###################################################
### code chunk number 7: eg.standard.2
###################################################
# example 2
# normalize mRNA and output a matrix of normalized counts
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	anno = NA,
	CodeCount = 'sum',
	Background = 'mean',
	SampleContent = 'housekeeping.sum',
	round.values = FALSE,
	take.log = FALSE,
	return.matrix.of.endogenous.probes = TRUE
	);


###################################################
### code chunk number 8: eg.standard.3
###################################################
# example 3
# this is the recommended method!
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	CodeCount = 'geo.mean',
	Background = 'mean.2sd',
	SampleContent = 'housekeeping.geo.mean',
	round.values = TRUE,
	take.log = TRUE
	);


###################################################
### code chunk number 9: eg.standard.4
###################################################
# setup a binary trait using rat strain
sample.names <- names(NanoString.mRNA)[-c(1:3)];
strain1 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain1[grepl('HW',sample.names)] <- 2;
strain2 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain2[grepl('WW',sample.names)] <- 2;
strain3 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain3[grepl('LE',sample.names)] <- 2;
trait.strain <- data.frame(
	row.names = sample.names,
	strain1 = strain1,
	strain2 = strain2,
	strain3 = strain3
	);


###################################################
### code chunk number 10: eg.standard.5
###################################################
# You can also input the gene annotation separately to allow flexibility
NanoString.mRNA.anno <- NanoString.mRNA[,c(1:3)];
NanoString.mRNA.data <- NanoString.mRNA[,-c(1:3)];
#NanoString.mRNA.anno$cool.genes <- vector.of.cool.genes;


###################################################
### code chunk number 11: eg.standard.6
###################################################
# example 3
# include a trait for differential expression and batch effect evaluation
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA.data,
	anno = NanoString.mRNA.anno,
	CodeCount = 'geo.mean',
	Background = 'mean.2sd',
	SampleContent = 'top.geo.mean',
	round.values = TRUE,
	take.log = TRUE,
	traits = trait.strain
	);


###################################################
### code chunk number 12: eg.other.1 (eval = FALSE)
###################################################
## # z-value transformation.
## # scale each sample to have a mean 0 and sd
## # by default all the other normalization methods are 'none'
## # you cannot apply a log because there are negative values
## # good for meta-analysis and cross platform comparison abstraction of effect size
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	OtherNorm = 'zscore',
## 	return.matrix.of.endogenous.probes = TRUE
## 	);


###################################################
### code chunk number 13: eg.other.2 (eval = FALSE)
###################################################
## # inverse normal transformation.
## # use quantiles to transform each sample to the normal distribution
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	OtherNorm = 'rank.normal',
## 	return.matrix.of.endogenous.probes = TRUE
## 	);


###################################################
### code chunk number 14: eg.other.3 (eval = FALSE)
###################################################
## # quantile normalization.
## # create an empirical distribution based on the median gene counts at the same
## # rank across sample. then transform each sample to the empirical distribution.
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	OtherNorm = 'quantile',
## 	return.matrix.of.endogenous.probes = FALSE
## 	);


###################################################
### code chunk number 15: eg.other.4 (eval = FALSE)
###################################################
## # vsn.
## # apply a variance stabilizing normalization.
## # fit and predict the model using 'all' genes i.e. 'controls' and 
## # 'endogenous' (this is the default)
## # note this is just a wrapper for the vsn package
## # you could even add strata for the controls vs. the endogenous to review 
## # systematic differences
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	OtherNorm = 'vsn',
## 	return.matrix.of.endogenous.probes = FALSE,
## 	genes.to.fit = 'all',
## 	genes.to.predict = 'all'
## 	);


###################################################
### code chunk number 16: eg.other.5 (eval = FALSE)
###################################################
## # vsn.
## # this time generate the parameters (fit the model) on the 'controls' 
## # and apply (predict) on the 'endogenous'
## # alternatively you may want to use the endogenous probes for both fitting and predicting
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	OtherNorm = 'vsn',
## 	return.matrix.of.endogenous.probes = FALSE,
## 	genes.to.fit = 'controls',
## 	genes.to.predict = 'endogenous'
## 	);


###################################################
### code chunk number 17: eg.other.6 (eval = FALSE)
###################################################
## # vsn.
## # apply standard NanoString normalization strategies as an alternative to 
## # the vsn affine transformation.
## # this effectively applies the glog2 variance stabilizing transformation 
## # to the adjusted counts
## NanoString.mRNA.norm <- NanoStringNorm(
## 	x = NanoString.mRNA,
## 	CodeCount = 'sum',
## 	Background = 'mean',
## 	SampleContent = 'top.geo.mean',
## 	OtherNorm = 'vsn',
## 	return.matrix.of.endogenous.probes = FALSE,
## 	genes.to.fit = 'endogenous',
## 	genes.to.predict = 'endogenous',
## 	calib = 'none'
## 	);


###################################################
### code chunk number 18: eg.static.plots.0
###################################################
# plot all the plots for use in figures
pdf("NanoStringNorm_Example_Plots_%03d.pdf", onefile=FALSE);
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = 'all'
	);
dev.off();


###################################################
### code chunk number 19: eg.static.plots.1 (eval = FALSE)
###################################################
## # plot all the plots as PDF report
## pdf('NanoStringNorm_Example_Plots_All.pdf');
## Plot.NanoStringNorm(
## 	x = NanoString.mRNA.norm,
## 	label.best.guess = TRUE,
## 	plot.type = 'all'
## 	);
## dev.off();


###################################################
### code chunk number 20: eg.static.plots.2 (eval = FALSE)
###################################################
## # publication quality tiff volcano plot
## tiff('NanoStringNorm_Example_Plots_Volcano.tiff', units = 'in',  height = 6, 
## 	width = 6, compression = 'lzw', res = 1200, pointsize = 10);
## Plot.NanoStringNorm(
## 	x = NanoString.mRNA.norm,
## 	label.best.guess = TRUE,
## 	plot.type = c('volcano'),
## 	title = FALSE
## 	);
## dev.off();


###################################################
### code chunk number 21: eg.static.plots.3 (eval = FALSE)
###################################################
## # all plots as separate files output for a presentation
## png('NanoStringNorm_Example_Plots_%03d.png', units = 'in',  height = 6,
## 	width = 6, res = 250, pointsize = 10);
## Plot.NanoStringNorm(
## 	x = NanoString.mRNA.norm,
## 	label.best.guess = TRUE,
## 	plot.type = c('cv', 'mean.sd', 'RNA.estimates', 'volcano', 'missing',
## 		'norm.factors', 'positive.controls', 'batch.effects')
## 	);
## dev.off();


###################################################
### code chunk number 22: eg.static.plots.4 (eval = FALSE)
###################################################
## # user specified labelling with optimal resolution for most digital displays
## png('NanoStringNorm_Example_Plots_Normalization_Factors.png', units = 'in', height = 6,
## 	width = 6, res = 250, pointsize = 10);
## Plot.NanoStringNorm(
## 	x = NanoString.mRNA.norm,
## 	label.best.guess = FALSE,
## 	label.ids = list(genes = rownames(NanoString.mRNA.norm$gene.summary.stats.norm), 
## 		samples = rownames(NanoString.mRNA.norm$sample.summary.stats)),
## 	plot.type = c('norm.factors')
## 	);
## dev.off();


###################################################
### code chunk number 23: eg.interactive.plots.1 (eval = FALSE)
###################################################
## # plot the sample summaries to your browser
## Plot.NanoStringNorm.gvis(
## 	x = NanoString.mRNA.norm,
## 	plot.type = c('gene.norm', 'sample'),
## 	save.plot = FALSE
## 	);


###################################################
### code chunk number 24: eg.interactive.plots.2 (eval = FALSE)
###################################################
## # plot the gene summaries to a directory for distribution and later viewing
## # note.  if you save.plot = TRUE the default for path to mongoose is "web" i.e. it tries to download from online.
## # alternatively, you can specify a path to the binary or use "none" as below.
## Plot.NanoStringNorm.gvis(
## 	x = NanoString.mRNA.norm,
## 	plot.type = c('gene.norm', 'sample'),
## 	save.plot = TRUE,
## 	path.to.mongoose = 'none',
## 	output.directory = "NanoStringNorm_Interactive_Plot"
## 	);


###################################################
### code chunk number 25: eg.norm.comp (eval = FALSE)
###################################################
## # compare normalization methods using coefficient of variation 
## # of controls (CV) and correlation of replicates (ICC)
## 
## # specifiy housekeeping genes in annotation 
## NanoString.mRNA[NanoString.mRNA$Name %in% 
##   c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';
## 
## # strain x experimental condition i.e. replicate.
## # this is only a small subset of the original data  used for the plot
## biological.replicates <- c("HW_1.5_0","HW_1.5_0","HW_1.5_0","HW_1.5_100","HW_1.5_100","HW_1.5_100",
##                      "HW_6_100","HW_6_100","HW_3_100","HW_3_100","HW_3_100","HW_3_100",
##                      "LE_19_0","LE_19_0","LE_19_0","LE_96_0","LE_96_0","LE_96_0","HW_10_100",
##                      "HW_10_100","HW_10_100","HW_10_100","HW_6_100","HW_6_100","HW_96_0");
## 
## norm.comp.results.test <- norm.comp(
##   x = NanoString.mRNA,
##   replicates = biological.replicates,
##   CodeCount.methods = 'none',
##   Background.methods = 'none', 
##   SampleContent.methods = c('none','housekeeping.sum', 'housekeeping.geo.mean',
##                             'top.mean', 'top.geo.mean'), 
##   OtherNorm.methods = 'none',
##   verbose = FALSE
##   );


