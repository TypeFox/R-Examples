### R code from vignette source 'Using_RCircos.Rnw'

###################################################
### code chunk number 1: intro
###################################################
library(RCircos)


###################################################
### code chunk number 2: histogramdata
###################################################
data(RCircos.Histogram.Data)
head(RCircos.Histogram.Data)


###################################################
### code chunk number 3: heatmapdata
###################################################
data(RCircos.Heatmap.Data)
head(RCircos.Heatmap.Data)


###################################################
### code chunk number 4: linkdata
###################################################
data(RCircos.Link.Data)
head(RCircos.Link.Data)


###################################################
### code chunk number 5: Using_RCircos.Rnw:102-110
###################################################
data(UCSC.HG19.Human.CytoBandIdeogram);
head(UCSC.HG19.Human.CytoBandIdeogram);

data(UCSC.Mouse.GRCm38.CytoBandIdeogram);
head(UCSC.Mouse.GRCm38.CytoBandIdeogram);

data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
head(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);


###################################################
### code chunk number 6: Using_RCircos.Rnw:128-134
###################################################
chr.exclude <- NULL;
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude, 
       tracks.inside, tracks.outside);


###################################################
### code chunk number 7: Using_RCircos.Rnw:218-221
###################################################
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.cyto <- RCircos.Get.Plot.Ideogram();
rcircos.position <- RCircos.Get.Plot.Positions();


###################################################
### code chunk number 8: Using_RCircos.Rnw:224-225
###################################################
RCircos.List.Parameters()


###################################################
### code chunk number 9: Using_RCircos.Rnw:233-236
###################################################
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 30000;
RCircos.Reset.Plot.Parameters(rcircos.params);


###################################################
### code chunk number 10: Using_RCircos.Rnw:239-240
###################################################
RCircos.List.Parameters();


###################################################
### code chunk number 11: Using_RCircos.Rnw:252-255 (eval = FALSE)
###################################################
## out.file <- "RCircosDemoHumanGenome.pdf";
## pdf(file=out.file, height=8, width=8, compress=TRUE);
## RCircos.Set.Plot.Area();


###################################################
### code chunk number 12: Using_RCircos.Rnw:260-263 (eval = FALSE)
###################################################
## par(mai=c(0.25, 0.25, 0.25, 0.25));
## plot.new();
## plot.window(c(-2.5,2.5), c(-2.5, 2.5));


###################################################
### code chunk number 13: Using_RCircos.Rnw:273-274 (eval = FALSE)
###################################################
## RCircos.Chromosome.Ideogram.Plot();


###################################################
### code chunk number 14: Using_RCircos.Rnw:284-296 (eval = FALSE)
###################################################
## data(RCircos.Gene.Label.Data);
## 
## name.col <- 4;
## side <- "in";
## 
## track.num <- 1;
## RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,
## 			track.num, side);
## 
## track.num <- 2;
## RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, 
## 			name.col,track.num, side);


###################################################
### code chunk number 15: Using_RCircos.Rnw:306-312 (eval = FALSE)
###################################################
## data(RCircos.Heatmap.Data);
## data.col <- 6;
## track.num <- 5;	
## side <- "in";
## RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, 
## 			track.num, side);


###################################################
### code chunk number 16: Using_RCircos.Rnw:315-322 (eval = FALSE)
###################################################
## data(RCircos.Scatter.Data);
## data.col <- 5;
## track.num <- 6; 
## side <- "in";
## by.fold <- 1;
## RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, 
## 			track.num, side, by.fold);


###################################################
### code chunk number 17: Using_RCircos.Rnw:325-331 (eval = FALSE)
###################################################
## data(RCircos.Line.Data);
## data.col <- 5;
## track.num <- 7;
## side <- "in";
## RCircos.Line.Plot(RCircos.Line.Data, data.col, 
## 			track.num, side);


###################################################
### code chunk number 18: Using_RCircos.Rnw:334-340 (eval = FALSE)
###################################################
## data(RCircos.Histogram.Data);
## data.col <- 4;
## track.num <- 8; 
## side <- "in";
## RCircos.Histogram.Plot(RCircos.Histogram.Data, 
## 			data.col, track.num, side);


###################################################
### code chunk number 19: Using_RCircos.Rnw:343-347 (eval = FALSE)
###################################################
## data(RCircos.Tile.Data);
## track.num <- 9;
## side <- "in";
## RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);


###################################################
### code chunk number 20: Using_RCircos.Rnw:357-366 (eval = FALSE)
###################################################
## data(RCircos.Link.Data);
## track.num <- 11;
## RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
## 
## data(RCircos.Ribbon.Data);
## RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, 
## 	track.num=11, by.chromosome=FALSE, twist=FALSE);
## 
## dev.off();


###################################################
### code chunk number 21: Using_RCircos.Rnw:378-381 (eval = FALSE)
###################################################
## library(RCircos);
## demo("RCircos.Demo.Human");
## demo("RCircos.Demo.Mouse.And.Rat");


###################################################
### code chunk number 22: sessionInfo
###################################################
sessionInfo()


