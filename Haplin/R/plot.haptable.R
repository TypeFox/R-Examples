plot.haptable <- function(x, separate.plots = F, filename, filetype = "png", use.dd, verbose = T, ...)
{
##
## PLOT A HAPTABLE
## MERK: DENNE HAR MYE FELLES MED plot.haplin, BURDE KANSKJE VAERT SAMKJOERTE
##
#
##
## MERK! HAR IKKE IMPLEMENTERT use.single, SOM BRUKES TIL AA PLOTTE BOYS ONLY 
## I X-CHROM. DEN ER HELLER IKKE IMPLEMENTERT I f.plot.effects. SE IMPLEMENTERING I 
## plot.haplin. ?? ER DEN BRUKBAR TIL NOE SOM HELST?


##
## HAR IKKE INNE ARGUMENTET "reference", SOM LIGGER I plot.haplin

.coef <- coef.haptable(x)
.markernames <- na.omit(x$marker)
.info <- attr(.coef, "info")
.haplos <- attr(.coef, "haplos")
.n.sel.haplo <- length(.haplos)
.maternal <- .info$model$maternal
.poo <- .info$model$poo
.comb.sex <- .info$model$comb.sex
.ref.cat <- .info$haplos$ref.cat
.reference.method <- .info$haplos$reference.method
#
if(missing(use.dd)){
	.use.dd <- 1:.n.sel.haplo
}else{
	.use.dd <- use.dd
}
## OVERRIDE
if(identical(.comb.sex, "males")) .use.dd <- F


#use.single!


.use.dd.mat <- .use.dd
#
.use.single.mat <- 1:.n.sel.haplo
#
## PRODUCE JPEG OR PNG, IF REQUESTED
.prod.file <- !missing(filename)

#
##
#
## PARAMETERS FOR CHILD PLOT
.params <- list(coeff = .coef, ref.cat = .ref.cat, reference.method = .reference.method, haplos = .haplos, markernames = .markernames, maternal = .maternal, poo = .poo, use.dd = .use.dd, verbose = verbose, ...)
#
## CHANGE SOME PARAMETERS FOR MATENAL PLOT
.params.mat <- .params
.params.mat$use.dd <- .use.dd.mat
.params.mat$use.single <- .use.single.mat
#
##
if(!.prod.file){
	## RETAIN OLD PARAMETERS
	.oldpar <- par(no.readonly = T)
	on.exit(par(.oldpar))
	#
	if(!.maternal){
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
	}
	if(.maternal & !separate.plots){
		par(mfrow = c(2,1), oma = c(2,0,0,0))
		.par <- par(no.readonly = T)
		.params$type <- 3
		invisible(do.call("f.plot.effects", .params))
		par(mar = .par$mar)
		.params.mat$type <- 4
		invisible(do.call("f.plot.effects", .params.mat))
	}
	if(.maternal & separate.plots){
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		.params.mat$type <- 2
		invisible(do.call("f.plot.effects", .params.mat))
	}
}# END !.prod.file


if(.prod.file){
	.jpeg.size <- c(440, 460)
	if(.n.sel.haplo > 4) .jpeg.size <- .jpeg.size + (.n.sel.haplo - 4) * c(30,0)
	#
	if(!.maternal){
		if(filetype == "png"){
			png(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		dev.off()
	}
	if(.maternal & !separate.plots){
		if(filetype == "png"){
			png(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = filename, width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		par(mfrow = c(2,1), oma = c(2,0,0,0))
		.par <- par(no.readonly = T)
		.params$type <- 3
		invisible(do.call("f.plot.effects", .params))
		par(mar = .par$mar)
		.params.mat$type <- 4
		invisible(do.call("f.plot.effects", .params.mat))
		dev.off()
	}
	if(.maternal & separate.plots){
		if(filetype == "png"){
			png(filename = paste("1", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = paste("1", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params$type <- 1
		invisible(do.call("f.plot.effects", .params))
		dev.off()
		if(filetype == "png"){
			png(filename = paste("2", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9)
		}
		if(filetype == "jpeg"){
			jpeg(filename = paste("2", filename, sep = ""), width = .jpeg.size[1], height = .jpeg.size[2], pointsize = 9, quality = 100)
		}
		.params.mat$type <- 2
		invisible(do.call("f.plot.effects", .params.mat))
		dev.off()	
	}
}# END .prod.file


}

