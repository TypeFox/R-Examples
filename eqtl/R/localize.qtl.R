#####################################################################
#
# localize.qtl.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Karl W Broman wrote the find.flanking function part of R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: localize.qtl
#
######################################################################

######################################################################
#
# localize.qtl: Calculate QTL physical positions from esimated
#                genetic positions
#
######################################################################

`localize.qtl` <-
function( cross, peak, data.gmap, round )
{

	require(qtl)

	if ( length(class(cross)) < 2 || class(cross)[2] != "cross")
		stop("Input should have class \"cross\".")
	if ( ! all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")

	scanone <- get( attr(peak,'scanone',exact=TRUE) )

	if ( ! all(levels(scanone$chr) == names(cross$geno)) )
		stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")
	if ( ! all(scanone$pos == pseudo.map(cross)) )
		stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which has been performed on cross object")
	if ( class(data.gmap)[1] != "data.frame" || any(!(names(data.gmap) %in% c('Marker','chr','PP'))) )
		stop("data.gmap should have class \"data.frame\" with columns names: 'Marker','chr','PP'")
	if ( missing(data.gmap) )
		stop("Argument data.gmap unspecified")
	if ( ! all( as.vector(data.gmap$Marker) %in%  mnames.map(cross) ) )
		stop("Arguments cross or data.gmap misspecified: data.gmap should have the same marker names as cross")
	if ( any( c("peak.bp","inf.bp","sup.bp") %in% attr(peak,'features',exact=TRUE)) ){
		 cat("WARNING: The physical location were already computed. They will be replaced.\n")
		 peak <- drop.peakfeat(peak,c('peak.bp','inf.bp','sup.bp'))
	}
	if ( !missing(round) && round < 0 && !is.numeric(round) )
		stop("Argument round misspecified: round should be an integer >= 0")

	res <- list(un=NA)
	for ( i in seq(length(peak)) ){

		trait <- names(peak[i])
		resbytrait <- list(un=NA)

		for ( y in 1:length(peak[[i]])){

			if ( ! is.na(peak[[i]][y][1])){

				# INFORMATIVE SCREEN MESSAGE
				# print('QTL found')

				p <- peak[[i]][[y]]$mname.peak

				inf_bp=NA;sup_bp=NA;peak_bp=NA;

				for ( z in seq(length(p)) ){ 

					flank <- function(cross,chr,pos_genetic){
						mflank <- find.flanking(cross,chr,pos_genetic)
						ml_g <- find.markerpos(cross,mflank$left)$pos
						mr_g <- find.markerpos(cross,mflank$right)$pos
						bool <- data.gmap$Marker %in% paste(mflank$left)
						ml_p <- data.gmap$PP[bool]
						bool <- data.gmap$Marker %in% paste(mflank$right)
						mr_p <- data.gmap$PP[bool]
						return(data.frame(left=c(ml_g,ml_p),right=c(mr_g,mr_p)))
					}

					pic <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$peak.cM[z])
					sup <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$sup.cM[z])
					inf <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$inf.cM[z])

					if (pic$left[1] == pic$right[1] && pic$right[1] == 0)
						pic <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$peak.cM[z]+1)
					if (sup$left[1] == sup$right[1] && sup$right[1] == 0)
						sup <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$sup.cM[z]+1)
					if (inf$left[1] == inf$right[1] && inf$right[1] == 0)
						inf <- flank(cross,names(peak[[i]][y]),peak[[i]][[y]]$inf.cM[z]+1)

					if (pic$left[1] == pic$right[1]) ratio_pic <- pic$left[1]/pic$left[2]
						else ratio_pic = (pic$right[1]-pic$left[1])/(pic$right[2]-pic$left[2])
					if (sup$left[1] == sup$right[1])	ratio_sup <- sup$left[1]/sup$left[2]
						else ratio_sup = (sup$right[1]-sup$left[1])/(sup$right[2]-sup$left[2])
					if (inf$left[1] == inf$right[1])	ratio_inf <- inf$left[1]/inf$left[2]
						else ratio_inf = (inf$right[1]-inf$left[1])/(inf$right[2]-inf$left[2])

					inf_p <- (peak[[i]][[y]]$inf.cM[z]-inf$left[1])/ratio_inf + inf$left[2]
					sup_p <- (peak[[i]][[y]]$sup.cM[z]-sup$left[1])/ratio_sup + sup$left[2]
					pic_p <- (peak[[i]][[y]]$peak.cM[z]-pic$left[1])/ratio_pic + pic$left[2]

					if ( !missing(round) ){
						inf_p <- round(inf_p,round)
						sup_p <- round(sup_p,round)
						pic_p <- round(pic_p,round)
					}

					inf_bp <- c(inf_bp,inf_p)
					sup_bp <- c(sup_bp,sup_p)
					peak_bp <- c(peak_bp,pic_p)

				}
				pos_bp <- data.frame(peak.bp=peak_bp[-1],inf.bp=inf_bp[-1],sup.bp=sup_bp[-1])
				pos_bp <- cbind(peak[[i]][[y]],pos_bp)
				pos_bp <- list(pos_bp)
			} else {
				# INFORMATIVE SCREEN MESSAGE
				# cat(trait,": no QTL found on chromosome ",names(peak[[i]][y]),"\n")
				pos_bp <- NA
			}
			attributes(pos_bp)$names<-names(peak[[i]][y])
			resbytrait<-c(resbytrait,pos_bp)
		}
		resbytrait<-list(resbytrait[-1])
		attributes(resbytrait)$names <- trait
		res<-c(res,resbytrait)
	}

	res <- res[-1]
	attributes(res)$class <- c('peak','list')
	attributes(res)$features <-  c(attributes(peak)$features,'peak.bp','inf.bp','sup.bp')
	attributes(res)$scanone <- attributes(peak)$scanone
	attributes(res)$lod.th <- attributes(peak)$lod.th
	attributes(res)$si <- attributes(peak)$si
	attributes(res)$window <- attributes(peak)$window

	return(as.list(res))
}

