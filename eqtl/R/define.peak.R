#####################################################################
#
# define.peak.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Karl W Broman wrote the plot.scanone function part of R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/eqtl package
# Contains: define.peak
#
######################################################################

######################################################################
#
# define.peak : Define QTL with support interval and exclusionary window
#
######################################################################

`define.peak` <-
function( scanone, lodcolumn=1, chr, th=2.3, si=1.5, graph=FALSE, window.size=20, round, save.pict=FALSE, phe.name, ...)
{
	require(qtl)

	if( !all(attr(scanone,'class',exact=TRUE) %in% c('scanone','data.frame')) )
		stop("Input should have class \"scanone\"")

	if( !missing(chr) ){
		if( !is.vector(chr) ) stop("Argument chr misspecified: expecting a vector")
		if( !all(chr %in% levels(scanone$chr)) ) stop("Could not find all chr \"",chr,"\" in scanone$chr")
		if( length(chr) > length(levels(scanone$chr)) )
			stop("Argument chr misspecified: chr could not be longer than the number of chr. in scanone",length(levels(scanone$chr)) )
	} else 	chr <- levels(scanone$chr)

	if( !is.numeric(th) | !is.vector(th) | length(th) != 1)
		stop("Argument th misspecified: expecting a numeric vector of length 1")
	if( !is.numeric(si) | !is.vector(si) | length(si) != 1)
		stop("Argument si misspecified: expecting a numeric vector of length 1")
	if( !is.logical(graph) & !is.vector(graph) & length(graph) != 1)
		stop("Argument graph misspecified: expecting a logical vector of length 1")
	if( !is.vector(window.size) | !is.numeric(window.size) | length(window.size)>1 )
		stop("Argument window.size misspecified: expecting a numeric vector of length 1")
	if( window.size < 0)
		stop("Argument window.size misspecified: wrong value (centiMorgan)")

	if( missing(lodcolumn) ) stop("Argument lodcolumn unspecified")
	if( !is.vector(lodcolumn) & (!is.numeric(lodcolumn) | !is.character(lodcolumn)) )
		stop("Argument lodcolumn misspecified: expecting numeric vector or character vector")
	if( length(lodcolumn) > length(names(scanone))-2)
		stop("Argument lodcolumn misspecified: cannot refered to more lodcolumn than scanone contains. The vector is too big.")
	if( length(names(scanone)) == 3 ) lodcolumn <- 'all'

	if( length(lodcolumn) == 1 & lodcolumn[1] == 'all' ){
		trait <- names(scanone[3:length(names(scanone))])

		if ( ncol(scanone) == 3 & names(scanone)[3] == 'lod' ) {
			if ( missing(phe.name) )
				stop("Warning: the lodcolumn name is 'lod', argument phe.name is necessary")
			trait <- phe.name
			if ( toupper(phe.name) == 'ID' ) stop("Warning: ID is reserved name for indexing the phenotypes")
		}

		num <- c(1:(length(names(scanone))-2))
	} else {
		trait <- 'NA'
		num <- NA
		for( i in seq(length(lodcolumn))){
			if( is.numeric(lodcolumn[i]) ){
				if( lodcolumn[i] < 1 | lodcolumn[i] >= length(names(scanone))-2 )
					stop("lodcolumn values should be between 1 and the no. of lod columns")
				else 	trait <- c(trait,names(scanone[lodcolumn[i]+2]))
			}
			if( is.character(lodcolumn[i]) ){
				trait <-  c(trait,lodcolumn[i])
				if( ! toupper(lodcolumn[i]) %in% toupper(names(scanone)) )
					stop("Could not identify trait lod column \"",lodcolumn[i],"\"")
				lodcolumn[i] <- as.numeric( grep( toupper(lodcolumn[i]),toupper(names(scanone)) ) ) - 2
				
			}
			num <- c(num,as.numeric(lodcolumn[i]))
		}
		lodcolumn<-as.numeric(lodcolumn)
		trait <- trait[-1]; num<-num[-1]
	}

        res<-list(init=NA)
	resbytrait<-list(init=NA)

	#SCREEN MESSAGE: NUMBER OF TRAITS, NAME, COLUMN NUMBER IN SCANONE
	cat("no. of traits: ",length(trait),"\n",sep="")
	for( i in seq(trait)) cat("trait: ",trait[i],"\tlodcolumn: ",num[i],"\n",sep="")
	cat("define.peak in process...\n")

	if ( !missing(round) && !is.numeric(round) && round< 0 )
		stop("Argument round misspecified: round should be an integer >= 0")

	#FUNCTION #1: LOD CURVE EXTRACTION FOR A GIVEN TRAIT
	define.curve <- function( scanone, num, chr ){
		curve<-list(un=NA)
		for (i in chr){
			bool <- scanone$chr %in% i
			res<-list(list(mname=row.names(scanone)[bool],pos=scanone$pos[bool],lod=scanone[bool,num]))
			attributes(res)$names<-paste(i)
			curve<-c(curve,res)
		}
		curve<-curve[-1]
		invisible(curve)
	}

	curve<-define.curve(scanone,num+2,chr)

	for( z in 1:length(trait)){
		for( i in seq(length(chr)) ){

			c <- grep(paste("^",chr[i],"$",sep=''),names(curve))

			if(length(trait) > 1) lod<-curve[[c]]$lod[[z]]
			else lod<-curve[[c]]$lod

			pos <- curve[[c]]$pos

			if( max(lod) < th ){

				# OPTIONAL INFORMATIVE SCREEN MESSAGE
				# cat(	"chr ",chr[i],", trait no. ",z," '",trait[[z]],"' : no LOD score has been found beyond ",
				#	th," LOD, max= ",max(lod),"\n",sep=''	)

				resbychr <- list(NA)
				attributes(resbychr)$names <- chr[i]
				resbytrait<-c(resbytrait,resbychr)
				next
			}

			# OPTIONAL INFORMATIVE SCREEN MESSAGE
			# else cat("chr ",chr[i],", trait no.",z," '",trait[[z]],"' : lod max= ",max(lod),"\n",sep='')

			#FUNCTION #2: CONSERVE THE HIGHEST PEAK IN A EXCLUSIONARY WINDOW
			filter.peak <- function( pos=pos, lod=lod, dmin=20, th ){

				if( ! any(lod >= th) ) return()

				mx <- NA
				for( i in seq(length(pos)) ){
					if( i == 1 && lod[i] > lod[i+1] && lod[i] >= th) mx<-c(mx,i)
					if( i == length(lod) && lod[i] > lod[i-1]  && lod[i] >= th) mx <- c(mx,i)
					if( i != 1 && i!= length(lod) && lod[i] > lod[i-1] && lod[i] > lod[i+1]  && lod[i] >= th) mx <- c(mx,i)
				}

				if( length(mx) == 1 ){

					# OPTIONAL INFORMATIVE SCREEN MESSAGE
					# cat("WARNING: no lod peak has been found. Unable to analyze the lod curve.\n")

					return()
				} else mx <- mx[-1]

				if( length(mx) == 1 ){

					# OPTIONAL INFORMATIVE SCREEN MESSAGE
					# cat("WARNING: only one maximum lod peak has been found. No filter applied\n")

					return(mx)
				}

				# OPTIONAL INFORMATIVE SCREEN MESSAGE
				# cat("Exclusionary window applied\n")

				step <- NA
				for( i in seq(length(mx)) ){
					if(i == 1) next
					step <- c( step, pos[mx[i]]-pos[mx[i-1]] )
				}
				step <- step[-1]
				step <- min(step)

				del_peak <- NA

				for( p in seq(1,length(pos),by=step) ){
					if( all(!(pos[mx] >= p &  pos[mx] <= p+dmin)) ) next
					cple <-  mx[ pos[mx] >= p & pos[mx] <= p+dmin ]

					if( all(cple %in% del_peak) ) next
					cple <- cple[ ! cple %in% del_peak ]

					del_peak <- c( del_peak, cple[lod[cple] !=  max(lod[cple])] )
				}
				del_peak <- del_peak[-1]

				mx <- mx[! mx %in% del_peak ]
				return(mx)
			}

			max <- filter.peak( pos=pos, lod=lod, dmin=window.size, th=th )

			if( length(max) < 1 ) next

			#FONCTION #3: CONFIDENCE INTERVAL DEFINITION
			if (length(trait)>1) trait_curve <- data.frame(pos=curve[[c]]$pos,lod=curve[[c]]$lod[[z]])
			else trait_curve <- data.frame(pos=curve[[c]]$pos,lod=curve[[c]]$lod)

			limit.peak <- function( curve, peak, s, m=10 ){

				# OPTIONAL INFORMATIVE SCREEN MESSAGE
				# cat("SI computing\n")

				inf<-NA; sup<-NA; infq<-NA; supq<-NA;
				for( i in seq(length(curve$lod)) ){
					for( y in curve$pos[peak]){
	
						if( curve$pos[i] == y){
							b <- i
							a <- i
							if( i == 1 ){
								for (b in i:(length(curve$lod)-1) ){
									if( curve$lod[b]<=0){ q <- '*'; b<-b-1; break }
									if( curve$pos[b]-curve$pos[i] > m & m != 0){ q <- '->'; break }
									if( curve$lod[b]<curve$lod[i]-s){ q <- '+'; break }
									b <- b+1
								}
								sup<-c(sup,b)
								inf<-c(inf,i)
								infq<-c(infq,'|')
								supq<-c(supq,q)
							}

							if( i==length(curve$lod) ){
								for (a in i:2){
									if( curve$lod[a]<=0){ q <- '*'; a<-a+1; break }
									if( curve$pos[i] - curve$pos[a] > m & m != 0){ q <- '<-'; break }
									if( curve$lod[a] < curve$lod[i]-s){ q <- '+'; break }
									a <- a - 1
								}
								inf<-c(inf,a)
								sup<-c(sup,b)
								supq<-c(supq,'|')
								infq<-c(infq,q)
							}

							if( i != 1 && i != length(curve$lod) ){

								for (a in i:2){
									if( curve$lod[a]<=0) { q <- '*'; a<-a+1; break }
									if( curve$pos[i] - curve$pos[a] > m & m != 0){ q <- '<-'; break }
									if( curve$lod[a] < curve$lod[i]-s){ q <- '+'; break }
									q <- '|'
									a <- a-1
								}
								inf<-c(inf,a)
								infq<-c(infq,q)

								for (b in i:(length(curve$lod)-1)){
									if( curve$lod[b]<=0) { q <- '*'; b<-b-1;break}
									if( curve$pos[b] - curve$pos[i] > m & m != 0){ q <- '->'; break }
									if( curve$lod[b] < curve$lod[i]-s){ q <- '+'; break }
									q <- '|'
									b <- b+1
								}
								sup<-c(sup,b)
								supq<-c(supq,q)
							}
							break
						}
					}
				}
				inf <- inf[-1]
				sup <- sup[-1]
				supq <- supq[-1]
				infq <- infq[-1]
				qual <- paste(infq,supq,sep='')

				invisible(data.frame(inf=inf,max=peak,sup=sup,qual=qual))
			}
			peak<-limit.peak(trait_curve,max,si,...)

			# OPTIONAL INFORMATIVE SCREEN MESSAGE
			# cat(	"SI position for the trait ",trait[z]," on chr ",chr[i],": \n",sep="")

			if ( !missing(round)) trait_curve$lod[peak$max] <- round(as.numeric(as.vector(trait_curve$lod[peak$max])),round)

			# OPTIONAL INFORMATIVE SCREEN MESSAGE
			# print(curve[[c]]$pos[peak$inf])
			# print(curve[[c]]$pos[peak$max])
			# print(curve[[c]]$pos[peak$sup])

			resbychr <- list(data.frame(	'lod' = trait_curve$lod[peak$max],
							'mname.peak' = curve[[c]]$mname[peak$max],
							'peak.cM' = curve[[c]]$pos[peak$max],
							'mname.inf' = curve[[c]]$mname[peak$inf],
							'inf.cM' = curve[[c]]$pos[peak$inf],
							'mname.sup' = curve[[c]]$mname[peak$sup],
							'sup.cM' = curve[[c]]$pos[peak$sup],
							'si.quality' = peak$qual
					)	)


			attributes(resbychr)$names<- chr[i]
			resbytrait<-c(resbytrait,resbychr)

			if(save.pict) png( filename=paste(trait[z],i,z,".png",sep="_"), width=1280, height=1024 )
			if (graph | save.pict) {
				qtl:::plot.scanone(scanone,lodcolumn=num[z],chr=chr[i],show.marker.names=TRUE,lwd=1,main=paste('chr',chr[i]))
				abline(h=th,col="pink",lwd=1)
				abline(v=curve[[c]]$pos[peak$inf],col="blue")
				abline(v=curve[[c]]$pos[peak$sup],col="blue")
				abline(v=curve[[c]]$pos[peak$max],col="red",lty=2)
			}
			if(save.pict)dev.off()
		}

		resbytrait <- list(resbytrait[-1])
		if (length(trait)>1) attributes(resbytrait)$names <- paste(names(curve[[c]]$lod[z]))
		else attributes(resbytrait)$names <- paste(trait)
		res <- c(res,resbytrait)
	}

	res <- res[-1]

	attributes(res)$class <- c('peak','list')
	attributes(res)$features <- c('lod','mname.peak','peak.cM','mname.inf','inf.cM','mname.sup','sup.cM','si.quality')
	attributes(res)$scanone <- deparse(substitute(scanone))
	attributes(res)$lod.th <- th
	attributes(res)$si <- si
	attributes(res)$window <- window.size

	# RETURN
	return(res)
}

