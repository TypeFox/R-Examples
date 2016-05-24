#' @export
#' @import ggplot2
plot_missings <- function(formule, data, show.ids=TRUE) {
	# parse formula
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");

	f0 <- all.vars(formule)
	grs <- long2matrix(formule, data, reduce=FALSE, skip3=FALSE)

	m1 <- reshape2::melt(grs)
	m1$value2 <- factor(ifelse(is.na(m1$value), 1, 0))
	colnames(m1)[1:2] <- c(f0[2], f0[3])
	m1[,f0[2]] <- factor(m1[,f0[2]], ordered=TRUE)
	m1[,f0[3]] <- factor(m1[,f0[3]])

	p2 <- ggplot(m1, aes_string(x=f0[3], y=f0[2])) + geom_point(aes_string(color="value2"), na.rm=TRUE) + facet_wrap(~L1, scales="free")
	p2 <- p2 + scale_colour_manual("Missing?", breaks=c(0,1), values=c("grey85", "black"), labels=c("no", "yes")) 

	p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust=1, size=7), axis.text.y = element_text(hjust=1, size=7), aspect.ratio=1) + ggtitle("Missing values")

	if (show.ids==FALSE) {
		p2 <- p2 + theme(axis.text.x = element_text(size=0), axis.text.y = element_text(size=0))
	}
	return(p2)
}




#' @export
#' @method plot RRbi
plot.RRbi <- function(x, ...) {
	plot.RRuni(x, ...)
}

#' @export
#' @method plot RRuni
plot.RRuni <- function(x, ..., measure=NA, geom="bar") {
	
	RRu <- x
	standardized <- type2 <- NULL
	if (is.na(measure)) {
		measure <- localOptions$style
	} else {
		measure <- match.arg(measure, c("behavior", "perception", "metaperception"))
	}
	
	if (measure=="behavior") {lab <- unilabels_b}
	if (measure=="perception") {lab <- unilabels_p}
	
	
	mode <- ifelse(length(RRu$univariate)==2,"bi","uni")
	
	if (mode=="uni") {
		df <- RRu$varComp
	} else {
		df <- rbind(cbind(RRu$univariate[[1]]$varComp, variable=1), cbind(RRu$univariate[[2]]$varComp, variable=2))
	}
	
	df$type2 <- lab
	df$standardized[is.na(df$standardized)] <- 0
	
	p1 <- ggplot(df[df$type2 %in% lab[1:3],], aes(x=factor(1), y=standardized, fill=as.character(type2))) + geom_bar(stat="identity", na.rm=TRUE) + scale_fill_discrete("Variance Component") + xlab("") + ylab("Standardized variance component")
	
	if (geom=="pie") p1 <- p1 + coord_polar(theta="y")
	
	if (mode=="bi") {p1 <- p1 + facet_wrap(~variable)}
	
	return(p1)
}




#' @export
#' @method plot RRmulti
plot.RRmulti <- function(x, ..., measure=NA, geom="scatter", conf.level=0.95, connect=FALSE) {
	RRm <- x
	
	# NULL out ggplot2 variables so that R CMD check is satisfied ...
	# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
	type <- jit.x <- type2 <- group.size <- estimate <- tcrit <- se <- group.id <- standardized <- NULL
	
	if (is.na(measure)) {
		measure <- localOptions$style
	} else {
		measure <- match.arg(measure, c("behavior", "perception", "metaperception"))
	}
	
	mode <- ifelse(length(RRm$univariate)==2,"bi","uni")
	
	if (mode=="uni") {
		df0 <- RRm$varComp.group
		grouplevel <- RRm$varComp
	} else {
		df0 <- rbind(RRm$univariate[[1]]$varComp.group, RRm$univariate[[2]]$varComp.group)
		grouplevel <- rbind(cbind(RRm$univariate[[1]]$varComp, variable=1), cbind(RRm$univariate[[2]]$varComp, variable=2))
	}
	
	
	if (measure=="behavior") {lab <- gsub(" ", "\n", unilabels_b, fixed=TRUE)}
	if (measure=="perception") {lab <- gsub(" ", "\n", unilabels_p, fixed=TRUE)}
	
	# set right order and factor labels
	df <- df0
	df$type <- factor(df$type, levels=unilabels_b, labels=lab)
	
	
	
	if (geom=="scatter") {
		#define deterministic jittering
		df$jit.x <- NA
		for (i in names(table(df$group.id))) df$jit.x[df$group.id==i] <- rnorm(1,0,0.1)
		
		grouplevel$type2 <- factor(grouplevel$type, levels=unilabels_b, labels=lab)
		grouplevel$tcrit <- abs(qt((1-conf.level)/2,length(RRm$groups)-1))

		
		p1 <- ggplot(df, aes_string(y="estimate"), na.rm=TRUE)

		p1 <- p1 + geom_point(aes(x=(as.numeric(type)+jit.x), size=group.size), alpha=0.6, color="grey60", na.rm=TRUE) + scale_size("Group size")
	
		p1 <- p1 + geom_point(aes(x=(as.numeric(type)+jit.x)), alpha=0.8, color="black", size=0.7, na.rm=TRUE)
		
		p1 <- p1 + scale_x_discrete()
		
		p1 <- p1 + geom_pointrange(data=grouplevel, aes(x=type2, y=estimate, ymin=estimate-tcrit*se, ymax=estimate+tcrit*se), color="darkgreen", size=1.1, na.rm=TRUE)
		
		
	
		# lines connecting the points (may be very cluttered)
		if (connect==TRUE) {
			p1 <- p1 +geom_line(data=df[as.numeric(df$type) <= 3,], aes(x=(as.numeric(type)+jit.x), y=estimate, group=group.id), color="grey40", alpha=0.7, na.rm=TRUE)
		}

	
		p1 <- p1 + theme(axis.text.x = element_text(angle = 0, vjust=1)) + xlab("") + ggtitle(paste("Multiple round robin groups:\nAbsolute (co-)variance estimates\nand",round(conf.level,2)*100,"%-CI (weighted for group size)"))
		
		if (mode=="bi") {
			p1 <- p1 + facet_wrap(~variable)
		}
		
		return(p1)
	}	
	
	
	
	if (geom=="bar") {

			
		# df3 <- df
		# x <- df3[df3$group.id==1,]
		# ddply(df3, .(group.id), function(x) {return(x$standardized[x$type=="actor\nvariance"])})
		# 
		# df3$group2 <- factor(
		# 	df3$group.id, levels=df3$group.id[df3$type=="actor\nvariance"][order(df3$standardized[df3$type=="actor\nvariance"], na.last=FALSE)],
		# 	ordered=TRUE)
		#df3$group2 <- factor(df3$group.id, levels=df3$group.id[order(df3$standardized[df3$type=="actor\nvariance"])])
		
		df3 <- na.omit(df)
		
		p2 <- ggplot(df3[as.numeric(df3$type)<=3,], aes(x=group.id, y=standardized, fill=as.character(gsub("\n", " ", type, fixed=TRUE)))) + geom_bar(aes(width=group.size), na.rm=TRUE)
		p2 <- p2 + scale_fill_discrete("Variance Component") + ylab("Standardized variances")
		
		if (mode=="bi") {p2 <- p2 + facet_wrap(~variable)}
		
		return(p2)
	}
	
}


