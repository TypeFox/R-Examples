plot.slidWin <- 
function(x, outliers = FALSE, ...){
	slidWin <- x
	distM <- function(){
		plot(slidWin$pos_out, slidWin$dist_mean_out, xlab = "Window position", ylab = "Distance", main="Mean K2P distance of each window")
		ylimzero <- range(slidWin$zero_out, slidWin$dat_zero_out, na.rm=TRUE)
		plot(slidWin$pos_out, slidWin$zero_out, ylim = ylimzero, xlab = "Window position", ylab = "Proportion", main="Proportion of zero cells in K2P distance matrix")
			abline(h=slidWin$dat_zero_out)
		#thres_main <- paste("Number of cells above ", slidWin$thresA, " (black) and below ", slidWin$thresB, " (blue)", sep="")
		#ylimthres <- range(slidWin$thres_above_out, slidWin$thres_below_out, na.rm=TRUE)
		#plot(slidWin$pos_out, slidWin$thres_above_out, ylim = ylimthres, xlab = "Position", ylab = "Number of cells", main = thres_main)
			#points(slidWin$pos_out, slidWin$thres_below_out, col = "blue")
		plot(slidWin$pos_out, slidWin$nd_out, xlab = "Window position", ylab = "Number", main = "Sum of diagnostic nucleotides", type="l")
		plot(slidWin$pos_out, slidWin$noncon_out, ylim=c(0, max(slidWin$noncon_out)), xlab = "Window position", ylab = "Proportion", main = "Proportion of zero non-conspecific K2P distances")
	}

	distT <- function(){
		ylimcomp <- range(slidWin$comp_out, slidWin$comp_depth_out, na.rm=TRUE)
		plot(slidWin$pos_tr_out, slidWin$comp_out, ylim = ylimcomp, xlab = "Window position", ylab = "Proportion", main="Congruence of NJ trees")
		#	points(slidWin$pos_tr_out, slidWin$comp_depth_out, col="blue")
		plot(slidWin$pos_tr_out, slidWin$win_mono_out, xlab = "Window position", ylab = "Proportion", main = "Proportion of species that are monophyletic")
	}

	if(slidWin$boxplot_out){
		layout(1)
		if(slidWin$method == "overall"){
			plot(1, type = "n", xlim = c(1,length(slidWin$bp_out)), ylim = slidWin$bp_range_out, xaxt = "n", main="Boxplot of distance matrix for each window", xlab="Window position", ylab="Distance")
			axis(1, at = 1:length(slidWin$bp_out), labels = as.character(slidWin$pos_out))
			for(i in 1:length(slidWin$bp_out)) bxp(slidWin$bp_out[[i]], at=i, add=TRUE, axes=FALSE, whisklty=1, whiskcol="red", staplelty=0, outpch=NA)
		}#x11()
		if(slidWin$method %in% c("interAll", "nonCon")){
			labs <- c("interAll", "nonCon")
			desc <- c("all inter-specific", "closest non-conspecific")
			plotTitle <- paste("Boxplots of", desc[which(slidWin$method==labs)], "distances (orange whiskers) and \nintra-specific distances (blue whiskers) in each window", sep=" ")
			plot(1, type = "n", xlim = c(1,length(slidWin$bp_InterSpp_out)), ylim = c(0, median(slidWin$bp_range_out)), xaxt = "n", main=plotTitle, xlab="Window position", ylab="Distance")
			axis(1, at = 1:length(slidWin$bp_InterSpp_out), labels = as.character(slidWin$pos_out))
			if(outliers){
				for(i in 1:length(slidWin$bp_InterSpp_out)) bxp(slidWin$bp_InterSpp_out[[i]], at=i, add=TRUE, axes=FALSE, whisklty=1, whiskcol="orange", staplelty=0, outpch=20, outcol="orange")
				for(i in 1:length(slidWin$bp_IntraSpp_out)) bxp(slidWin$bp_IntraSpp_out[[i]], at=i, add=TRUE, axes=FALSE, whisklty=1, whiskcol="blue", staplelty=0, outpch=20, outcol="blue", outcex=0.75)
			}
			else{
				for(i in 1:length(slidWin$bp_InterSpp_out)) bxp(slidWin$bp_InterSpp_out[[i]], at=i, add=TRUE, axes=FALSE, whisklty=1, whiskcol="orange", staplelty=0, outpch=NA)
				for(i in 1:length(slidWin$bp_IntraSpp_out)) bxp(slidWin$bp_IntraSpp_out[[i]], at=i, add=TRUE, axes=FALSE, whisklty=1, whiskcol="blue", staplelty=0, outpch=NA)
			}
		}
	}

	if(slidWin$distMeasures && slidWin$treeMeasures){
		layout(matrix(1:6, ncol=2))
		distM()
		distT()
	} else {
			if(slidWin$distMeasures) {
				layout(1:4)
				distM() }
			if(slidWin$treeMeasures){
				layout(1:2)
				distT()}				
		}
layout(1)
}
