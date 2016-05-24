###
topoplot<-function(erpobj, startmsec=-200, endmsec=1200, win.ini, win.end, exclude = NULL,
elec.coord=NULL, projection="orthographic", palette.col="jet", palette.steps=10, return.coord = FALSE,
zlim=NULL, interpolation = "cubicspline", extrap = TRUE, interp.points = 500, return.notfound=FALSE, mask = TRUE,  contour=TRUE, x.rev=FALSE,
draw.elec.pos=TRUE,  draw.nose=FALSE, draw.elec.lab=TRUE, elec.lab.adj=c(0.5, NA), head.col="black", head.lwd=1, ...)

{
	requireNamespace("akima")
	#  initial checks
	#projection checks
	if (!projection%in%c("orthographic", "equalarea")){
		stop("Available projectios are: orthographic, equalarea
		The projection specified is: ", projection, call.=F)
	}
	
	# palette col checks
	available.palettes=c("jet", "heat")
	if ((!palette.col[1]%in%available.palettes)&length(palette.col)==1){
		
		available.palettes=paste("\"",available.palettes,"\"", collapse=", ", sep="")
		stop("To create a palette for the plot, at least two colors must be specified.
		Type colors() to see a list of available colors.
		Two default palettes are available:", available.palettes, ".", call.=F)
		
	}
	
	# check zlims (SEE AFTER) to avoid code duplication I moved this section after iterpolation.
	# In the case of cubic interpolation, interpolated data may be higher or lower than original data.	
	# for this reason I check the limits afterwards.
	
	
	
	
	##################
	# RETRIEVE ELECTRODE POSITIONS
	########################
	if (is.null(elec.coord)){
		elec.coord=structure(list(el.name = structure(c(81L, 117L, 82L, 68L, 70L, 
69L, 10L, 11L, 3L, 1L, 2L, 4L, 7L, 5L, 6L, 8L, 45L, 43L, 41L, 
39L, 77L, 40L, 42L, 44L, 46L, 66L, 64L, 62L, 60L, 61L, 63L, 65L, 
67L, 74L, 72L, 51L, 49L, 47L, 59L, 48L, 50L, 52L, 73L, 71L, 75L, 
57L, 55L, 53L, 54L, 56L, 58L, 76L, 121L, 119L, 17L, 15L, 13L, 
14L, 16L, 18L, 120L, 118L, 128L, 23L, 21L, 19L, 20L, 22L, 24L, 
129L, 125L, 123L, 30L, 28L, 26L, 38L, 27L, 29L, 31L, 124L, 122L, 
126L, 36L, 34L, 32L, 33L, 35L, 37L, 127L, 99L, 97L, 95L, 93L, 
90L, 115L, 92L, 94L, 96L, 98L, 91L, 113L, 111L, 112L, 114L, 105L, 
103L, 101L, 110L, 102L, 104L, 100L, 106L, 108L, 109L, 85L, 89L, 
86L, 107L, 87L, 88L, 78L, 80L, 79L, 12L, 83L, 84L, 9L, 25L, 116L, 
130L), .Label = c("AF3", "AF4", "AF7", "AF8", "AFF1h", "AFF2h", 
"AFF5h", "AFF6h", "AFp10h", "AFp3", "AFp4", "AFp9h", "C1", "C2", 
"C3", "C4", "C5", "C6", "CCP1h", "CCP2h", "CCP3h", "CCP4h", "CCP5h", 
"CCP6h", "Centroid", "CP1", "CP2", "CP3", "CP4", "CP5", "CP6", 
"CPP1h", "CPP2h", "CPP3h", "CPP4h", "CPP5h", "CPP6h", "CPz", 
"F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "FC1", "FC2", 
"FC3", "FC4", "FC5", "FC6", "FCC1h", "FCC2h", "FCC3h", "FCC4h", 
"FCC5h", "FCC6h", "FCz", "FFC1h", "FFC2h", "FFC3h", "FFC4h", 
"FFC5h", "FFC6h", "FFT7h", "FFT8h", "Fp1", "Fp2", "Fpz", "FT10", 
"FT7", "FT8", "FT9", "FTT7h", "FTT8h", "Fz", "I1", "I2", "Iz", 
"Left", "Nasion", "NFp1h", "NFp2h", "O1", "O2", "OI1h", "OI2h", 
"Oz", "P1", "P10", "P2", "P3", "P4", "P5", "P6", "P7", "P8", 
"P9", "PO10", "PO3", "PO4", "PO7", "PO8", "PO9", "POO1", "POO10h", 
"POO2", "POO9h", "POz", "PPO1h", "PPO2h", "PPO3h", "PPO4h", "Pz", 
"Ref", "Right", "T10", "T7", "T8", "T9", "TP10", "TP7", "TP8", 
"TP9", "TPP7h", "TPP8h", "TTP7h", "TTP8h", "Cz"), class = "factor"), 
    d = c(108, 114, 110, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 
    69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 99, 120, 69), y = c(0, 
    0, 0.923, 0.9511, 1, 0.9511, 0.9565, 0.9565, 0.809, 0.892, 
    0.8919, 0.809, 0.7777, 0.8289, 0.8289, 0.7777, 0.5878, 0.6343, 
    0.6726, 0.6979, 0.7067, 0.6979, 0.6726, 0.6343, 0.5878, 0.4741, 
    0.5107, 0.5384, 0.5533, 0.5533, 0.5384, 0.5107, 0.4741, 0.2852, 
    0.309, 0.3373, 0.3612, 0.377, 0.3826, 0.377, 0.3612, 0.3373, 
    0.309, 0.2852, 0.164, 0.1779, 0.1887, 0.1944, 0.1944, 0.1887, 
    0.1779, 0.164, -1e-04, 0, 1e-04, 1e-04, 2e-04, 1e-04, 1e-04, 
    1e-04, 0, 0, -0.1639, -0.1778, -0.1883, -0.194, -0.194, -0.1884, 
    -0.1778, -0.1639, -0.2852, -0.309, -0.3372, -0.3609, -0.3767, 
    -0.3822, -0.3767, -0.3608, -0.3372, -0.309, -0.2853, -0.474, 
    -0.5106, -0.5384, -0.5532, -0.5532, -0.5384, -0.5106, -0.474, 
    -0.5429, -0.5878, -0.6342, -0.6724, -0.6975, -0.7063, -0.6975, 
    -0.6724, -0.6342, -0.5878, -0.5429, -0.8108, -0.8284, -0.8284, 
    -0.8108, -0.7467, -0.809, -0.8918, -0.923, -0.8918, -0.809, 
    -0.7467, -0.9739, -0.9739, -0.873, -0.9511, -1, -0.9511, 
    -0.873, -0.9679, -0.9679, -0.8785, -0.923, -0.8785, 0.8732, 
    0.9601, 0.9601, 0.8732, 0, 0, 0), x = c(-0.9237, 0.9237, 
    0, -0.309, 0, 0.3091, -0.2508, 0.2508, -0.5878, -0.3554, 
    0.3553, 0.5878, -0.5417, -0.122, 0.122, 0.5417, -0.809, -0.721, 
    -0.5399, -0.2888, 0, 0.2888, 0.5399, 0.721, 0.809, -0.8642, 
    -0.7218, -0.4782, -0.1672, 0.1672, 0.4782, 0.7218, 0.8642, 
    -0.8777, -0.9511, -0.8709, -0.6638, -0.3581, 0, 0.3581, 0.6638, 
    0.8709, 0.9511, 0.8777, -0.9669, -0.8184, -0.5466, -0.1919, 
    0.1919, 0.5466, 0.8184, 0.9669, -0.9237, -1, -0.9237, -0.7066, 
    -0.3824, 0.3824, 0.7066, 0.9237, 1, 0.9237, -0.9669, -0.8185, 
    -0.5465, -0.1918, 0.1918, 0.5465, 0.8185, 0.9669, -0.8777, 
    -0.9511, -0.8712, -0.6635, -0.358, 0, 0.358, 0.6635, 0.8712, 
    0.9511, 0.8777, -0.8639, -0.722, -0.4786, -0.1673, 0.1673, 
    0.4786, 0.722, 0.8638, -0.7472, -0.809, -0.7211, -0.5401, 
    -0.2889, 0, 0.2889, 0.5401, 0.7211, 0.809, 0.7472, -0.352, 
    -0.122, 0.122, 0.3519, -0.5425, -0.5878, -0.3549, 0, 0.3549, 
    0.5878, 0.5425, -0.1285, 0.1286, -0.4448, -0.309, 0, 0.309, 
    0.4448, -0.1533, 0.1533, -0.2854, 0, 0.2854, -0.4449, -0.1564, 
    0.1564, 0.4449, 0, 0.9237, 0), z = c(-0.3826, -0.3826, -0.3824, 
    1e-04, 1e-04, 0, 0.1438, 0.1437, 0, 0.2782, 0.2783, 0, 0.3163, 
    0.5452, 0.5452, 0.3163, 0, 0.2764, 0.5043, 0.6542, 0.7067, 
    0.6542, 0.5043, 0.2764, 0, 0.1647, 0.4651, 0.6925, 0.8148, 
    0.8148, 0.6925, 0.4651, 0.1647, -0.3826, 0, 0.3549, 0.6545, 
    0.8532, 0.9233, 0.8532, 0.6545, 0.3549, 0, -0.3826, 0.1915, 
    0.5448, 0.8154, 0.9615, 0.9615, 0.8154, 0.5448, 0.1915, -0.3826, 
    0, 0.3826, 0.7066, 0.9231, 0.9231, 0.7066, 0.3826, 0, -0.3826, 
    0.1915, 0.5449, 0.8153, 0.9611, 0.9611, 0.8153, 0.5449, 0.1915, 
    -0.3826, -1e-04, 0.3552, 0.6543, 0.8534, 0.9231, 0.8534, 
    0.6543, 0.3552, -1e-04, -0.3826, 0.1646, 0.4653, 0.6933, 
    0.8155, 0.8155, 0.6933, 0.4653, 0.1646, -0.3826, -1e-04, 
    0.2764, 0.5045, 0.6545, 0.7065, 0.6545, 0.5045, 0.2764, -1e-04, 
    -0.3826, 0.4659, 0.5453, 0.5453, 0.4659, -0.3825, 0, 0.2776, 
    0.3824, 0.2776, 0, -0.3825, 0.1822, 0.1822, -0.195, 0, 0, 
    0, -0.195, -0.195, -0.195, -0.3824, -0.3823, -0.3824, -0.1949, 
    -0.1949, -0.1949, -0.1949, 0, -0.5773, 1)), .Names = c("el.name", 
"d", "y", "x", "z"), row.names = c(NA, 130L), class = "data.frame")
			}
	
	# additional check if external electrode coordinate is supplied
	necessary.el.coord.cols=c("el.name", "y", "x")
	if (!all(necessary.el.coord.cols%in%names(elec.coord))) {
		colsnotfound=necessary.el.coord.cols[!necessary.el.coord.cols%in%names(elec.coord)]
		colsnotfound=paste("\"",colsnotfound,"\"", collapse=", ", sep="")
		stop("The electrode coordinate object should contain at least these columns: 
		1) \"el.name\": a column containing electrode names.
		2) \"y\": a column containing the y of electrode coordinates.
		3) \"x\": a column containing the x of electrode coordinates.
		The following column(s) have not been found in the object supplied: ", colsnotfound, call.=FALSE)
	}
	
	# show coordinates of electrodes already contained in the function and the function breaks
		if (return.coord == TRUE){
		return(elec.coord)
	} else {
	
	
	### EXCLUDE ELECTRODES
	
	erpobj=erpobj[,!names(erpobj)%in%exclude]
	
	# create a vector with the electrode to be included in the future steps
	curr.el=names(erpobj[,!names(erpobj)%in%exclude])
	
	up.curr.el=toupper(curr.el)
	up.elec.coord.names=toupper(elec.coord$el.name)	
	
	found=names(erpobj)[up.curr.el%in%up.elec.coord.names]
	notfound=names(erpobj)[!up.curr.el%in%up.elec.coord.names]
	# notice that I use the indices with upper case, but I retrive the original names
	
	if (length(notfound)>0&return.notfound==FALSE){
		notfound.list=paste(notfound, "\n", sep="")
		stop("the following electrodes have not been found\n", notfound.list, call.=F)
		}
	
	if (length(notfound)>0&return.notfound==TRUE&is.null(exclude)){
		return(notfound)
	}
	
	}
	
	# CRUCIAL!! I reorder both the names of elec.coord and the names of erpobj for matching
	# i have to do so, because I'm not taking into account the letter case.
	curr.elec.coord=elec.coord[up.elec.coord.names%in%up.curr.el,]
	curr.elec.coord=curr.elec.coord[order(as.character(curr.elec.coord$el.name)),]
	
	erpobj=erpobj[,order(names(erpobj))]
	
	# CHANGE X DIRECTIONS IF REQUIRED
	if (x.rev==TRUE){
		curr.elec.coord$x=-curr.elec.coord$x
	}
		
	### IMPORTANT restore the correct order of x and y.
	
		
	##################
	# PROJECT ELECTRODE COORDINATES FROM 3D TO 2D
	####################################
	
	if(projection=="orthographic"){
		x=curr.elec.coord$x
		y=curr.elec.coord$y
	}
	
	if (projection=="equalarea"){
	x=curr.elec.coord$x *sqrt(1/(1 + curr.elec.coord$z))
	y=curr.elec.coord$y *sqrt(1/(1 + curr.elec.coord$z))

	}
		
	
	#build palette
	if (palette.col[1]=="jet"){
		 mypalette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	
	} else {
		if (palette.col[1]=="heat"){
			 mypalette <- colorRampPalette(heat.colors(palette.steps))
	
		} else {
			mypalette <-colorRampPalette(palette.col)
		}
	}
	
	
	
	##################
	# RETRIEVE AMPLITUDE. If necessary compute means across timepoints. 
	####################################

	
	
	if (win.ini == win.end){
		ampl = as.numeric(erpobj[round(msectopoints(win.ini, dim(erpobj)[1], startmsec, endmsec)), ])
	} else {
		ampl = colMeans(erpobj[round(msectopoints(win.ini, dim(erpobj)[1], startmsec, endmsec)):
		round(msectopoints(win.end, dim(erpobj)[1], startmsec, endmsec)), ])
	}
	
	###
	# SET GRAPHIC MARGINS. They are not currently manipolable by user
	#####
	# I do this in this step because I interpolate also according to margins
	xlim=c(-1.3, 1.3)
	ylim=c(-1.3, 1.3)
	


	
	###############
	# INTERPOLATE DATA
	########################
	if (interpolation=="cubicspline"){
		interp.linear=F
	}
	if (interpolation=="linear"){
		interp.linear=T
		extrap=FALSE #notice that extrapolation is not possible with linear interp
	}
	
	interp.data=interp(x,y, ampl, xo=seq(xlim[1], xlim[2], length = interp.points), yo=seq(ylim[1], ylim[2], length = interp.points), linear=interp.linear, extrap= extrap)


		
	
	
	################
	# TOPOPLOT
	################
	
		
	
	# reduce margin to avoid in plotting that interpolation goes too much near the margins.
	interp.xlim.up=which(interp.data$x>(xlim[2]-0.1)) #o
	interp.xlim.inf=which(interp.data$x<(xlim[1]+0.1)) #

	interp.ylim.up=which(interp.data$y>(ylim[2]-0.1)) #
	interp.ylim.inf=which(interp.data$y<(ylim[1]+0.1))#
	
	# data trimming with the new limits imposed
	interp.data$z[c(interp.xlim.up, interp.xlim.inf),]=NA
	interp.data$z[,c(interp.ylim.up, interp.ylim.inf)]=NA
	
	
	#### CHECK PROBLEMS WITH ZLIM
	
	range.used=round(range(interp.data$z, na.rm=T),2)
	
	### create zlims if not specified
	# this lims are created to be symmetric (centered on 0) .
	if (is.null(zlim)){
		newlim=round(max(abs(range.used))) # of the absolute value of observed data I take the max. Then, I round with 0 digits.
		zlim=c(-newlim, newlim)
	}
	
	
	if ((min(zlim)>min(range.used))|(max(zlim)<max(zlim))){
		cat("WARNING: your data (after interpolation) are out of range as compared to the zlims specified.\n",
	"Your data range is:", paste(range.used, collapse= ", "), "\n",
		"Your zlims are:", paste(zlim, collapse=", "), "\n")
	} # 

	
	
		
	
	# plot
	# notice that xlim and ylim are reversed. But this doesn't matter.
	
	image(interp.data, col=mypalette(palette.steps), xlim=c(xlim[1], xlim[2]), ylim=c(ylim[1], ylim[2]), zlim=c(zlim[1], zlim[2]), frame.plot=FALSE, axes=FALSE, ...)

if (contour == TRUE){ 	
	
	cont_levels=seq(zlim[1], zlim[2], dist(zlim)/palette.steps)
	contour(interp.data, ylim=c(xlim[1], xlim[2]), xlim=c(ylim[1], ylim[2]), zlim=c(zlim[1], zlim[2]), frame.plot=FALSE, axes=FALSE, add=TRUE, drawlabels=FALSE, levels=cont_levels)
	
	}
	
	
		
	
	#####
	# ADD MASK
	#################
	
if (mask == TRUE){ 	
	#add mask (optional)
	
	plotcircle <- function (x, y, r, ...) 
	{
		#x^2+y^2 = r^2
	
		angle = seq(0, 2*pi, length = 200)
		xc = x + (r * cos(angle))
		yc = y + (r * sin(angle))
		polygon(xc, yc, ...)
		return(list(x = xc, y = yc))	
	}

	circ.coord=plotcircle(0,0,1, border = head.col, lwd = head.lwd)
	
	circle.points=length(circ.coord$x)

	#nota il grafico è costruito da dx verso sx, perchè i punti in circ.coord sono ordinati da 1 a -1 come x.
	pol.x=c(xlim[2], circ.coord$x[1:120], xlim[1], xlim[1], xlim[2], xlim[2])
	pol.y=c(0.4, -circ.coord$y[1:120], 0.4, ylim[1], ylim[1], 0.4)
	# nota che metto -0.2, ma sarebbe 0. non metto esattamente 0, perché altrimenti le due metà della maschera si sfiorano e rimane una piccola riga colorata.

	polygon(pol.x, pol.y, col="white", lty="blank")
	polygon(pol.x, -pol.y, col="white", lty="blank")
	}
	
	
		### ADD ELECTRODES AND TEXT
	#add electrodes locations and labels (optional)
	# TO BE ADJUSTED!!! devi permettere un controllo più fine di dimensioni elettrodi, caratteri, etc.
	
	if (draw.elec.pos==TRUE){
		points(x,y, pch=19, col="black")
	}
	
	if (draw.elec.lab==TRUE){
	text(x,y+0.1, labels=names(erpobj), adj=elec.lab.adj)
	}

	
	#### DRAW NOSE #####
	if (draw.nose==TRUE){
	
	Nose.left=structure(list(x = c(-0.165709488227992, -0.170777055719452, 
	-0.160641920736532, -0.120101380804849, -0.0542230034158648, 
	0), y = c(0.991010207442792, 1.00031714102429, 1.02885836681128, 
	1.06310783775568, 1.13731502480188, 1.16)), .Names = c("x", "y"
	), row.names = 2:7, class = "data.frame")

	Nose.right=structure(list(x = c(0.165709488227992, 0.170777055719452, 0.160641920736532, 
	0.120101380804849, 0.0542230034158648, 0), y = c(0.991010207442792, 
	1.00031714102429, 1.02885836681128, 1.06310783775568, 1.13731502480188, 
	1.16)), .Names = c("x", "y"), row.names = 2:7, class = "data.frame")
	
	lines(Nose.left, col=head.col, lwd=head.lwd)
	lines(Nose.right, col=head.col, lwd=head.lwd)
	}


	


	
	##############
	### RETURN	SOME PARAMETERS
	#####################
	invisible(list(zlim=zlim, palette=mypalette(palette.steps)))
	}
