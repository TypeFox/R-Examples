#---------------------------- GUI function ------------------------------------#
# calcAssess (calculate assessment):
# Purpose:     Calculates all the data needed for assessment and requests that it
#                be plotted appropriately
# Parameters:  Assessment parameters as50, as95, am50, am95, winf, linf, b, t0,
#                k, h1, h2, sigma1, sigma2, tau1, tau2, A, T, M, R, q, gamma1,
#                chk, unitType, plotType, percent, maxB, powr
# GUI inputs:  Assessment parameters as50, as95, am50, am95, winf, linf, b, t0,
#                k, h1, h2, sigma1, sigma2, tau1, tau2, A, T, M, R, q, gamma1,
#                chk, unitType, plotType, percent, maxB, powr
# Returns:     A list containing the calculated data betaa, ma, wa, Ft, Rt, RBt,
#                DetN, Nbar, N, Pt, DetPt, PBt, DetPBt, Ct, DetCt, CBt, DetCBt,
#                Ibar, StoI, IBbar, StoBI, uat, Detuat, uBat, DetuBat, Pat, 
#                PBat, DetNB, NB, St, SBt, Tt, TBt
# GUI outputs: None
# Notes:       Usually this will be called with no parameters, then it will grab
#                the information from the GUI
#              The parameter chk is a logical saying whether or not we should
#                do error checking
#-------------------------------------------------------------------------------
local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

calcAssess <- function(as50, as95, am50, am95, winf, linf, b, t0, k, h1, h2, 
    sigma1, sigma2, tau1, tau2, A, T, M, R, q, gamma1, unitType, plotType, 
    percent, maxB, powr, chk=FALSE, act=NULL)
{     
	if(is.null(act)) act <- getWinAct()[1];
  if(!is.null(act) && act == "interact") {
    print("Select a point on the R Graphics window to display its x and y coordinates")
    pt <- locator(1)
    print(paste("x is  ", round(pt$x, digits=4)))
    print(paste("y is  ", round(pt$y, digits=4)))
    return(invisible())
  }
	if(missing(as50) && missing(as95) && missing(winf)) getWinVal(scope="L")
	if((chk) && (!validParam()))  return(invisible())

# Calculate random numbers if there are none, or if the user wants to recalculate them
	if(!is.null(act) && (act != "disp") || (!exists(".PBSpopsim",where=.PBSmodEnv))) {
		.PBSpopsim <- list(colPre="green",colHigh="red",colLow="deepskyblue")
		nAT <- length((2-A):T); RtOffset <- A - 1;
		delta <- .getRnorm(nA=A, nY=nAT);
		epsilon <- .getRnorm(nA=A, nY=T)
		epsilon0 <- .getRnorm(nA=1, nY=nAT);
		.PBSpopsim <- c(.PBSpopsim,list(RtOffset=RtOffset,
			delta=delta,epsilon=epsilon,epsilon0=epsilon0));

		betaa <- getSelectivity(as50, as95, A)                 # selectivity
		ma    <- getMaturity(am50, am95, A)                    # maturity
		wa    <- getWeight(winf, linf, b, k, t0, A)            # weight
		Ft    <- getFishing(M, T, h1, h2)                      # fishing rate
		#tput(.PBSpopsim)
		Rt    <- getRecruitmentN(R, gamma1, sigma1, delta, A, T)      # recruitment numbers
		RBt   <- getRecruitmentB(R, gamma1, sigma1, delta, A, T, wa)  # recruitment biomass
		DetN  <- N <- matrix(nrow=A, ncol=T)                   # initial number of fish

		DetN[,1] <- getDetInitial(M, A, Rt);  Nbar <- DetN;
		N[,1]    <- getStoInitial(M, sigma2, delta, A, Rt, Nbar)

# initialize so that I can use the indexing in the loop
		Pt <- DetPt <- PBt <- DetPBt <- Ct <- DetCt <- CBt <- DetCBt <- Ibar <- StoI <- IBbar <- StoBI <- 0
		uat <- Detuat <- uBat <- DetuBat <- Pat <- PBat <- matrix(nrow=A, ncol=T)

		for(yr in 1:T) {
			Pt[yr]       <- getSelectedN(A, betaa, N, yr)        # selected population numbers
			DetPt[yr]    <- getSelectedN(A, betaa, DetN, yr) 
			PBt[yr]      <- getSelectedB(A, betaa, N, yr, wa)    # selected population biomass
			DetPBt[yr]   <- getSelectedB(A, betaa, DetN, yr, wa)
			uat[,yr]     <- getUN(A, betaa, N, yr, Pt)           # actual proportion numbers
			if (doDir) 
				.PBSpopsim$epsilon[,yr] <- .getRdir(pvec=uat[,yr],nY=1, N=dirN) # replace random normals with Dirichlets
			Detuat[,yr]  <- getUN(A, betaa, DetN, yr, DetPt)
			uBat[,yr]    <- getUB(A, betaa, N, yr, PBt, wa)      # actual proportion biomass
			DetuBat[,yr] <- getUB(A, betaa, DetN, yr, DetPt, wa)
			Ct[yr]       <- getCatch(Ft, Pt, yr)                 # catch numbers
			DetCt[yr]    <- getCatch(Ft, DetPt, yr)
			CBt[yr]      <- getCatch(Ft, PBt, yr)                # catch biomass
			DetCBt[yr]   <- getCatch(Ft, DetPBt, yr)
			Ibar[yr]     <- getDetIndex(q, yr, Pt, Ct)           # index numbers
			StoI[yr]     <- getStoIndex(tau1, yr, Ibar, A, epsilon0)
			IBbar[yr]    <- getDetIndex(q, yr, PBt, CBt)         # index biomass
			StoBI[yr]    <- getStoIndex(tau1, yr, IBbar, A, epsilon0)
			Pat[1:A,yr]  <- getProportion(A, tau2, epsilon, uat, yr, doDir)  # observed proportion numbers
			PBat[1:A,yr] <- (Pat[1:A,yr]*wa[1:A]) / sum(Pat[1:A,yr]*wa[1:A]) # obs proportion biomass
			#PBat[,yr]    <- getProportion(A, tau2, epsilon, uBat, yr, doDir) # observed proportion biomass

			if(yr < T) {                                  # number of fish
				DetN[,yr+1] <- getDetPredict(M, A, Rt, DetN, yr+1, Detuat, DetCt)
				Nbar[,yr+1] <- getDetPredict(M, A, Rt, N, yr+1, uat, Ct)      
				N[,yr+1]    <- getStoPredict(M, A, Rt, N, Nbar, yr+1, sigma2, delta)
			}
		}
		DetNB <- getBioPredict(DetN, wa, A)      # biomass of fish
		NB    <- getBioPredict(N, wa, A)
		St    <- getSpawnersN(ma, N, A, T)       # spawner numbers
		SBt   <- getSpawnersB(ma, N, A, T, wa)   # spawner biomass
		Tt    <- getTotal(N)                     # total numbers
		TBt   <- getTotal(NB)                    # total biomass

	.PBSpopsim <- c(.PBSpopsim,list(betaa=betaa, ma=ma, wa=wa, Ft=Ft, Rt=Rt, RBt=RBt, 
         DetN=DetN, Nbar=Nbar, N=N, Pt=Pt, DetPt=DetPt, PBt=PBt, 
         DetPBt=DetPBt, Ct=Ct, DetCt=DetCt, CBt=CBt, DetCBt=DetCBt, 
         Ibar=Ibar, StoI=StoI, IBbar=IBbar, StoBI=StoBI, uat=uat, 
         Detuat=Detuat, uBat=uBat, DetuBat=DetuBat, Pat=Pat, PBat=PBat, 
         DetNB=DetNB, NB=NB, St=St, SBt=SBt, Tt=Tt, TBt=TBt));
   tput(.PBSpopsim)
	} else tget(.PBSpopsim)
#browser()
	unpackList(.PBSpopsim,scope="L")
	if(!is.null(act) && act == "save") {
	# save data to an Excel file
		saveExcel(betaa, ma, wa, Ft, RBt, DetNB, NB, PBt, uBat, CBt, StoBI, PBat, SBt, TBt, A, T)
	} else {
	# plotResults(betaa, ma, wa, Ft, Rt, DetN, N, Pt, uat, Ct, StoI, Pat, St, Tt, A, T, unitType, plotType, powr, maxB, percent)
		plotResults() }
	return(invisible());
}

# ------------------------------ Plotting function ----------------------------#

# saveExcel (save data in Excel format):
# Purpose:    To save the given data in an .xls file
# Parameters: Calculated data betaa, ma, wa, Ft, Rt/RBt, DetN/DetNB, N/NB,
#               Pt/PBt, uat/uBat, Ct/CBt, StoI/StoBI, Pat/PBat, St/SBt, Tt/TBt
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years  
# Returns:    NULL
saveExcel <- function(betaa, ma, wa, fish, recruit, dbubble, sbubble, sel, uat, catch, index, pbubble, spawner, total, A, T)
{
	tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L")
	fname <- promptSaveFile(filetype=list(c(".xls", "Microsoft Excel Spreadsheet")))
	# no name provided, or they pushed cancel
	if(fname == "") { return(invisible()) }
	# if that file name already exists, delete it first
	if(file.exists(fname)) { file.remove(fname) }

	conn <- RODBC::odbcConnectExcel(fname, readOnly=FALSE)
	.excelTable(conn, betaa, "Selectivity", "Selectivity", sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, ma, "Maturity", "Maturity", sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, wa, "Weight", "Weight", sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, fish, "Fishing_mortality", "Fishing_mortality", sprintf("Y%*.3i", 3, 1:T))
	.excelTable(conn, recruit, "Recruits", "Recruits", sprintf("Y%*.3i", 3, (2-A):T))
	.excelTable(conn, dbubble, "Theoretical_ages", sprintf("Y%*.3i", 3, 1:T), sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, sbubble, "True_ages", sprintf("Y%*.3i", 3, 1:T), sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, sel, "Selected_population", "Selected_population", sprintf("Y%*.3i", 3, 1:T))
	.excelTable(conn, uat, "Selected_proportion", sprintf("Y%*.3i", 3, 1:T), sprintf("A%*.3i", 3, 1:A))    
	.excelTable(conn, catch, "Catch", "Catch", sprintf("Y%*.3i", 3, 1:T))
	.excelTable(conn, index, "Index", "Index", sprintf("Y%*.3i", 3, 1:T))
	.excelTable(conn, pbubble, "Observed_ages", sprintf("Y%*.3i", 3, 1:T), sprintf("A%*.3i", 3, 1:A))
	.excelTable(conn, spawner, "Spawners", "Spawners", sprintf("Y%*.3i", 3, 1:T))     
	.excelTable(conn, total, "Total_population", "Total_population", sprintf("Y%*.3i", 3, 1:T)) 
	.excelTable(conn, t(delta), "Random_delta", sprintf("A%*.3i", 3, 1:A), sprintf("Y%*.3i", 3, (2-A):T))
	.excelTable(conn, t(epsilon), "Random_epsilon", sprintf("A%*.3i", 3, 1:A), sprintf("Y%*.3i", 3, 1:T))
	.excelTable(conn, epsilon0[1,], "Random_epsilon0", "Epsilon0", sprintf("Y%*.3i", 3, (2-A):T)) 
	odbcClose(conn)
	return()
}

# plotResults (plot the results):
# Purpose:    Display the calculated data in the requested way.
#             If "dbubble" then display a bubble plot of the deterministic
#               number or biomass of fish in each age class
#             If "sbubble" then display a bubble plot of the stochastic number
#               or biomass of fish in each age class
#             If "pbubble" then display a bubble plot of the observed proportion 
#               of fish (in numbers or biomass) in the catch
#             If "catch" then display a bar plot of the catch in either numbers
#               or biomass
#             If "total" then display a bar plot of the total number of fish or
#               total biomass
#             If "index" then display a bar plot of the observed index
#             If "recruit" then dispaly a bar plot of the number of recruits or
#               the biomass of recruits
#             If "spawner" then display a bar plot of the spawner stock biomass
#               or the number of spawners
#             If "spawnVSrt" then display a scatter plot of the spawner stock 
#               for year t-1 against the recruitment for year t in either 
#               numbers or biomass
#             If "compare" then display a plot with lines for total, spawner,
#               and selected, and bars for catch in either numbers or biomass
#             If "fish" then display a line graph of the fishing mortality
#             Highlights a given top percentage of recruitment years in a 
#               different colour in all the graphs except the fishing mortality
# Parameters: Calculated data betaa, ma, wa, Ft, Rt/RBt, DetN/DetNB, N/NB,
#               Pt/PBt, uat/uBat, Ct/CBt, StoI/StoBI, Pat/PBat, St/SBt, Tt/TBt
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years  
#             unitType is a character to indication of the units the data is in;
#               options are "numbers" and "biomass"
#             plotType is a character to indicate which plot to show; options are
#               dbubble, sbubble, pbubble, catch, total, index, select,
#               recruit, spawner, spawnVSrt, compare, fish, save
#             powr is the powr parameter to send to the plotBubble() function
#             maxB is the maximum bubble size parameter to send to the 
#               plotBubble() function
#             percent is a number between 0 and 100 indicating the top
#               percentage to highlight in a different colour
# Returns:    NULL (invisible)
#plotResults <- function(betaa, ma, wa, fish, recruit, dbubble, sbubble, sel, uat, 
#	catch, index, pbubble, spawner, total, A, T, unitType, plotType, powr, maxB, percent)
plotResults <- function() {
	getWinVal(scope="L"); tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L");
	RtOffset <- RtOffset-1; # =A-2
	if(unitType == "biomass") {
		fish<-Ft; recruit<-RBt; dbubble<-DetNB; sbubble<-NB; sel<-PBt; uat<-uBat;
		catch<-CBt; index<-StoBI; pbubble<-PBat; spawner<-SBt; total<-TBt;
    ylabel <-      "Biomass (kg)"
    catchMain <-   "Biomass of fish caught"
    selMain <-     "Selected biomass"
    totalMain <-   "Total biomass"
    indexMain <-   "Biomass index"
    recruitMain <- "Recruitment biomass"
    spawnerMain <- "Spawner stock biomass"
    spawnVSrtX <-  "Spawner biomass in year t - 1 (kg)"
    spawnVSrtY <-  "Recruitment biomass in year t (kg)"
    compareMain <- "Biomass comparison chart"
  } else {
		fish<-Ft; recruit<-Rt; dbubble<-DetN; sbubble<-N; sel<-Pt; uat<-uat;
		catch<-Ct; index<-StoI; pbubble<-Pat; spawner<-St; total<-Tt;
    ylabel <-      "Number of fish"
    catchMain <-   "Number of fish caught"
    selMain <-     "Number of fish selected"
    totalMain <-   "Total number of fish"
    indexMain <-   "Observed index in numbers"
    recruitMain <- "Number of recruits"
    spawnerMain <- "Number of spawners"
    spawnVSrtX <-  "Number of spawners in year t - 1"
    spawnVSrtY <-  "Number of recruits in year t"
    compareMain <- "Comparison chart in numbers of fish"
  }
  
  # the logic vector that decides when a year has had high recruitment
  if(percent == 0) {
    logVecR <- rep(FALSE, times=length(recruit))
  } else if (percent == 100) {
    logVecR <- rep(TRUE, times=length(recruit))
  } else {
    logVecR <- recruit > sort(recruit)[ceiling((100-percent) / 100 * length(recruit))]
  }
	Aborn <- ((-(A-2):T)[logVecR])-1; # age born
	bline <- function(A,b=1,...) { # Assumes birth at Year A, age 0.
		n  <- length(A); x1 <- A; x2 <- rep(par()$usr[2],n); y1 <- rep(0,n);
		y2 <- b*(x2-x1) + y1; y0 <- rep(par()$usr[3],n);
		for (i in 1:n)
			lines(x=c(x1[i],x1[i],x2[i]),y=c(y0[i],y1[i],y2[i]),...) }

  logVec <- logVecR[A:(T + A - 1)]          # removes the prehistory
  ltyVec <- ifelse(logVecR, 1, 0)
  colVec <- ifelse(logVec, colHigh, colLow)
    
	if(plotType == "dbubble") {
		dbubble<-round(dbubble,4); z0<-dbubble==0; dbubble[z0]<- -0.0001;
		plotBubbles(dbubble, xlab="Year", ylab="Age class",clrs=c("grey40","white"),
			main=paste("Theoretical age of fish (",substring(ylabel,1,1),") at the start of the year",sep=""), 
			powr=powr, size=maxB, xlim=extendrange(1:T, f=0.02), ylim=extendrange(1:A, f=0.02))  
		bline(Aborn,col=colHigh)
#    for(i in -(A-2):T) {        # i is the year
#      abline(a=-(i-1), b=1, col=colHigh, lty=ltyVec[[A+i-1]], lwd=1) }

  } else if(plotType == "sbubble") {
		sbubble<-round(sbubble,4); z0<-sbubble==0; sbubble[z0]<- -0.0001;
		plotBubbles(sbubble, xlab="Year", ylab="Age class", clrs=c("grey40","white"),
			main=paste("True age of fish (",substring(ylabel,1,1),") at the start of the year",sep=""),
			powr=powr, size=maxB, xlim=extendrange(1:T, f=0.02), ylim=extendrange(1:A, f=0.02))
		bline(Aborn,col=colHigh)
#    for(i in -(A-2):T) {        # i is the year
#      abline(a=-(i-1), b=1, col=colHigh, lty=ltyVec[[A+i-1]], lwd=1) }

	} else if(plotType == "pbubble") {
		pbubble<-round(pbubble,4); z0<-pbubble==0; pbubble[z0]<- -0.0001;
		plotBubbles(pbubble, xlab="Year", ylab="Age class",clrs=c("grey40","white"),
			main=paste("Observed age of fish (",substring(ylabel,1,1),") in the catch",sep=""), 
			powr=powr, size=maxB, xlim=extendrange(1:T,f=0.02), ylim=extendrange(1:A,f=0.02))
		bline(Aborn,col=colHigh)
#    for(i in -(A-2):T) {        # i is the year
#      abline(a=-(i-1), b=1, col=colHigh, lty=ltyVec[[A+i-1]], lwd=1) }

  } else if(plotType == "catch") {  
    barplot(catch, names.arg=1:T, xlab="Year", ylab=ylabel, main=catchMain, col=colVec)
    legend("topright", legend=c(paste("Max:", round(max(catch), digits=4)), paste("Min:", round(min(catch), digits=4))), bty="n")

  } else if(plotType == "select") {
    barplot(sel, names.arg=1:T, xlab="Year", ylab=ylabel, main=selMain, col=colVec)
    legend("topright", legend=c(paste("Max:", round(max(sel), digits=4)), paste("Min:", round(min(sel), digits=4))), bty="n")
    
  } else if(plotType == "total") {
    barplot(total, names.arg=1:T, xlab="Year", ylab=ylabel, main=totalMain, col=colVec)
    legend("topright", legend=c(paste("Max:", round(max(total), digits=4)), paste("Min:", round(min(total), digits=4))), bty="n")

  } else if(plotType == "index") {
    barplot(index, names.arg=1:T, xlab="Year", ylab=ylabel, main=indexMain, col=colVec)
    legend("topright", legend=c(paste("Max:", round(max(index), digits=4)), paste("Min:", round(min(index), digits=4))), bty="n")

  } else if(plotType == "recruit") {
    colVecR <- ifelse(logVecR, colHigh, colLow)
    colVecR[1:(A-1)] <- ifelse(colVecR[1:(A-1)] == colLow, colPre, colHigh)
    barplot(recruit, names.arg=(2-A):T, xlab="Year", ylab=ylabel, main=recruitMain, col=colVecR)
    legend("topright", legend=c(paste("Max:", round(max(recruit), digits=4)), paste("Min:", round(min(recruit), digits=4))), bty="n")

  } else if(plotType == "spawner") {
    barplot(spawner, names.arg=1:T, xlab="Year", ylab=ylabel, main=spawnerMain, col=colVec)
    legend("topright", legend=c(paste("Max:", round(max(spawner), digits=4)), paste("Min:", round(min(spawner), digits=4))), bty="n")
    
  } else if(plotType == "spawnVSrt") {
    plot(spawner[1:(T-1)], recruit[(A+1):(T+A-1)], col=colVec[2:T], type="p", xlab=spawnVSrtX, ylab=spawnVSrtY, main="Recruitment VS spawner relationship", lwd=2)    

  } else if(plotType == "compare") {
    ltyVec <- ifelse(logVec, 3, 0)
    ymax <- max(total, catch, recruit, spawner)
    plot(total, type="n", ylim=c(0,ymax), xlab="Year", ylab=ylabel, main=compareMain)
    abline(v=1:T, lty=ltyVec, col=colHigh)
    lines(total, col="purple", type="l", lwd=2)
    lines(spawner, col="seagreen", type="l", lwd=2)
    lines(sel, col="darkorange", type="l", lwd=2)
    drawBars(1:T, catch, col="deepskyblue", lwd=2,width=1)
    legend("topright", legend=c("Total", "Spawner", "Selected", "Catch"), col=c("purple", "seagreen", "darkorange", "deepskyblue"), lwd=2)  
    
  } else if(plotType == "fish") {
    plot(fish, type="l", xlab="Year", ylab="F", main="Fishing mortality")
  }
  
  invisible()
    
}       

#------------------------------ Helper functions ------------------------------#
# .getRnorm (get random normal numbers):
# Purpose:   Creates and returns a matrix of random normal numbers 
#            with nA (# ages) rows and nY (# years) columns
# Arguments: nA  - numeric scalar for the number of rows (ages)
#            nY  - numeric scalar for the number of columns (years)
# Note:      In this program, the number of columns uses the same index as Rt,
#            so when accessing these values, remember to add RtOffset
.getRnorm <- function(nA, nY, nI=1) {
	if (nI==1) Rnorm <- matrix(rnorm(nA*nY), nrow=nA, ncol=nY)
	else Rnorm <- array(rnorm(nA*nY*nI), dim=c(nA,nY,nI))
	return(Rnorm) };

#----------------------------------------------------------------
# Dirichlet random sample
# Purpose:   Creates and returns a matrix of Dirichlet random proportions 
#            with nA (# ages) rows and nY (# years) columns
# Arguments: pvec - numeric vector of observed proportions (ages)
#            nY   - numeric scalar for the number of columns (years)
#            N    - effective sample size
#   Output: matrix (length(pvec) x nY)
#           each column is a Dirichlet sample that sums to 1
.getRdir <- function(pvec, nY, N) {
	# Composition operator: converts vector v to proportions
	comp <- function(v) {
		v <- v[v>0 & !is.na(v)]
		y <- abs(v)/sum(abs(v)); y; };
	pvec <- comp(pvec); # forces sum to 1
	nA   <- length(pvec);
	yvec <- rgamma(nA*nY,shape=N*pvec,scale=N);
	ymat <- matrix(yvec,nrow=nA,ncol=nY,byrow=FALSE);
	ysum <- apply(ymat,2,"sum");
	Rdir <- sweep(ymat,2,ysum,"/");	
	return(Rdir); };

#----------------------------------------------------------------
# .excelTable (Microsoft Excel table):
# Purpose:    For a given connection and a given set of data, create a dataframe
#               and save that dataframe as a table in the database. If an
#               .xls file of that name already exists then it deletes that file
#               and creates a new one, otherwise it just creates the new one
# Parameters: channel is a RODBC connection string with the .xls file name
#             dat is a list, vector, or matrix that you want to be in the 
#               table
#             tablename is a character containing the string you want to appear
#               on the Excel tabs between worksheets (tables)
#             colnam is a vector of length equal to the number of columns in dat
#               containing the column names
#             rownam is a vector of length equal to the number of rows in dat
#               containing the row names
# Returns:    NULL
.excelTable <- function(channel, dat, tablename, colnam, rownam) {
  dframe <- as.data.frame(dat)
  names(dframe) <- colnam
  rownames(dframe) <- rownam
  sqlSave(channel, dframe, tablename=tablename)
  
  return()
}

#----------------------------------------------------------------
# .workingDir (working directory):
# Purpose:    Change the working directory to a temporary directory and copy all
#               files necessary for the GUI into the new working directory from
#               the .../library/PBSref/files directory
# Parameters: None
# Returns:    Nothing
# changes the working directory to a temporary directory and copy necessary 
# files with it
.workingDir <- function () {
  # .workingDirQuit (quit the working directory)
  # Purpose:    At exiting the GUI, change the working directory back to how
  #               it was before
  # Parameters: None
  # Returns:    Nothing
	.workingDirQuit <- function() {
		closeWin(c("window","runE"))
		setwd(tget(.cwd))
		remove(list = setdiff(ls(pos=locale), tget(.cls)), pos=locale)
		return()
	}
	tput(.workingDirQuit)
	
	# save the current working directory 
	.cls <- ls(pos=locale, all.names=TRUE); tput(.cls)
	.cwd <- getwd(); tput(.cwd)
	
	# the files you wish to copy to the new working directory
	pckg <- "PBSassess"
	dnam <- "files"
	
	# find the path to the R directory on the user's computer
	rdir <- system.file(package = pckg)
	
	# get the names of all the files in the PBSref/inst directory
	wdir <- paste(rdir, "/", dnam, sep = "")
	
	# the files you want to copy to the new working directory
	fnam <- paste(wdir, list.files(wdir), sep = "/")
	
	# find the user's default temp directory
	rtmp <- tempdir()
	
	# copy all the files into the temp directory
	file.copy(fnam, rtmp, overwrite = TRUE)
	
	# set the working directory to the temp directory
	setwd(rtmp)
}

#------------------------------ Calculations ----------------------------------#

# getSelectivity (get the selectivity vector):
# Purpose:    Calculate the selectivity proportion for each age class
# Parameters: as50 is a numeric for the age at 50% selectivity
#             as95 is a numeric for the age at 95% selectivity
#             A is a numeric for the total number of age groups
# Returns:    A numeric vector of length A containing the corresponding 
#               selectivity proportions
# Note:       Formulae (2) and (3) from PBSassessDoc2.pdf
getSelectivity <- function(as50, as95, A) {
	aRange <- 1:A
	d <- log(19) / (as95 - as50)                    # Formula (2)
	betaVec <- 1 / (1 + exp(-d * (aRange - as50)))  # Formula (3)
	return(betaVec) };

#----------------------------------------------------------------
# getMaturity (get the maturity vector):
# Purpose:    Calculate the maturity proportion for each age class
# Parameters: am50 is a numeric for the age at 50% maturity
#             am95 is a numeric for the age at 95% maturity
#             A is a numeric for the total number of age groups
# Returns:    A numeric vector of length A containing the corresponding 
#               maturity proportions
# Note:       Formulae (4) and (5) from PBSassessDoc2.pdf
getMaturity <- function(am50, am95, A) {
	aRange <- 1:A
	g = log(19) / (am95 - am50)                    # Formula (4)
	mVec <- 1 / (1 + exp(-g * (aRange - am50)))    # Formula (5)
	return(mVec) };

#----------------------------------------------------------------
# getWeight (get the weight vector):
# Purpose:    Calculate the weight for each age class
# Parameters: winf is a numeric for w infinity
#             linf is a numeric for L infinity
#             b is a numeric for the VonBertlanaffey b
#             k is a numeric for the VonBertlanaffey K
#             t0 is a numeric for the VonBertlanaffey t0
#             A is a numeric for the total number of age classes
# Returns:    A numeric vector of length A containing the corresponding weights
# Note:       Formulae (6) and (7) in PBSassessDoc2.pdf
getWeight <- function(winf, linf, b, k, t0, A) {
	aRange <- 1:A
	lVec <- linf * (1 - exp(-k * (aRange - t0)))  # Formula (6) Length vector 
	wVec <- winf * (lVec[aRange] / linf) ^ b      # Formula (7)
	return(wVec) };

#----------------------------------------------------------------
# getRecruitmentN (get the recruitment vector in numbers):
# Purpose:    Calculate the number of recruits into age class 1 in each year
#               starting from year 2-A 
# Parameters: R is a numeric for the mean recruitment 
#             gamma1 is a numeric for the recruitment autocorrelation
#             sigma1 is a numeric for standard error
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years
# Returns:    A numeric vector of length (T+A-1) containing the corresponding 
#                 number of recruits
# Notes:      The range of t is from 2-A <= t <= T, but the indicies are from
#               1 <= t-(A-1) <= T-(A-1). Therefore, the written t is the index 
#               t+(A-1), call this A-1 the RtOffset
getRecruitmentN <- function(R, gamma1, sigma1, delta, A, T) {
	rVec <- getRecruitmentB(R, gamma1, sigma1, delta, A, T, 1)
	return(rVec) };

#----------------------------------------------------------------
# getRecruitmentB (get the recruitment vector in biomass):
# Purpose:    Calculate the biomass of recruits into age class 1 in each year
#               starting from year 2-A 
# Parameters: R is a numeric for the mean recruitment 
#             gamma1 is a numeric for the recruitment autocorrelation
#             sigma1 is a numeric for standard error
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years
#             wa is a numeric vector of length A containing the corresponding
#               weight
# Returns:    A numeric vector of length (T+A-1) containing the corresponding 
#                 number of recruits or their biomass
# Notes:      Formulae (9) - (11) in PBSassessDoc2.pdf
#             The range of t is from 2-A <= t <= T, but the indicies are from
#               1 <= t-(A-1) <= T-(A-1). Therefore, the written t is the index 
#               t+(A-1), call this A-1 the RtOffset
getRecruitmentB <- function(R, gamma1, sigma1, delta, A, T, wa) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L");
	RtOffset = A-1
	rVec <- R * exp((sigma1 / sqrt(1 - (gamma1 ^ 2))) * delta[1,2-A+RtOffset])         # Formula (9)
	for(i in 2:(T+RtOffset)) {
		rVec[[i]] <- R ^ (1 - gamma1) * rVec[[i-1]] ^ gamma1 * exp(sigma1 * delta[1,i]) # Formula (10)
	}
	rVec <- rVec * wa[1]   # Formula (11)
	return(rVec) };

#----------------------------------------------------------------
# getDetInitial (get the deterministic initial numbers of fish):
# Purpose:    Calculate the deterministic number of fish in each age group in 
#               year 1
# Parameters: M is a numeric for the natural mortality rate
#             A is a numeric for the total number of age classes
#             Rt is a numeric vector of length (T+A-1) containing the 
#               corresponding number of recruits
# Returns:    A numeric vector of length A containing the corresponding number
#               of fish
# Note:       Formulae (12) and (13) from PBSassessDoc2.pdf
getDetInitial <- function(M, A, Rt) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	RtOffset = A-1
	aRange <- 1:A
	NVec <- Rt[2 - aRange + RtOffset] * exp(-M * (aRange - 1))  # Formula (12)
	NVec[[A]] <- NVec[[A]] / (1 - exp(-M))                      # Formula (13)
	return(NVec) };

#----------------------------------------------------------------
# getStoInitial (get stochastic initial numbers of fish):
# Purpose:    Calculate the stochastic number of fish in each age group in
#               year 1
# Parameters: M is a numeric for the natural mortality rate
#             sigma2 is a numeric for standard error
#             A is a numeric for the total number of age classes
#             Rt is a numeric vector of length (T+A-1) containing the corresponding 
#               number of recruits
#             Nbar is a A by T numeric matrix containing the number of fish in 
#               year t of age group a
# Returns:    A numeric vector of length A conatining the corresponding number
#               of fish
# Note:       Formulae (14) and (15) from PBSassessDoc2.pdf
getStoInitial <- function(M, sigma2, delta, A, Rt, Nbar) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	RtOffset = A-1
	aRange <- 2:A;
	NVec <- Rt[1 + RtOffset]               # Formula (14)
	for(i in aRange) {                     # Formula (15)
		NVec[i] <- 1
		for(j in 2:i) {
			NVec[i] <- NVec[i] * exp(sigma2 * delta[j, j-i+1+RtOffset]) / (1 - exp(-M) + exp(-M) * exp(sigma2 * delta[j, j-i+1+RtOffset]))
		}
		NVec[i] <- Nbar[i,1] * NVec[i]
	}
	return(NVec) };

#----------------------------------------------------------------
# getSelectedN (get the number of selected fish):
# Purpose:    Calculate the selected numbers of fish at the start of the given 
#               year
# Parameters: A is a numeric for the total number of age classes
#             betaa is a numeric vector of length A containing the corresponding 
#               selectivity proportions
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             yr is a numeric for the year
# Returns:    A numeric for the selected number of fish in the given year
getSelectedN <- function(A, betaa, N, yr) {
	sVec <- getSelectedB(A, betaa, N, yr, rep(1, times=A))
	return(sVec) };

#----------------------------------------------------------------
# getSelectedB (get the biomass of selected fish):
# Purpose:    Calculate the selected biomass at the start of the given yaer
# Parameters: A is a numeric for the total number of age classes
#             betaa is a numeric vector of length A containing the corresponding
#               selectivity proportions
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             yr is a numeric for the year
#             wa is a numeric vector of length A containing the corresponding
#               weight
# Returns:    A numeric for the selected biomass in the given year
# Note:       Formula (17) from PBSassessDoc2.pdf
getSelectedB <- function(A, betaa, N, yr, wa) {
	aRange <- 1:A
	sVec <- sum(betaa[aRange] * wa[aRange] * N[aRange,yr])   # Formula (17)
	return(sVec) };

#----------------------------------------------------------------
# getUN (get the exploitable proportion in numbers):
# Purpose:    Calculate the exploitable proportion of fish in numbers for the
#               given year
# Parameters: A is a numeric for the total number of age groups
#             betaa is a numeric vector of length A containing the corresponding
#               selectivity proportions
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             yr is a numeric for the year
#             Pt is a numeric vector of length T containing the selected number
#               of fish
# Returns:    A numeric vector containing the exploitable proportion moment for
#               the corresponding age group in the given year in numbers
getUN <- function(A, betaa, N, yr, Pt) {
	uVec <- getUB(A, betaa, N, yr, Pt, rep(1, times=A))
	return(uVec) }

#----------------------------------------------------------------
# getUB (get the exploitable proportion in biomass):
# Purpose:    Calculate the exploitable proportion of fish in biomass for the
#               given year
# Parameters: A is a numeric for the total number of age groups
#             betaa is a numeric vector of length A containing the corresponding
#               selectivity proportions
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             yr is a numeric for the year
#             PBt is a numeric vector of length T containing the selected 
#               biomass
#             wa is a numeric vector of length A containing the corresponding
#               weight
# Returns:    A numeric vector containing the exploitable proportion moment for
#               the corresponding age group in the given year in biomass
# Note:       Formula (19) from PBSassessDoc2.pdf
getUB <- function(A, betaa, N, yr, PBt, wa) {
	aRange <- 1:A
	uVec <- betaa[aRange] * N[aRange, yr] * wa[aRange] / PBt[yr]  # Formula (19)
	return(uVec) }

#----------------------------------------------------------------
# getFishing (get the fishing mortality):
# Purpose:    Calculate the fishing mortailty for each year
# Parameters: M is a numeric for the natural mortality rate
#             T is a numeric for the total number of years
#             h1 is a numeric for the first fishing rate scalar
#             h2 is a numeric for the second fishing rate scalar
# Returns:    A vector of length T containing the corresponding fishing 
#               mortality
# Note:       Formula (8) from PBSassessDoc2.pdf
getFishing <- function(M, T, h1, h2) {
	yr <- 1:T
	Ft <- M / (2 * T) * (2 * h1 * (T - abs(2 * yr - T)) + h2 * ( abs(2 * yr - T) + 2 * yr - T))  # Formula (5)
	return(Ft) }

#----------------------------------------------------------------
# getCatch (get the catch):
# Purpose:    Calculate the number of fish caught in the given year or their
#               biomass depending on the units of Pt
# Parameters: Ft is a vector of length T containing the corresponding fishing 
#               mortality
#             Pt is a numeric vector of length T containing the corresponding
#               selected population numbers or biomass
#             yr is a numeric for the year
# Returns:    A numeric for the number of fish caught in the given year or their
#               biomass
# Note:       Formulae (20) and (21) from PBSassessDoc2.pdf
getCatch <- function(Ft, Pt, yr) {
	Ct <- (1 - exp(-Ft[yr])) * Pt[yr]   # Formula (20) and (21)
	return(Ct) }

#----------------------------------------------------------------
# getDetIndex (get the deterministic abundance index):
# Purpose:    Calculate the deterministic abundance index for the given year
#               in numbers or biomass depending on the units of Pt and Ct
# Parameters: q is a numeric for the catchability
#             yr is a numeric for the year
#             Pt is a numeric vector of length T containing the corresponding
#               selected number of fish or their biomass
#             Ct is a numeric vector of length T containing the corresponding
#               number of fish caught or their biomass
# Returns:    A numeric for the deterministic abundance index for the given year
#               in numbers or in biomass
# Note:       Formula (22) and (24) from PBSassessDoc2.pdf
getDetIndex <- function(q, yr, Pt, Ct) {
	It <- q * (Pt[yr] - Ct[yr] / 2)   # Formula (22) and (24)
	return(It) }

#----------------------------------------------------------------
# getStoIndex (get the stochastic index):
# Purpose:    Calculate the stochastic abundance index for the given year in 
#               either numbers or biomass depending on the units of Ibar
# Parameters: tau1 is a numeric for standard error
#             yr is a numeric for the year
#             Ibar is a numeric vector of length T containing the corresponding
#               deterministic index in either numbers or biomass
#             A is the total number of age classes
# Returns:    A numeric for the stochastic index for the given year in either
#               numbers or biomass
# Note:       Formula (23) and (25) from PBSassessDoc2.pdf
getStoIndex <- function(tau1, yr, Ibar, A, epsilon0) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	RtOffset = A - 1
	It <- Ibar[yr] * exp(tau1 * epsilon0[yr + RtOffset])   # Formula (23) and (25)
	return(It) }

#----------------------------------------------------------------
# getDetPredict (get the deterministic prediction of number of fish):
# Purpose:    Calculate the deterministic number of fish for the given year
# Parameters: M is a numeric for the natural mortality rate
#             A is a numeric for the total number of age classes
#             Rt is a numeric vector of length T containing the corresponding
#               recruitment
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             yr is a numeric for the year
#             uat is a A by T numeric matrix containing the exploitable 
#               proportion moment for age group a in a year t
#             Ct is a numeric vector of length T containing the corresponding
#               number of fish caught
# Returns:    A numeric vector of length A containing the deterministic number
#               of fish for the given year
# Note:       Formulae (D.7), (D.8), and (D.9) from PBSassessDoc1.pdf
getDetPredict <- function(M, A, Rt, N, yr, uat, Ct) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	RtOffset = A - 1
	aRange <- 2:(A-1);
	nVec <- Rt[yr + RtOffset]   # Formula (D.7)
	nVec[aRange] <- exp(-M) * (N[aRange-1, yr-1] - uat[aRange-1, yr-1] * Ct[yr-1])  # Formula (D.8)
	nVec[A] <- exp(-M) * (N[A-1, yr-1] + N[A, yr-1] - (uat[A-1,yr-1] + uat[A, yr-1]) * Ct[yr-1])  # Formula (D.9)
	return(nVec) }

#----------------------------------------------------------------
# getStoPredict (get the stochastic prediction of number of fish):
# Purpose:    Calculate the stochastic number of fish for the given year
# Parameters: M is a numeric for the natural mortality
#             A is a numeric for the total number of age classes
#             Rt is a numeric vector of length T containing the corresponding
#               recruitment
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             Nbar is a A by T numeric matrix containing the deterministic 
#               number of fish in year t of age group a
#             yr is a numeric for the year
#             sigma2 is a numeric for standard error
# Returns:    A numeric vector of length A containing the stochastic number
#               of fish for the given year
# Note:       Formulae (32) and (33) from PBSassessDoc2.pdf
getStoPredict <- function(M, A, Rt, N, Nbar, yr, sigma2, delta) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	RtOffset = A - 1
	aRange <- 2:A;
	N[1, yr] <- Rt[yr + RtOffset]   # Formula (32)
	# Formula (33)
	N[aRange,yr] <- Nbar[aRange,yr] * exp(sigma2 * delta[aRange,yr+RtOffset]) / (1 - exp(-M) + exp(-M) * exp(sigma2 * delta[aRange,yr+RtOffset]))
	return(N[,yr]) }

#----------------------------------------------------------------
# getBioPredict (get the biomass of the predicted ages):
# Purpose:    Given the matrix of number of fish at age a in year t, calculate
#               the biomass of the fish of age a in year t
# Parameters: N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             wa is a numeric vector of length A containing the corresponding
#               weights
#             A is a numeric for the number of age classes
# Returns:    An A by T numeric matrix containing the corresponding biomasses
# Note:       Formula (36) from PBSassessDoc2.pdf
getBioPredict <- function(N, wa, A) { 
	mult <- function(N, wa, A) {
		aRange <- 1:A; nVec <- N * wa; }
	nMat <- apply(N, 2, mult, wa=wa, A=A)   # Formula (36)
	return(nMat) }

#----------------------------------------------------------------
# getProportion (get the observed proportion of fish in the catch):
# Purpose:    Calculate the stochastic observed proportion of fish in the catch 
#               in the given year in numbers or biomass depending on the units
#               of uat
# Parameters: A is a numeric for the total number of age classes
#             tau2 is a numeric for standard error
#             uat is a A by T numeric matrix containing the exploitable 
#               proportion moment for the corresponding age group in the given 
#               year in numbers or in biomass
#             yr is a numeric for the year
# Returns:    A numeric vector of length A containing the corresponding 
#               stochastic observed proportion of fish in the catch in the given
#               year in numbers or biomass
# Note:       Formulae (26)-(29)
getProportion <- function(A, tau2, epsilon, uat, yr, doDir) {
	#tget(.PBSpopsim); unpackList(.PBSpopsim,scope="L"); 
	aRange <- aRange2 <- 1:A;
	# Formula (26) and (28)
	if (doDir)
		pVec <- epsilon[aRange,yr]
	else {
		xVec <- log(uat[aRange,yr]) + tau2 * epsilon[aRange,yr] - 
			     (1/A) * sum(log(uat[aRange2,yr]) + tau2 * epsilon[aRange2,yr])
		pVec <- exp(xVec[aRange]) / sum(exp(xVec[aRange])) }  # Formula (27) and (29)
	return(pVec) }

#----------------------------------------------------------------
# getSpawnersN (get the spawners in numbers):
# Purpose:    Calculate the number of spawners for each year
# Parameters: ma is a numeric vector of length A containing the corresponding
#               maturity
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years
# Returns:    A numeric vector of length T containing the corresponding number 
#               of spawners
getSpawnersN <- function(ma, N, A, T) {
	sVec <- getSpawnersB(ma, N, A, T, rep(1, times=A))
	return(sVec) };

#----------------------------------------------------------------
# getSpawnersB (get the spawners in biomass):
# Purpose:    Calculate the spawner stock biomass for each year
# Parameters: ma is a numeric vector of length A containing the corresponding
#               maturity
#             N is a A by T numeric matrix containing the number of fish in year
#               t of age group a
#             A is a numeric for the total number of age classes
#             T is a numeric for the total number of years
#             wa is a numeric vector of length A containing the corresponding
#               weight
# Returns:    A numeric vector of length T containing the corresponding spawner
#               stock biomass
# Note:       Formula (36)
getSpawnersB <- function(ma, N, A, T, wa) {
	aRange <- 1:A;  yr <- 1:T;  sVec <- 0;
	for(i in yr) {
		sVec[i] <- sum(ma[aRange] * wa[aRange] * N[aRange,i])   # Formula (36)
	}
	return(sVec) }

#----------------------------------------------------------------
# getTotal (get fish totals):
# Purpose:    Calculate the total number of fish or the total biomass for each
#               year depending on the units of N
# Parameters: N is a A by T numeric matrix containing the number of fish in year
#               t of age group a or their biomass
# Returns:    A numeric vector of length T containing the corresponding total
#               number of fish or their biomass
# Note:       Formula (37) and (38)
getTotal <- function(N) {
	# Formula (37) and (38)
	return(apply(N, 2, sum)) }


#------------------------------ Error checking --------------------------------#
# validParam (valid parameters):
# Purpose:     Check whether the parameters supplied in the GUI are valid or not.
#              If invalid, display an appropriate error message in the R console
#                and clear the invalid field
#              If it is a correctable field, corrects it in the GUI and does not
#                flag it as invalid
# Parameters:  None
# GUI inputs:  Assessment parameters as50, as95, am50, am95, winf, linf, b, t0,
#                k, h1, h2, sigma1, sigma2, tau1, tau2, A, T, M, R, q, gamma1
# Returns:     TRUE if the parameters were valid
#              FALSE otherwise
# GUI outputs: Clears invalid fields
#              Corrects correctable fields
validParam <- function()
{
  getWinVal(scope="L")
                    
  isValid <- TRUE                                                     
  changes <- list()
  
  if(is.na(percent)) {
    print(paste("Hightlight top percentage has been changed to 0"))
    changes$percent <- 0
    assign("percent", 0, pos=parent.frame(1))
  } else if((percent > 100) || (percent < 0)) {
    print("Highlight top percentage must be between 0% and 100%")
    changes$percent <- NA
    isValid <- FALSE
  }

  if((am50 < 1) || (am50 > am95) || (is.na(am50)) || (is.na(am95))) {
    print("Age at 95% maturity must be greater than age at 50% maturity")
    isValid <- FALSE
  }
  if((as50 < 1) || (as50 > as95) || (is.na(as50)) || (is.na(as95))) {
    print("Age at 95% selectivity must be greater than age at 50% selectivity")
    isValid <- FALSE
  }
  if(is.na(winf)) {
    print("w infinity cannot be blank")
    changes$winf <- NA
    isValid <- FALSE
  }
  if((is.na(linf)) || (linf < 0)) {
    print("L infinity must be greater than 0")
    changes$linf <- NA
    isValid <- FALSE
  }
  if((is.na(b)) || (b <= 2) || (b >= 4)) {
    print("b must be between 2 and 4")
    changes$b <- NA
    isValid <- FALSE
  }  
  if(is.na(t0)) {
    print("t0 cannot be blank")
    changes$t0 <- NA
    isValid <- FALSE
  }
  if(is.na(k)) {
    print("K cannot be blank") 
    changes$k <- NA
    isValid <- FALSE
  }
  if((is.na(h1)) || (h1 < 0)) {
    print(paste("h1 value", h1, "has been changed to 0"))
    changes$h1 <- 0
    assign("h1", 0, pos=parent.frame(1))
  }
  if((is.na(h2)) || (h2 < 0)) {
    print(paste("h2 value", h2, "has been changed to 0"))
    assign("h2", 0, pos=parent.frame(1))
    changes$h2 <- 0
  }
  if(is.na(sigma1)) {
    assign("sigma1", 0, pos=parent.frame(1))
    changes$sigma1 <- 0
  }
  if(is.na(sigma2)) {
    assign("sigma2", 0, pos=parent.frame(1))
    changes$sigma2 <- 0
  }
  if(is.na(tau1)) {
    assign("tau1", 0, pos=parent.frame(1))
    changes$tau1 <- 0
  }
  if(is.na(tau2)) {
    assign("tau2", 0, pos=parent.frame(1))
    changes$tau2 <- 0
  }
  if((is.na(A)) || (A <= 2)) {
    print("Number of age groups must be greater than 2")
    changes$A <- NA
    isValid <- FALSE
  }
  if((is.na(T)) || (T <= 2)) {
    print("Number of years must be greater than 2")
    changes$T <- NA
    is.Valid <- FALSE
  }
  if((is.na(M)) || (M <= 0) || (M >= 1)) {
    print("Natural mortality rate must be between 0 and 1")
    changes$M <- NA
    is.Valid <- FALSE
  }
  if((is.na(R)) || (R <= 0)) {
    print("Recruitment must be greater than 0")
    changes$R <- NA
    is.Valid <- FALSE
  }
  if(is.na(q)) {
    print("q cannot be blank")
    changes$q <- NA
    isValid <- FALSE
  }
  if((is.na(gamma1)) || (abs(gamma1) >= 1)) {
    print("gamma must be between -1 and 1")
    changes$gamma1 <- NA
    isValid <- FALSE
  }
  if((is.na(maxB)) || (maxB <= 0)) {
    print("Max bubble size must be greater than 0")
    changes$maxB <- NA
    isValid <- FALSE
  }
  if(is.na(powr)) {
    print("Bubble power cannot be blank")
    changes$powr <- NA
    isValid <- FALSE
  }
  setWinVal(changes)
  return(isValid)
  
}

#===================================================================================================
if (!require(PBSmodelling, quietly=TRUE)) stop("The `PBSmodelling` package is required for this example")
if (!require(RODBC, quietly=TRUE)) stop("The `RODBC` package is required for this example")
  resetGraph(); if (exists(".PBSpopsim",where=.PBSmodEnv)) rm(.PBSpopsim,pos=.PBSmodEnv)
  createWin("PopSimWin.txt"); calcAssess();

}) # end local scope
