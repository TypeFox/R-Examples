SPACECAP <- function() {
# Modified by Mike Meredith: lik1, if(Mb), GOF moved,  misc. speed tweaks,
#   scaling and neighbours, MMDM for intial sigma, corrected p2,
#   reports in info txt file.
# 4/25/2014 -- Andy changed 12.5 to 1 should resolve sigma scaling problem
# 27 Apr 2014 -- Mike changed scaling of output sigma
# 26 Jun 2014 -- Arjun added pdf of help files. Isolated "windows" and "unix-alike" functions 
#   calling pdf of the HELP file. Details on priors in help files. Other minor comments. 

  # NB: SPACECAP "Depends" on tcltk and coda (see the DESCRIPTION file).
  # That means that tcltk and coda are automatically loaded before SPACECAP package.
  #  Individual functions should not try to reload packages which are already loaded,
  #  and the CRAN checks generate a "Note" if they do.
  # You only need to load tcltk and coda packages if you are running this code
  #  without first loading the SPACECAP package.
  # require(tcltk)
  # library(coda)

  tclRequire("Tktable", warn=FALSE)
    # if loading fails, only consequence is that summary table does not appear in the GUI

#################### GUI related START ######################
  # Create a top level window from a tkwidget
  tt <- tktoplevel()

  # Define Global variables
  locidso <- NULL
  locso <- NULL
  grid2500 <- NULL
  Global_iterNum <- 2500
  Global_burn <- 500
  Global_skip <- 4
  Global_nz <- 250
  Global_acSize <- 0.336
  statusText <- tclVar(" ")
  shouldIStop <- FALSE

  ## Action function to read animal capture details file
  readFile1 <- function ()
  {
    file <- tclvalue(tkgetOpenFile())
    if (!length(file)) return()
    locidso <<- read.csv(file)
    tkconfigure(capt, state="normal")
    tkdelete(capt, "0", "insert")
    tkinsert(capt, "0", file)
    tkconfigure(capt, state="disable")
    tkfocus(tt)
  }

  ## Action function to read trap deployment details file
  readFile2 <- function ()
  {
    file <- tclvalue(tkgetOpenFile())
    if (!length(file)) return()
    locso <<- read.csv(file)
    tkconfigure(trap, state="normal")
    tkdelete(trap, "0", "insert")
    tkinsert(trap, "0", file)
    tkconfigure(trap, state="disable")
    tkfocus(tt)
  }

  ## Action function to read potential activity centers file
  readFile3 <- function ()
  {
    file <- tclvalue(tkgetOpenFile())
    if (!length(file)) return()
    grid2500 <<- read.csv(file)
    tkconfigure(grids, state="normal")
    tkdelete(grids, "0", "insert")
    tkinsert(grids, "0", file)
    tkconfigure(grids, state="disable")
    tkfocus(tt)
  }

  ## Define the frame structure for SPACECAP
  one <- tkframe(tt, width = 800, height = 250, bg= "light blue",relief="groove",borderwidth=2)
  tkpack(one, side = "top")
  two <- tkframe(tt, width = 800, height = 300, bg= "white",relief="groove",borderwidth=2)
  tkpack(two, side = "top")
  three <- tkframe(tt, width = 800, height = 40, bg= "white",relief="groove",borderwidth=2)
  tkpack(three, side = "bottom")

  # Sub-Frames
  fileFrame <- tkframe(one, width = 300, height = 250, bg= "light blue",relief="groove",borderwidth=2)
  tkpack(fileFrame, side = "left")
  choiceFrame <-tkframe(one, width = 250, height = 250, bg= "light blue",relief="groove",borderwidth=2)
  tkpack(choiceFrame, side = "left")
  valueFrame <-tkframe(one, width = 250, height = 250, bg= "light blue",relief="groove",borderwidth=2)
  tkpack(valueFrame, side = "right")

  # Define status window
  statusWin <- tktext(three,height=5)
  scr <- tkscrollbar(three, command=function(...) tkyview(statusWin,...))
  tkconfigure(statusWin, yscrollcommand=function(...) tkset(scr,...))
  tkpack(scr, side="right", fill="y")
  tkpack(statusWin, fill="both", expand=TRUE)
#################### GUI related END ######################

  ## Function for reading input data and creating 3D data structures from it
  readData <- function() {
    nID = max(locidso[,2])
    nSO = (dim(locso)[2]) - 3
    nLOC = dim(locso)[1]

    # Sort locso by location ID, as following code assumes row no. = location ID
    locso <- locso[order(locso[, 1]), ]

    ### Function for getting capture values in a 3 dimensional ID x SO x LOC format
    makeData3d = function()     {
      data3d = structure(rep(0, times=nID*nSO*nLOC), .Dim=c(nID, nSO, nLOC))
      len = length(locidso[,1])
      for (i in 1:len)
      {
        loc = locidso[i,1]
        id = locidso[i,2]
        so = locidso[i,3]
        data3d[id,so,loc] = 1
      }
      data3d
    }
    ### Function for getting the deployment values in a LOC x SO format
    makeMask3d = function() {
      ## In the deployment file, the SO start from column 4
      posFirstSO = 4
      posLastSO = nSO+3
      mask3d = as.matrix(locso[,posFirstSO:posLastSO])

      mask3d
    }

    tiger3dData <- structure(list(makeData3d(), makeMask3d()),
      .Names = c("data3d", "mask3d"))
    # grid2500 <- structure(c(grid2500[,1], grid2500[,2], grid2500[,3]), .Dim = c(length(grid2500[,1]), 3), .Dimnames = list(c(1:length(grid2500[,1])), c("X_Coord", "Y_Coord", "HABITAT")))  # a complicated way to change names
    ctLocs <- locso[, 2:3]
    colnames(ctLocs) <- c("x", "y")
    rownames(ctLocs) <- c(1:nLOC)

    list(tiger3dData=tiger3dData, grid2500=grid2500, ctLocs=ctLocs)
  } ## End of readData function

  ####  SCRd.fn  ##################################################################
  ## The main script for running the Spacecap analysis
  ####
  SCRd.fn <- function(ni=52000, burn=2000, skip=50,
                      bsigma=1, Mb=0, nz=450, dexp=2)
  {
    # Feb 23 -- andy tinkered with activity center updating
    # Feb 24 -- changed behavioral response to regression variable
    # Sept 26 -- work on including Mb.

    # ni = number of iterations, total
    # burn = number to discard
    # skip = thin rate (i.e., keep one in skip iterations)
    # bsigma = 0 fits non-spatial model
    # nz = number of "all zero" encounter histories to add
    # dexp = 2 for halfnormal, 1 for exponential detection function

    # Required data objects are the encounter array Y (bindary observations)
    # which is nind x T x ntraps
    # MASK which is ntraps x reps (it gets transposed though)
    # traplocs = coordinates of trap locations

#################### GUI related START ######################
    cat("Starting Analysis...", fill=TRUE)
    statusText <<- "\nStarting analysis\n"
    tkinsert(statusWin, "end", statusText)
#################### GUI related END #######################

    # Read data from input files
    hold <- readData()

    tiger3dData <- hold$tiger3dData
    grid2500 <- hold$grid2500
    ctLocs <- hold$ctLocs

    Y <- tiger3dData$data3d          # nind x nT x ntraps
    MASK <- t(tiger3dData$mask3d)    # want this to be rep x traps
    # traplocs <- as.matrix(ctLocs$grid)
    traplocs <- as.matrix(ctLocs)
    GRID <- grid2500

################### Output related START ####################
    gridsub=subset(grid2500, grid2500[,3]>0)
    # nhrc=sum(grid2500[,3])
    nhrc <- nrow(gridsub)
################### Output related END ######################

    # what are the dimensions of the problem
    nind <- dim(Y)[1]
    nT <- dim(Y)[2]
    M <- nind + nz   # total size of data set
    ntraps <- nrow(traplocs)

    ## following lines scale coordinate system to range from 0 to max. 10
    ## Using raw numbers in metres produces huge values for squared distances, and
    ## exp(dist squared) risks hitting Inf.
    ## Changed by MM
    G0 <- GRID[GRID[,3] == 1, -3]  # XY coords for good-habitat pixels
    mincoord <- apply(G0, 2, min) # lower left corner of study area
    Range <- max(apply(G0, 2, max) - mincoord)/10 # Change 10 for a different max.
    G <- scale(G0, center=mincoord, scale=rep(Range,2))
    traplocs <- scale(traplocs, center=mincoord, scale=rep(Range,2))
    nG <- nrow(G)

    ### Data processing -- this block of code determines a neighborhood
    ### for every pixel. That information is used in the MCMC updating
    ### of activity centers. Some tuning would be required

    NN <- matrix(NA, nrow=nG, ncol=400)

    # RAD is a number related to the MCMC for the activity centers.
    # Pixels of good habitat within a radius of RAD are counted as neighbours
    # Here, RAD is based on width of pixels so that each has approx 310 neighbours.
    # Some latitude for error is needed, as rounding of the entry of Global_acSize
    # by the user can have a large effect on number of neighbours. MM's modification

    gridres <- sqrt(Global_acSize)*1000/Range #  pixel width in scaled units
    RAD <- gridres*10 # Change 10 to get more or fewer neighbours
    for(i in 1:nG){
      od <- sqrt( (G[i,1]-G[,1])^2  +  (G[i,2]-G[,2])^2  )
      od <- (1:length(od))[od < RAD]
      NN[i,1:length(od)] <- od
    }
    numnn <- apply(!is.na(NN),1,sum)
    NN <- NN[, 1:max(numnn)]

    ### does the data augmentation
    Yaug <- array(0, dim=c(nind+nz, nT, ntraps))
    for(j in 1:nind){
      Yaug[j, 1:nT, 1:ntraps] <- Y[j, 1:nT, 1:ntraps]
    }

    # `e2dist` <- function (x, y) ....
    # Moved to its own file so now accessible outside SCRd.fn

    ## this block of code picks starting coordinates for each individual
    ## and MMDM as starting value for sigma ## MM
    centers1 <- rep(NA,nind)
    mdm <- NULL
    for(i in 1:nind){
      tt <- Yaug[i,,]
      tt <- col(tt)[tt == 1]   # which traps was he captured in
      xxx <- traplocs[tt,]  ## coordinates of those traps  ### add drop=FALSE
      av.coord <- apply(matrix(xxx, ncol=2), 2, mean) ## takes average coordinate
      dvec <- as.vector(e2dist(matrix(av.coord,ncol=2), G))  # finds closest grid pt
      centers1[i] <- (1:length(dvec))[dvec==min(dvec)][1]   # that is initial loc
      tt1 <- unique(tt)
      if(length(tt1) > 1)
        mdm <- c(mdm, max(e2dist(traplocs[tt1,], traplocs[tt1,])))
    }
    mmdm <- mean(mdm) # Mean Maximum Distance Moved
    # uncaptured guys need centers too.....
    centers2 <- sample(1:nG, M-nind, replace=TRUE)
    centers <- c(centers1, centers2)
    S <- G[centers,]   # initial locations for all M individuals

    # create "Data" vector but with trap mask information
    msk <- MASK
    msk2 <- array(NA, c(nind+nz, nT, ntraps))
    for(i in 1:(nind+nz)){
      msk2[i, 1:nT, 1:ntraps] <- msk[1:nT, 1:ntraps]
    }
    msk2 <- as.vector(msk2)

    # create covariate of previous capture
    prevcap <- array(0,c(nind+nz,nT,ntraps))
    for(i in 1:(nind)){
      for(j in 1:ntraps){
        tmp <- Yaug[i, 1:nT, j]
        if(any(tmp == 1)){
          fst <- min( (1:nT)[tmp == 1] )
          if(fst < nT)
            prevcap[i, (fst+1):nT, j] <- 1
        }
      }
    }
    prevcap <- as.vector(prevcap)

    ## vectorize all the data objects
    arr.trues <- array(TRUE, c(nind+nz,nT,ntraps))
    idx <- which(arr.trues, arr.ind = TRUE)
    y <- as.vector(Yaug)
    y <- y[msk2==1]
    prevcap <- prevcap[msk2==1]   #### + 1   ### add 1 if want 1/2 otherwise dummy
    indid <- idx[msk2==1,1] ## [AMG] Individual IDs
    repid <- idx[msk2==1,2] ## [AMG] Replicate/Sampling Occasion IDs
    trapid <- idx[msk2==1,3] ## [AMG] Trap IDs

    ## starting values of parameters and other stuff needed, including utility funcs
    clogloginvpart1<-function(lp){
    # goes from linear predictor, lp, to log(1-p)
      -1*exp(lp)
    }
    ##clogloginvpart2 function never used ## MM

    trapgridbig <- traplocs[trapid,]   # stretches out the trap coord matrix
    y1 <- y == 1
    c1 <- (S[indid,1] - trapgridbig[,1])^2
    c2 <- (S[indid,2] - trapgridbig[,2])^2

    # sigma <- mmdm/2   ## MM
    sigma <- ifelse(bsigma == 1, mmdm/2, 1) ## bsigma == 0 is nonspatial model
    loglam0 <- log(0.018)
    beta <- 0
    lam0 <- exp(loglam0)
    psi <- 0.6  # not a good starting values
    pz <- 0.5
    z <- c(rep(1, nind), rbinom(nz, 1, pz)) ## [AMG] pz will be replaced by probz later during updating. I am 
          									## simply sticking with the definition of psi. To not be confused with 
          									##the paper Royle et al (2009) 

################ Output related START ###################
    out<-matrix(NA,nrow=(ni-burn)/skip,ncol=5)
    dimnames(out)<-list(NULL,c("sigma","lam0","beta","psi","N"))
    zout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
    Sout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
    gof.new<-gof.data<-rep(NA,(ni-burn)/skip)
######### Output related END ###############
######### GUI related START ################
    if(shouldIStop==TRUE) {
      shouldIStop <<- FALSE;
      statusText <<- "Analysis stopped\n"
      tkinsert(statusWin, "end", statusText)
      stop(call.=FALSE, "Analysis stopped")
    }
    #Progress Bar
    pb <- tkProgressBar(title = "SPACECAP Progress Bar", min = 0, max = ni, width = 600)
    statusText <<- "Burn-in in progress\n"
    tkinsert(statusWin, "end", statusText)
######### GUI related END ##################

    # m <- 1 # moved to just before main iterations loop
    ### Headings for console output: ## MM
    cat ("    iter    sigma   lam0   beta    psi N\n")

    ts <- format(Sys.time(), "%y%m%d_%H%M%S")  ## moved here by MM
    folderName <- paste("output_", ts, sep="")
    dir.create(path=folderName)
    info.file <- file.path(folderName, paste("Info_",ts,".txt", sep=""))  ## info file added by MM
    # SPACECAP does not keep track of data file names so can't add those
    start.time <- Sys.time()
    info.lines <- c("Analysis with SPACECAP 1.1.0",
      format(Sys.time(), "%a %d %b %Y"), "",
      paste("Area of habitat pixel:", Global_acSize, "sq km"),
      # paste("Distance scaling: 1 unit =", Range, "m"),
      "\nModel selected:",
      paste("\tTrap response", ifelse(Mb, "present,", "absent,")),
      paste(ifelse(bsigma, "\tSpatial", "\tNonspatial"), "Capture-Recapture,"),
      paste(ifelse(dexp==2, "\tHalf-normal", "\tNegative Exponential"),
        "detection function,"),
      "\tBernoulli detection process",
      "\nMCMC simulation settings:",
      paste("\tIterations:", ni),
      paste("\tBurn-in:", burn),
      paste("\tThinning:", skip),
      paste("\tNumber of values saved:", (ni - burn) / skip),
      paste("\tData augmentation:", nz, "Total size of dataset:", M),
      "\nStarting values:",
      paste("sigma = ", ifelse(bsigma == 1, sigma*Range/sqrt(dexp), NA), ", lam0 = ", lam0,
        ", beta = ", beta, ", psi = ", psi, sep=""),
      paste("\nStarted at:", start.time))
    writeLines(info.lines, info.file)

## ================= MAIN LOOP ========================#####
    m <- 1 # This is the counter for rows of output.
    for(i in 1:ni)  {

######### GUI related START ##################
      if(shouldIStop==TRUE) {
        shouldIStop <<- FALSE;
        statusText <<- "Analysis stopped\n"
        tkinsert(statusWin, "end", statusText)
        stop(call.=FALSE, "Analysis stopped")
      }
      if(i%%(ni/100)==0)  ## no need to update progress bar every iteration ## MM
        setTkProgressBar(pb, i, label=paste(floor((i*100)/ni),"% completed"))
######### GUI related END ##################

      # PART 1 OF THE MCMC ALGORITHM UPDATES THE REGRESSION PARAMETERS. FOR THIS MODEL
      # THE REGRESSION PARAMETERS ARE (1) INTERCEPT (2) EFFECT OF PREVIOUS CAPTURE
      # (BEHAVIORAL RESPONSE) (3) THE SPATIAL PARAMETER "sigma"
      ### Updating parameters here should only involve guys with z = 1 (i.e., members of the population)

      lp <- loglam0 + Mb*beta*prevcap - ((1*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5)) 
                                         ## [AMG] lp<-This is cloglog(pi(ij))

      lik1 <- log( expm1(exp(lp[y1]))) ## replacing lik1<- log( expm1(exp(lp))) ## MM
      lik2 <- clogloginvpart1(lp)
      llvector <- lik2
      llvector[y1] <- llvector[y1]+ lik1 ## replacing llvector[y1]<- llvector[y1]+ lik1[y1] ## MM # [AMG]This is the 
      									 ## log(argument), where "argument" is the original probability of 	
      									 ##	observation, pi(ij) in the paper Royle et al 2009. The argument is the
      									 ## likelihood and llvector represents the loglikelihood vector. 

      ## 4/27/2014 changed the tuning parameter to be log( mmdm/4 ) here
      if(bsigma == 1) {     # only update sigma for spatial models
        sigmac <- rnorm(1, sigma, mmdm/4)
        if(sigmac > 0) {
          lpc<-  loglam0 + Mb*beta*prevcap - ((1*bsigma)/(sigmac^2))*((c1+c2)^(dexp*0.5))
          lik1c<- log(expm1(exp(lpc[y1])))  ## lik1c<- log(expm1(exp(lpc)))##MM
          lik2c<- clogloginvpart1(lpc)
          llvector.new<- lik2c
          llvector.new[y1]<- llvector.new[y1]+ lik1c
          if(runif(1)<exp(sum( rowsum(llvector.new-llvector,indid)[z==1]))) {
            sigma <- sigmac
            llvector <- llvector.new
          }
        }
      }
      ### May 3 2014 made loglam0 update independently of sigma
      ### regarding choice of tuning parameter, Andy thinks loglam0 should be in the ballpark of -3 to -2 and
      ## uses a tuning parameter about 10% of the magnitude (this is completely subjective/arbitrary)
      loglam0c <- rnorm(1,loglam0,.2)
      lpc <- loglam0c + Mb*beta*prevcap - ((1*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5))
      lik1c <- log(expm1(exp(lpc[y1])))  ## lik1c<- log(expm1(exp(lpc)))##MM
      lik2c <- clogloginvpart1(lpc)
      llvector.new <- lik2c
      llvector.new[y1] <- llvector.new[y1]+ lik1c  ## llvector.new[y1]<- llvector.new[y1]+ lik1c[y1] ## MM
      if(  runif(1) < exp(sum(rowsum(llvector.new - llvector, indid)[z==1]))) {
        loglam0<-loglam0c
        lam0<-exp(loglam0)
        llvector<-llvector.new
      }
      #### done updating loglam0 and sigma

      ### Sept 26 2009 added block of code below
      ### to deal with model Mb.
      if(Mb == 1) {   # Only update beta for behavioral response models
        betac<- rnorm(1,beta,.05)
        lpc<-  loglam0 + Mb*betac*prevcap - ((1*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5))
        #lik1c<- log( (1/exp(-exp(lpc)))  -1 )
        lik1c<- log( expm1(exp(lpc[y1]))) ## lik1c<- log( expm1(exp(lpc))) ## MM

        lik2c<- clogloginvpart1(lpc)
        llvector.new<- lik2c
        llvector.new[y1]<- llvector.new[y1]+ lik1c  ## llvector.new[y1]<- llvector.new[y1]+ lik1c[y1] ## MM
        if(runif(1)<exp(sum( rowsum(llvector.new-llvector,indid)[z==1])))
        {
           beta<- betac
           llvector<-llvector.new
        }
      }
      ########################
      # PART 2 OF THE MCMC ALGORITHM UPDATES THE DATA AUGMENTATION PARAMETERS
      # THIS INCLUDES THE LATENT "z" VARIABLES AS WELL AS THE CRITICAL
      # PARAMETER "psi"
      ########################
      # This is the data augmentation part. This code updates each z[i]
      # for all i=1,2,...M individuals. z[i]=1 is a "real" animal that lives
      # in S, whereas z[i]=0 are excess zeros in the data set
      # this is an application of Bayes rule to get Pr(z=1| y[i,,]=0)

      probz <- exp(rowsum(llvector[indid > nind], indid[indid > nind])) # only for nind+1...M
      probz <- (probz*psi) / (probz*psi + (1-psi)) 
      ## [AMG] This can be thought of as follows: 
      ## psi is probability of ## being a real individual 
      ## (detected+undetected) in the inflated population M. probz is
      ## probability of being an ## undetected but real individual. So,
      ## probz*psi is the probability of an individual being real and undetected. 
      z[(nind+1):M] <- rbinom(M-nind, 1, probz)
      psi <- rbeta(1, 1+sum(z), 1+M-sum(z))

      ########################
      ## PART III OF THE ALGORITHM -- THIS BLOCK OF CODE UPDATES THE
      ## ACTIVITY CENTERS
      ###
      # in practice only have to do this for guys with z = 1
      # from the data augmentation.......
      # guys with z = 0 should be drawn uniformly
      # and their acceptance probability should be set to 1.0

      if(bsigma == 1) {       # Only update centers for spatial model
        newcenters<-ceiling(runif(M,0,numnn[centers]))
        newcenters<- NN[cbind(centers,newcenters)]

        # these probabilities are needed in the Metropolis "acceptance probability" calculation
        qnew <- 1/numnn[centers]
        qold <- 1/numnn[newcenters]

        Sc <- G[newcenters,]
        c1c <- (Sc[indid,1]-trapgridbig[,1])^2
        c2c <- (Sc[indid,2]-trapgridbig[,2])^2

        # in theory xx can be 0 which evaluates log(xx/(1-xx)) to -Inf
        # never a problem in practice
        xx <- 1-exp(-(exp(loglam0+ Mb*beta*prevcap[y1]) )* exp(-((1*bsigma)/(sigma^2))*((c1c[y1]+c2c[y1])^(dexp*0.5))))
        xx <- log( xx/(1-xx)  )
        zz <-  -exp(loglam0+Mb*beta*prevcap)*exp(-((1*bsigma)/(sigma^2))*((c1c+c2c)^(dexp*0.5)))
        llvector.tmp <- zz
        llvector.tmp[y1] <- llvector.tmp[y1] + xx
        likdiff <- rowsum(llvector.tmp-llvector, indid)
        likdiff[z==0] <- 0   # this line was in wrong place if using local proposal as above

        likdiff <- likdiff + log(qold/qnew)
        accept <- runif(M)<exp(likdiff)
        S[accept, ] <- Sc[accept, ]
        centers[accept] <- newcenters[accept]
        c1 <- (S[indid,1]-trapgridbig[,1])^2
        c2 <- (S[indid,2]-trapgridbig[,2])^2
      }

      ### should update llvector right here.................
      ## This is taken care of at the beginning of the loop.

      ## Do GOF stuff and save output
      if( (i > burn) & (i%%skip == 0) ) {
        zout[m,] <- z
        Sout[m,] <- centers  # Sout is an index into G, NOT the actual coordinates
        if(bsigma == 1)  {
          out[m,] <- c(sigma*Range/sqrt(dexp), lam0, beta, psi, sum(z))
        } else {
          out[m,] <- c(NA, lam0, beta, psi, sum(z))
        }
        cat(sprintf("%8d %8.3f %5.4f %5.4f %5.4f %6.0f\n", i, out[m,1], out[m,2], out[m,3], out[m,4],out[m,5]))

########## GUI related START ############
        if(Mb == 0)  {   #Rashmi
          statusText <<- paste("sigma\tlam0\tpsi\tN\n",round(out[m,1],2), "\t",round(out[m,2],2), "\t",round(out[m,4],2), "\t",round(out[m,5],2),"\n")
        } else {
          statusText <<- paste("sigma\tlam0\tbeta\tpsi\tN\n",round(out[m,1],2), "\t",round(out[m,2],2),"\t", round(out[m,3],2), "\t",round(out[m,4],2),"\t",round(out[m,5],2),"\n")
        }
        tkinsert(statusWin, "end", statusText)
########## GUI related END #############

        ## Jan 25 GoF stuff for JWM paper ## Moved inside reporting section by MM
        ## 09/19/2011
        #####
        ## Performance issue (GOF is slow, mostly because aggregate is slow):
        ## we're calculating stats for phantoms, then throwing them away.
        ## Changed MM 31 May 2014

        logmu <-loglam0 + Mb*beta*prevcap - ((1*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5))
        # mu <- ( 1-exp(-exp(logmu)))*z[indid]  # zeros out the z=0 guys so they contribute nothing
        # newy <- rbinom(length(mu),1,mu)
        # gof.stats <- cbind(y, newy, mu)
        # gof.stats <- aggregate(gof.stats, list(indid), sum)
        # gof.data[m] <- sum((sqrt(gof.stats[,2])-sqrt(gof.stats[,4]))[z==1]^2)
        # gof.new[m] <- sum((sqrt(gof.stats[,3])-sqrt(gof.stats[,4]))[z==1]^2)

        real <- z[indid] == 1  # TRUE/FALSE vector
        mu <- ( 1 - exp( -exp(logmu[real])))
        newy <- rbinom(sum(real), 1, mu)
        gof.stats <- cbind(y[real], newy, mu)
        gof.stats <- sqrt(rowsum(gof.stats, indid[real], reorder=FALSE)) # rowsum faster than aggregate
        gof.data[m] <- sum((gof.stats[,1] - gof.stats[,3])^2)
        gof.new[m]  <- sum((gof.stats[,2] - gof.stats[,3])^2)

        m <- m+1
      }
    } ############ End of the iterations loop #################################

    end.time <- Sys.time()  ## Add end time and duration to info file
    dur <- end.time - start.time
    #infocon <- file(info.file, "a") ; on.exit(close(infocon))
    sink(info.file, append=TRUE)
      cat(paste("Finished at:", end.time, "\n"))
      cat("Duration:", dur, attr(dur, "units"), "\n\n")
    sink()

    # MM Calculate derived Density, p1 and p2 here;
    #    include p1 & p2 even for non-Mb models
    derivedDen <- (out[,5]/(nG * as.numeric(Global_acSize))) ## *100
    p1 <- 1-exp(-1*out[,2]) # cloglog(p1) = loglam0
    # p2 <- 1-exp(-1*out[,3]) # cloglog(p2) = loglam0 + beta ## this is wrong ## MM
    p2 <- 1 - exp(-exp(log(out[,2])+out[,3]))
    newout <- cbind(out, density=derivedDen, p1=p1, p2=p2)

    
    resSummary <- cbind(colMeans(newout), apply(newout, 2, sd),
      t(apply(newout, 2, hdi)))
    colnames(resSummary) <- c("Posterior_Mean", "Posterior_SD", "95%_Lower_HPD_Level", "95%_Upper_HPD_Level")

########## Output related START ###########
      fname <- paste(folderName, "/param_val_", ts,".csv", sep="")
      write.csv(file=fname, newout)
      cat("Analysis Complete. Parameter estimates written to ", getwd(), "/", fname, sep="", fill=TRUE)
      fname <- paste(folderName, "/summaryStats_", ts,".csv", sep="")
      write.csv(file=fname, resSummary)

      sink(file=info.file, append=TRUE, split=TRUE)  ## Summary added to info file ## MM
      cat("\nSummary of results:\n")
      print(resSummary, digits=3)
      sink()


      ### Addition on October 21, 2011 to generate pixel densities ###
      if(bsigma == 1) {   # only for spatial model
        nSaved <- nrow(Sout) # Sout contains pixel IDs, not coordinates

        # Set pixel ID of home range centers for phantom animals to zero
        indlocs <- Sout * zout
        # Count the proportion of times each pixel was a home range centre,
        #   convert to animals per sq km
        densVec <- tabulate(indlocs, nbins=nG) / nSaved / as.numeric(Global_acSize)
        # Combine densVec with the original habitat data
        #  (unsuitable habitat gets density zero)
        pixelDensity <- grid2500
        pixelDensity$`Pixel Density` <- grid2500[, 3]
        pixelDensity$`Pixel Density`[grid2500[, 3] > 0] <- densVec

        #Generate csv file for pixel densities
        nameoffile3 = paste(folderName,"/pixeldensities_val_",ts,".csv", sep="")
        write.csv(pixelDensity, file =nameoffile3)
      }

############ GUI related START #############
      statusText <<- paste("Analysis Complete\nParameter estimates written to ", getwd(), "/", fname, "\n", sep="")
      tkinsert(statusWin, "end", statusText)
      close(pb)
############ GUI related END ################
############ Output related START ######
      resTable = matrix(NA, nrow=8, ncol=4)  # resTable appears in the GUI
      resTable[,1] = round(apply(newout, 2, mean), digits=4)
      resTable[,2] = round(apply(newout, 2, sd), digits=4)
      for (i in 1:8){
        resTable[i,3:4] <- round(hdi(newout[,i]), digits=4)
      }
      # naming of first col and first row changed to avoid confusion with row/colnames functions # MM
      param.names <- c("_","sigma", "lam0", "beta", "psi", "N", "Density", "p1", "p2")
      stat.names <- c("Posterior_Mean", "Posterior_SD", "95%_Lower_HPD_Level", "95%_Upper_HPD_Level")
      resTable <- cbind(param.names, rbind(stat.names, resTable))

      ### Calculating MCMC diagnostics ###
      showDiag <- rep(TRUE, ncol(out))
      if(bsigma == 0) showDiag[1] <- FALSE  # skip sigma
      if(Mb == 0)     showDiag[3] <- FALSE  # skip beta
      outMCMC <- as.mcmc(out[, showDiag])

      gewekeDiag <- tryCatch(geweke.diag(outMCMC), error=function(e) NULL)
      sink(file=info.file, append=TRUE, split=TRUE)  ## Geweke stat added to info file ## MM
      cat("\nResults of the Geweke Diagnostic:\n")
      print(gewekeDiag)
      sink()

      effective.n <- tryCatch( effectiveSize(outMCMC) , error=function(e) NULL)   ## MM suggests add effectiveSize calculation.... requires coda...
      sink(file=info.file, append=TRUE, split=TRUE)  ## andy  added effective.n
      cat("Effective posterior sample size: \n")
      print(effective.n)
      sink()


########### GUI START ################
    {
      statusText <<- paste("Geweke diagnostic for assessing MCMC convergence is computed and results stored in:", getwd(), "/", folderName, "\n", sep="")
    }
    tkinsert(statusWin, "end", statusText)
########### GUI END ##################

    ### Calculation of Bayesian P-value ###
    sink(file=info.file, append=TRUE, split=TRUE)  ## add GOF stats to info file ## MM
    cat("\nBayesian p-value based on individual encounters: ",
      mean(gof.data > gof.new), "\n\n")
    sink()

########### GUI START ################
      {
      statusText <<- paste("Bayesian P-value for assessing model fit is computed and results stored in:", getwd(), "/", folderName, "\n", sep="")
      }
    tkinsert(statusWin, "end", statusText)
########### GUI END ##################

    ### Display of prior and posterior densities as graphs ###
    lnames <- c('Prior distribution', 'Posterior distribution')
    for ( i in 1:5 ) {
      if(i == 1 && bsigma == 0) next  # skip the plot of sigma for nonspatial model
      if(i == 3 && Mb == 0) next      # skip the plot of beta if no behavioral response
      nameoffile = paste("./", folderName, "/density_", param.names[i+1], "_", ts, ".jpeg", sep="")
      jpeg(filename=nameoffile)
      xlim <- NULL
      if(i == 4)  # psi
        xlim <- c(0, 1)
      if(i == 5)  # N
        xlim <- c(0, nind+nz)
      plot(density(out[,i]), main=param.names[i+1], col="blue", xlim=xlim)
      if(i == 4)  # psi
        segments(c(0, 0, 1), c(0, 1, 1),
                  c(0, 1, 1), c(1, 1, 0), col='red', lty=3)
      if(i == 5)  # N
        segments(c(0, 0, nind+nz), c(0, rep(1/(nind+nz), 2)),
                  c(0,nind+nz, nind+nz), c(rep(1/(nind+nz), 2), 0),
                  col='red', lty=3)
      if(i > 3)
        legend('topleft', lnames, col = c("red","blue"), lty = c(3, 1), bty='n')
      dev.off()
    }

    ### Do graph of detection function
    if(bsigma == 1) {     # only for spatial model
      nameoffile = paste0("./", folderName, "/detectionFunction_", ts, ".jpeg")
      jpeg(filename=nameoffile)
      plotDetFunc(resSummary[2, 1], resSummary[1,1], as.matrix(grid2500[, 1:2]),
        as.matrix(ctLocs), detFunc=c("NE", "HN")[dexp])
      dev.off()
    }

############ Output related END ##############
############ GUI related START ###############
      if(Mb == 0)   #Rashmi
      {
      statusText <<- paste("Density plots of parameters sigma, lam0, psi, N saved in jpg format to ", getwd(), "/", folderName, "\n", sep="")
      }
      else
      {
      statusText <<- paste("Density plots of parameters sigma, lam0, beta, psi, N saved in jpg format to ", getwd(), "/", folderName, "\n", sep="")
      }

      tkinsert(statusWin, "end", statusText)

      tclarray <- tclArray()
      if(tclvalue(rbValue78)=="1")  {
        nrowsTclarray <- 8; nrowsTable <- 9
      } else {     #Rashmi
        nrowsTclarray <- 6; nrowsTable <- 7
      }

      statusText <<- paste("Summary statistics written to ", getwd(), "/", fname, "\n", sep="")
      tkinsert(statusWin, "end", statusText)

      if (Mb == 0)            #Rashmi
      {
      for (i in (0:2))
        for (j in (0:4))
           tclarray[[i,j]] <- resTable[i+1,j+1]

      for (i in (4:6))    #should be 4:6
        for (j in (0:4))
           tclarray[[i-1,j]] <- resTable[i+1,j+1]      #i is made to be i-1
        table1 <- tkwidget(two,"table",variable=tclarray,rows=nrowsTable-1,cols=5,titlerows=1,titlecols=1,selectmode="extended",colwidth=20,background="white",state="disabled",borderwidth=2)
      }
      else
      {
      for (i in (0:nrowsTclarray))
        for (j in (0:4))
           tclarray[[i,j]] <- resTable[i+1,j+1]
           table1 <- tkwidget(two,"table",variable=tclarray,rows=nrowsTable,cols=5,titlerows=1,titlecols=1,selectmode="extended",colwidth=20,background="white",state="disabled",borderwidth=2)
      }

      mod1<-""
      mod2<-""
      if(tclvalue(rbValue34)=="1") {
        mod2<-"Spatial Capture-Recapture"
      } else { mod2<-"Non-Spatial Capture-Recapture" }

        if(tclvalue(rbValue78)=="1")
        { mod1<-"Trap response present"
        }else{ mod1<-"Trap response absent" }

    if(tclvalue(rbValue56)=="1")
      { mod3<-"Half-normal detection function"
      }else{ mod3<-"Negative exponential function" }

        mod4<-"Bernoulli detection process"

        txt1 <- paste("Input Summary\nArea of each pixel representing a potential home-range center:", tclvalue(acSize), "sq km")
        txt2 <- paste("Model selected: ", mod1, ", ", mod2, ", ", mod3, ", ", mod4, sep="")
        txt3 <- paste("MCMC simulation settings: Iterations -", Global_iterNum, "Burnin -", Global_burn, "Thinning -", Global_skip, "Data Augmentation -", Global_nz)
        lab1 <- tklabel(two, text=txt1, background="white")
        lab2 <- tklabel(two, text=txt2, background="white")
        lab3 <- tklabel(two, text=txt3, background="white")
        tkgrid(lab1, sticky="w")
        tkgrid(lab2, sticky="w")
        tkgrid(lab3, sticky="w")

      tkgrid(table1, sticky="w")
      #tkpack(img,fill="x",pady=4)
      #tkpack(s,fill="x",pady=4)
      tkpack(two,fill="both",expand="yes")
########### GUI related END ##################

    }  ## End of SCRd.fn  #################################################


    run <-function()
    {
      if( tclvalue(tkcget(ok1,"-state"))=="normal" | tclvalue(tkcget(ok2,"-state"))=="normal" | tclvalue(tkcget(ok3,"-state"))=="normal") {
          tkmessageBox(message="Please complete input data, model definition and MCMC simulation settings before starting the analysis!",icon="error",type="ok")
          statusText <<- "Please complete input data, model definition and MCMC simulation settings before starting the analysis\n"
          tkinsert(statusWin, "end", statusText)
          tkfocus(tt)
        }else{
          tkdestroy(two)
          two <<- tkframe(tt, width = 800, height = 300, bg= "white",relief="groove",borderwidth=2)
          tkpack(two, side = "top")

          ni <-as.integer(Global_iterNum)
          burn <- as.integer(Global_burn)
          skip <- as.integer(Global_skip)
          nz <- as.integer(Global_nz)

          mod1<-""
          mod2<-""
      mod3<-""
          if(tclvalue(rbValue34)=="1")
          { bsigma=1; mod2<-"Spatial Capture-Recapture"
          }else{ bsigma=0; mod2<-"Non-Spatial Capture-Recapture" }

          if(tclvalue(rbValue78)=="1")
          {Mb=1; mod1<-"Trap response present"
          }else{ Mb=0; mod1<-"Trap response absent" }

      if(tclvalue(rbValue56)=="1")
          {dexp=2; mod3<-"Half-normal detection function"
          }else{ dexp=1; mod3<-"Negative exponential detection function" }

          mod4<-"Bernoulli detection process"

          statusText <<- paste("Input Summary\n------------\nArea of each pixel representing a potential home-range center:", tclvalue(acSize), "sq km")
          statusText <<- paste(statusText, "\nModel selected: ", mod1, ", ", mod2, ", ", mod3, ", ", mod4, sep="")
          statusText <<- paste(statusText, "\nMCMC simulation settings: Iterations -", ni, "Burnin -", burn, "Thinning -", skip, "Data Augmentation -", nz)
          tkinsert(statusWin, "end", statusText)

          SCRd.fn(ni, burn, skip, bsigma, Mb, nz, dexp)  }
    }

    stopit <- function()      {
      shouldIStop <<- TRUE    }

    exit <- function()       {
      tkdestroy(tt)          }

    helpme <- function()
    {
      message <- "Sorry, I can't open the help file. Please Exit the GUI and use help(\"SPACECAP-package\") in the R Console."
      filemessage <- "Sorry, help file does not exist. Please Exit the GUI and use help(\"SPACECAP-package\") in the R Console."
      # tkmessageBox(message=message, icon="error",type="ok") ##### Commented temporarily 12:29 Friday 13 June 

      # help.start() # that works fine, but doesn't get what we want
      # help("SPACECAP") # does nothing
      # help("SPACECAP", help_type="html") # does nothing
      # do.call("help", args=list("SPACECAP-package"), envir=.GlobalEnv) # Does nothing

	r_helpfile_location <- paste(.Library, "/SPACECAP/doc/SPACECAP110_Manual.pdf", sep="")
	if(file.exists(r_helpfile_location)=="TRUE") {
    	if(.Platform$OS.type=="windows"){
	# if(!is.null(try(
	  		tryCatch(shell.exec(r_helpfile_location), warning=function(w) tkmessageBox(message=message,icon="error",type="ok") ,error=function(e) tkmessageBox(message=message,icon="error",type="ok"))
      		}
		else {
	# if(!is.null(try(
	     	tryCatch(system(paste("open", r_helpfile_location)), warning=function(w) tkmessageBox(message=message,icon="error",type="ok"), error=function(e) tkmessageBox(message=message,icon="error",type="ok"))
	     	}
	  	}
	  else {tkmessageBox(message=filemessage, icon="error", type="ok")}   
    }

    tkwm.title(tt,"SPACECAP Ver 1.1.0")

    topMenu <- tkmenu(tt)
    tkconfigure(tt, menu=topMenu)
    fileMenu <- tkmenu(topMenu, tearoff=FALSE)

    tkadd(topMenu, "command", label="Run", command=run)
    tkadd(topMenu, "command", label="Stop", command=stopit)
    tkadd(topMenu, "command", label="Exit", command=exit)
    tkadd(topMenu, "command", label="Help", command=helpme)

    # fileFrame contents
    head1 <- tklabel(fileFrame, text="Input Data", background="light blue", pady=5)
    tkgrid(head1)

    lab1 <- tklabel(fileFrame, text="Select animal capture details file", background="light blue")
    capt <-tkentry(fileFrame, background="white", width=40, state="disable")
    button1 <- tkbutton(fileFrame, text="Browse", command=readFile1, width=10)

    lab2 <- tklabel(fileFrame, text="Select trap deployment details file", background="light blue")
    trap <-tkentry(fileFrame, background="white", width=40, state="disable")
    button2 <- tkbutton(fileFrame, text="Browse", command=readFile2, width=10)

    lab3 <- tklabel(fileFrame, text="Select potential home-range centers data file", background="light blue")
    grids <-tkentry(fileFrame, background="white", width=40, state="disable")
    button3 <- tkbutton(fileFrame, text="Browse", command=readFile3, width=10)

    tkgrid(lab3, sticky="w", padx=5)
    tkgrid(grids, button3, sticky="w", padx=5)
    tkgrid(lab2, sticky="w", padx=5)
    tkgrid(trap, button2, sticky="w", padx=5)
    tkgrid(lab1, sticky="w", padx=5)
    tkgrid(capt, button1, sticky="w", padx=5)

    ac_label1 <- tklabel(fileFrame, text="Specify the area of each pixel (in sq. km)", background="light blue")
    ac_label2 <- tklabel(fileFrame, text="that represents a potential home-range center", background="light blue")
    acSize <- Global_acSize
    ac_size <-tkentry(fileFrame, textvariable = acSize, background="white", width=10)
    tkgrid(ac_label1, sticky="w", padx=5)
    tkgrid(ac_label2, sticky="w", padx=5)
    tkgrid(ac_size, sticky="w", padx=5)

    OnOk1 <- function()
    {
        ## Error checks for input files incorporated here
        ## 1. data is numeric with no NA values
        ## 2. Number of cols and column order is correct
        ## 3. captures cross checked with traps
        errorFlag = 0
        if(is.null(locidso) | is.null(locso) | is.null(grid2500)) {
          print("Error - Input files not selected")
          statusText <<- "Error - Input files not selected\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error - Input files not selected!",icon="error",type="ok")
          errorFlag = -1 }

        ##Check locidso
        if(dim(locidso)[1]==0)
        {
          print("Error in animal capture details file - No data or bad file format, csv file expected")
          statusText <<- "Error in animal capture details file - No data or bad file format, csv file expected\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in animal capture details file - No data or bad file format, csv file expected!",icon="error",type="ok")
          errorFlag = -1
        }
        if(dim(locidso)[2]==3)
        {
          if (is.integer(locidso[,1]) & is.integer(locidso[,2]) & is.integer(locidso[,3]) & sum(is.na(locidso))==0)
          {
            if(toupper(names(locidso)[1])!="LOC_ID" | toupper(names(locidso)[2])!="ANIMAL_ID" | toupper(names(locidso)[3])!="SO") {
               print("Error in animal capture details file - incorrect column sequence or header")
               statusText <<- "Error in animal capture details file - incorrect column sequence or header\n"
               tkinsert(statusWin, "end", statusText)
               tkmessageBox(message="Error in animal capture details file - incorrect column sequence or header!",icon="error",type="ok")
               errorFlag = -1  }
          }else {
            print("Error in animal capture details file - non-integer or missing values")
            statusText <<- "Error in animal capture details file - non-integer or missing values\n"
            tkinsert(statusWin, "end", statusText)
            tkmessageBox(message="Error in animal capture details file - non-integer or missing values!",icon="error",type="ok")
            errorFlag = -1 }
        }else {
          print("Error in animal capture details file - incorrect number of columns")
          statusText <<- "Error in animal capture details file - incorrect number of columns\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in animal capture details file - incorrect number of columns!",icon="error",type="ok")
          errorFlag = -1 }
        ### Check locso
        if(dim(locso)[1]==0)
        {
          print("Error in trap deployment details file - No data or bad file format, csv file expected")
          statusText <<- "Error in trap deployment details file - No data or bad file format, csv file expected\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in trap deployment details file - No data or bad file format, csv file expected!",icon="error",type="ok")
          errorFlag = -1
        }
        if(is.integer(locso[,1]) & is.numeric(locso[,2]) & is.numeric(locso[,3]) & sum(is.na(locso))==0)
        {
          if(toupper(names(locso)[1])!="LOC_ID" | toupper(names(locso)[2])!="X_COORD" | toupper(names(locso)[3])!="Y_COORD") {
            print("Error in trap deployment details file - incorrect column sequence or header")
            statusText <<- "Error in trap deployment details file - incorrect column sequence or header\n"
            tkinsert(statusWin, "end", statusText)
            tkmessageBox(message="Error in trap deployment details file - incorrect column sequence or header!",icon="error",type="ok")
            errorFlag = -1 }
        }else {
          print("Error in trap deployment details file - non-numeric or missing values")
          statusText <<- "Error in trap deployment details file - non-numeric or missing values\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in trap deployment details file - non-numeric or missing values!",icon="error",type="ok")
          errorFlag = -1 }
        nso = (dim(locso)[2]) - 3
        nloc = dim(locso)[1]
        for(i in 4:(nso+3)) {
          if (sum(locso[,i]>=0)<nloc | sum(locso[,i]<=1)<nloc) {
            statusText <<- paste("Error in trap deployment details file - incorrect trap status in column ", i,", should be 0 or 1\n", sep="")
            print(statusText)
            tkinsert(statusWin, "end", statusText)
            tkmessageBox(message=statusText,icon="error",type="ok")
            errorFlag = -1 } }
        ### Check grids
        if(dim(grid2500)[1]==0)
        {
          print("Error in potential home range center data file - No data or bad file format, csv file expected")
          statusText <<- "Error in potential home range center data file - No data or bad file format, csv file expected\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in potential home range center data file - No data or bad file format, csv file expected!",icon="error",type="ok")
          errorFlag = -1
        }
        if(dim(grid2500)[2]==3)
        {
          if(is.numeric(grid2500[,1]) & is.numeric(grid2500[,2]) & is.integer(grid2500[,3]) & sum(is.na(grid2500))==0)
          {
            if(toupper(names(grid2500)[1])!="X_COORD" | toupper(names(grid2500)[2])!="Y_COORD" | toupper(names(grid2500)[3])!="HABITAT")  {
               print("Error in potential home range center data file - incorrect column sequence or header")
               statusText <<- "Error in potential home range center data file - incorrect column sequence or header\n"
               tkinsert(statusWin, "end", statusText)
               tkmessageBox(message="Error in potential home range center data file - incorrect column sequence or header!",icon="error",type="ok")
               errorFlag = -1  }
          }else {
            print("Error in potential home range center data file - non-numeric or missing values")
            statusText <<- "Error in potential home range center data file - non-numeric or missing values\n"
            tkinsert(statusWin, "end", statusText)
            tkmessageBox(message="Error in potential home range center data file - non-numeric or missing values!",icon="error",type="ok")
            errorFlag = -1   }
        }else  {
          print("Error in potential home range center data file - incorrect number of columns")
          statusText <<- "Error in potential home range center data file - incorrect number of columns\n"
          tkinsert(statusWin, "end", statusText)
          tkmessageBox(message="Error in potential home range center data file - incorrect number of columns!",icon="error",type="ok")
          errorFlag = -1   }

        #### Check captures vs traps
        len = length(locidso[,1])
        for (i in 1:len)
        {
          loc = locidso[i,1]
          so = locidso[i,3]

          if(locso[loc,so+3]==0) {
            cat("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so, fill=TRUE)
            statusText <<- paste("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so)
            tkinsert(statusWin, "end", statusText)
            tkmessageBox(message=paste("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so),icon="error",type="ok")
            errorFlag = -1   }
        }

        Global_acSize <<- as.numeric(tclvalue(acSize))

        if (is.na(Global_acSize)) {
          tkmessageBox(message="Error - potential home-range center area should be numeric!",icon="error",type="ok")
          tkfocus(ac_size)
          statusText <<- "Error - potential home-range center area should be numeric\n"
          tkinsert(statusWin, "end", statusText)
          errorFlag = -1   }

      ##########
      if(errorFlag==-1) {
        statusText <<- "Input data could not be read\n"
        tkinsert(statusWin, "end", statusText)
      }else {
        statusText <<- "Input data read successfully\n"
        tkinsert(statusWin, "end", statusText)
        tkconfigure(button1,state="disable")
        tkconfigure(button2,state="disable")
        tkconfigure(button3,state="disable")
        tkconfigure(ac_size,state="disable")
        tkconfigure(ok1,state="disable")  }

        #tkentryconfigure(topMenu,1,state="active")
    }

    OnReset1 <- function()
    {
      tkconfigure(button1,state="active")
      tkconfigure(button2,state="active")
      tkconfigure(button3,state="active")
      tkconfigure(ac_size,state="normal")
      tkconfigure(ok1,state="active")
    }

    ok1 <-tkbutton(fileFrame,text="   OK   ",command = OnOk1, width=8)
    reset1 <-tkbutton(fileFrame,text="  Edit  ",command = OnReset1, width=8)
    tkgrid(ok1, reset1, padx=5,pady=40)
    tkgrid.configure(ok1, sticky="e")
    tkgrid.configure(reset1, sticky="w")

    #choiceFrame contents
    head2 <- tklabel(choiceFrame, text="Model Definition", background="light blue", pady=5)
    tkgrid(head2)

    rb78_label <- tklabel(choiceFrame,text="Trap response", background = "light blue", padx=10, pady=1)
    rb7 <- tkradiobutton(choiceFrame, text = "Trap response present", background = "light blue", padx=10, pady=1)
    rb8 <- tkradiobutton(choiceFrame, text = "Trap response absent", background = "light blue", padx=10, pady=1)
    rbValue78 <- tclVar("0")
    tkconfigure(rb7,variable=rbValue78,value="1")
    tkconfigure(rb8,variable=rbValue78,value="0")
    tkgrid( rb78_label,sticky="w" )
    tkgrid( rb7,sticky="w" )
    tkgrid( rb8,sticky="w" )

    rb34_label <- tklabel(choiceFrame,text="Capture-Recapture model", background = "light blue", padx=10, pady=1)
    rb3 <- tkradiobutton(choiceFrame, text = "Spatial Capture-Recapture", background = "light blue", padx=10, pady=1)
    rb4 <- tkradiobutton(choiceFrame, text = "Non-Spatial Capture-Recapture", background = "light blue", padx=10, pady=1)
    rbValue34 <- tclVar("1")
    tkconfigure(rb3,variable=rbValue34,value="1")
    tkconfigure(rb4,variable=rbValue34,value="0")
    tkgrid( rb34_label,sticky="w" )
    tkgrid( rb3,sticky="w" )
    tkgrid( rb4,sticky="w" )

    rb56_label <- tklabel(choiceFrame,text="Detection function", background = "light blue", padx=10, pady=1)
    rb5 <- tkradiobutton(choiceFrame, text = "Half Normal", background = "light blue",padx=10, pady=1)
    rb6 <- tkradiobutton(choiceFrame, text = "Negative Exponential", background = "light blue",padx=10, pady=1)
    rbValue56 <- tclVar("1")
    tkconfigure(rb5,variable=rbValue56,value="1")
    tkconfigure(rb6,variable=rbValue56,value="0")
    tkgrid( rb56_label,sticky="w" )
    tkgrid( rb5,sticky="w")
    tkgrid( rb6,sticky="w")

    rb12_label <- tklabel(choiceFrame,text="Capture encounters", background = "light blue", padx=10, pady=1)
    rb1 <- tkradiobutton(choiceFrame, text = "Bernoulli process", background = "light blue", padx=10, pady=1)
    rb2 <- tkradiobutton(choiceFrame, text = "Poisson process", background = "light blue", padx=10, pady=1)
    rbValue12 <- tclVar("1")
    tkconfigure(rb1,variable=rbValue12,value="1",state="disabled")
    tkconfigure(rb2,variable=rbValue12,value="0",state="disabled")
    tkgrid( rb12_label,sticky="w" )
    tkgrid( rb1,sticky="w" )
    tkgrid( rb2,sticky="w")

    onOk2 <- function()
    {
      tkconfigure(rb3,state="disable")
      tkconfigure(rb4,state="disable")
    tkconfigure(rb5,state="disable")
    tkconfigure(rb6,state="disable")
      tkconfigure(rb7,state="disable")
      tkconfigure(rb8,state="disable")
      tkconfigure(ok2,state="disable")
      statusText <<- "Model definition complete\n"
      tkinsert(statusWin, "end", statusText)
    }
    onReset2 <- function()
    {
      tkconfigure(rb3,state="normal")
      tkconfigure(rb4,state="normal")
    tkconfigure(rb5,state="normal")
    tkconfigure(rb6,state="normal")
      tkconfigure(rb7,state="normal")
      tkconfigure(rb8,state="normal")
      tkconfigure(ok2,state="active")
    }
    ok2 <-tkbutton(choiceFrame,text="     OK      ",command=onOk2, width=8)
    reset2 <-tkbutton(choiceFrame,text="  Edit  ",command=onReset2, width=8)
    tkgrid(ok2, reset2, padx=5,pady=10)
    #tkgrid(reset2, padx=10,pady=10, sticky="w")
    #tkgrid(ok1, reset1, padx=5,pady=40)
    tkgrid.configure(ok2, sticky="e")
    tkgrid.configure(reset2, sticky="w")


    #valueFrame
    head3 <- tklabel(valueFrame, text="MCMC simulations settings", background="light blue", pady=5)
    tkgrid(head3)

    rb_7_label1 <- tklabel(valueFrame, text="Specify number of MCMC iterations", background="light blue", padx=5, pady=1)
    rb_7_label2 <- tklabel(valueFrame, text="(usually a large number [>50000])", background="light blue", padx=5, pady=1)
    ItrNum <- Global_iterNum
    Iterations <-tkentry(valueFrame, textvariable = ItrNum, background="white", width=10)
    tkgrid( rb_7_label1,sticky="w")
    tkgrid( rb_7_label2,sticky="w")
    tkgrid(Iterations,sticky="w", padx=5)

    rb_8_label1 <- tklabel(valueFrame, text="Specify the burn-in period", background="light blue", padx=5, pady=1)
    rb_8_label2 <- tklabel(valueFrame, text="(no. of initial iterations to be discarded)", background="light blue", padx=5, pady=1)
    burnIn <- Global_burn
    BurnIn <-tkentry(valueFrame, textvariable = burnIn, background="white", width=10)
    tkgrid(rb_8_label1,sticky="w")
    tkgrid(rb_8_label2,sticky="w")
    tkgrid( BurnIn,sticky="w", padx=5)


    rb_9_label1 <- tklabel(valueFrame, text="Specify thinning rate", background="light blue", padx=5, pady=1)
    rb_9_label2 <- tklabel(valueFrame, text="(no. of iterations to be skipped when", background="light blue", padx=5, pady=1)
    rb_9_label3 <- tklabel(valueFrame, text="reporting summary statistics)", background="light blue", padx=5, pady=1)
    skip <- Global_skip
    Skip <-tkentry(valueFrame, textvariable = skip, background="white", width=10)
    tkgrid( rb_9_label1,sticky="w")
    tkgrid( rb_9_label2,sticky="w")
    tkgrid( rb_9_label3,sticky="w")
    tkgrid( Skip,sticky="w" , padx=5)


    rb_10_label1 <- tklabel(valueFrame, text="Specify data augmentation value", background="light blue", padx=5, pady=1)
    rb_10_label2 <- tklabel(valueFrame, text="(an upper limit to the number of individuals", background="light blue", padx=5, pady=1)
    rb_10_label3 <- tklabel(valueFrame, text="in the area)", background="light blue", padx=5, pady=1)
    nz <- Global_nz
    NZ <-tkentry(valueFrame, textvariable = nz, background="white", width=10)
    tkgrid( rb_10_label1,sticky="w")
    tkgrid( rb_10_label2,sticky="w")
    tkgrid( rb_10_label3,sticky="w")
    tkgrid(NZ,sticky="w", padx=5)
    #
    tkpack(one,two,three,expand=TRUE, fill="both")
    tkpack(fileFrame,choiceFrame,valueFrame,expand=TRUE, fill="both")

    OnOk3 <- function()
    {
      errorFlag <- 0
      Global_iterNum <<- as.numeric(tclvalue(ItrNum))
      if (is.na(Global_iterNum)) {
          tkmessageBox(message="Number of iterations should be an integer!",icon="error",type="ok")
          tkfocus(Iterations)
          statusText <<- "Error - Number of iterations should be numeric\n"
          tkinsert(statusWin, "end", statusText)
          errorFlag <- -1        }

      Global_burn <<- as.numeric(tclvalue(burnIn))
      if (is.na(Global_burn))   {
         tkmessageBox(message="Burn-in value should be an integer!",icon="error",type="ok")
         tkfocus(BurnIn)
         statusText <<- "Error - Burn-in value should be an integer\n"
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1       }

      Global_skip <<- as.numeric(tclvalue(skip))
      if (is.na(Global_skip))   {
         tkmessageBox(message="Thinning rate should be an integer!",icon="error",type="ok")
         tkfocus(Skip)
         statusText <<- "Error - Thinning rate should be an integer\n"
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1        }

      Global_nz <<- as.numeric(tclvalue(nz))
      if (is.na(Global_nz))     {
         tkmessageBox(message="Data augmentation value should be an integer!",icon="error",type="ok")
         tkfocus(NZ)
         statusText <<- "Error - Data augmentation value should be an integer\n"
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1        }

      if(Global_burn >= Global_iterNum) {
         tkmessageBox(message="Burn-in value should be less than number of iterations!",icon="error",type="ok")
         tkfocus(BurnIn)
         statusText <<- "Error - Burn-in value should be less than number of iterations\n"
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1               }

      if(Global_skip >= Global_iterNum) {
         tkmessageBox(message="Thinning should be less than number of iterations!",icon="error",type="ok")
         tkfocus(BurnIn)
         statusText <<- "Error - Thinning should be less than number of iterations\n"
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1               }

      if(errorFlag==-1) {
        statusText <<- "MCMC simulation settings incorrect\n"
        tkinsert(statusWin, "end", statusText)
      }else {
        tkconfigure(Iterations, state="disabled")
        tkconfigure(BurnIn, state="disabled")
        tkconfigure(Skip, state="disabled")
        tkconfigure(NZ, state="disabled")
        tkconfigure(ok3,state="disabled")
        statusText <<- "MCMC simulation settings complete\n"
        tkinsert(statusWin, "end", statusText)  }
    }

    onReset3 <- function()
    {
      tkconfigure(Iterations, state="normal")
      tkconfigure(BurnIn, state="normal")
      tkconfigure(Skip, state="normal")
      tkconfigure(NZ, state="normal")
      tkconfigure(ok3,state="active")
      #tkentryconfigure(topMenu,1,state="disabled")
    }
    ok3 <- tkbutton(valueFrame,text="     OK      ",command=OnOk3, width=8)
    reset3 <-tkbutton(valueFrame,text="  Edit  ",command=onReset3, width=8)
    tkgrid(ok3, reset3, padx=10,pady=10)
    #tkgrid(reset3, padx=10,pady=10, sticky="w")
    tkgrid.configure(ok3, sticky="e")
    tkgrid.configure(reset3, sticky="w")
############################ GUI related END ################################

}
