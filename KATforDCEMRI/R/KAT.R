##    KATforDCEMRI: a Kinetic Analysis Tool for DCE-MRI
##    Copyright 2014 Genentech, Inc.
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com

KAT <- function(file="concatenate.KAT.with.KAT.checkData.RData", results_file="my_results", method.optimization="L-BFGS-B", show.rt.fits=FALSE, param.for.avdt="Ktrans", range.map=1.5, cutoff.map=0.85, export.matlab=TRUE, export.RData=TRUE, verbose=FALSE, show.errors = FALSE, try.silent=TRUE, fracGTzero=0.75, AIF.shift="", Force.AIF.peak=FALSE, tlag.Tofts.on=FALSE, ...){

    ## SPECIFY LOWER BOUND FOR PARAMETER ESTIMATES
    lo <- 0

    options(show.error.messsages = show.errors)

    ftype <- strsplit(file, split="a")[[1]]
    ftype <- ftype[length(ftype)-1]
    ftype <- strsplit(ftype, split="")[[1]]
    ftype <- ftype[length(ftype)]

    if(ftype =="m")
        file.format <- "matlab"
    if(ftype =="D")
        file.format <- "RData"

    file_short <- file

    KAT.version <- "0.740"
    ptm_total <- proc.time()[3]

    modeltype1 <- "xTofts"
    modeltype2 <- "Tofts"


    ## ########  ########  ######   #### ##    ##
    ## ##     ## ##       ##    ##   ##  ###   ##
    ## ##     ## ##       ##         ##  ####  ##
    ## ########  ######   ##   ####  ##  ## ## ##
    ## ##     ## ##       ##    ##   ##  ##  ####
    ## ##     ## ##       ##    ##   ##  ##   ###
    ## ########  ########  ######   #### ##    ##

    ## ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
    ## ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
    ## ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
    ## ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
    ## ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
    ## ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
    ## ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######


    ## ############################################################
    ## #####FUNCTION TO DEAL WITH FUSSY NUMBERS####################
    ## ############################################################
    funcz <- function(x){
        as.numeric(as.character(x))
    }


    ## ############################################################################
    ## #####SHIFT THE AIF/VIF######################################################
    ## ############################################################################
    aif.shift.func <- function(t, cp, time_shift){
        ## GF add peak AIF value into vector
        x <- cbind(t,cp)
        tmax <- subset(x[,1], x[,2]==max(x[,2]))
        if(AIF.shift=="ARTERY")
            tmax <- tmax + time_shift
        if(AIF.shift=="VEIN")
            tmax <- tmax - time_shift

        ## GF Approximate the AIF as a function
        cpFUNC <- approxfun(t, cp, rule=2)

        ## GF Shift the time vector
        if(AIF.shift=="ARTERY")
            tshift <- t - time_shift
        if(AIF.shift=="VEIN")
            tshift <- t + time_shift

        ## GF Calculate the time shifted AIF
        cp.shift <- cpFUNC(tshift)

        ## GF Replace peak value of cp.shift with peak value of cp
        if(Force.AIF.peak==TRUE)
            cp.shift[cp.shift==max(cp.shift)] <- max(cp)

        return(cp.shift)
    }

    ## ##################################################################################
    ## #####SPECIFY THE ONE COMPARTMENT EXTENDED TOFTS MODEL WITH ESTIMATED TIME DELAY###
    ## ##################################################################################
    roi.modelT <- function(p, t, dt, cp, zing=0) {
        if(zing==0){

            Ktrans <- p[1]
            kep <- p[2]
            ## GF Add a time shift parameter
            ## GF The fix.tlag argument allows tlag to be estimated based on median roi values and fixed to that estimated value for subsequest per-voxel fits
            if(tlag.Tofts.on == TRUE){
                if(fix.tlag != TRUE & AIF.shift != "NONE")
                    time_shift <- p[3]

                if(fix.tlag == TRUE & AIF.shift != "NONE")
                    time_shift <- param.est.whole.roi.Tofts$tlag

            }

            if(tlag.Tofts.on == FALSE)
                time_shift <- param.est.whole.roi.xTofts$tlag

            ## GP ODE: ct'(t) = cp(t)*ktrans - ct(t)*kep ; c(0) = 0
            ## GP closed form solution: ct(t) = exp(-kep*t) * int(ktrans*exp(kep*tau)*cp(tau), tau=1..t)
            ## GP we are here approximating int(f(t)dt) by the the non-uniform trapezoidal rule:
            ## GP int(f(t)dt, t=a..b) = sum((t_{k+1} - t_{k})*(f(x_{k+1})+f(x_{k})), k=1..n)

            if(AIF.shift=="ARTERY" | AIF.shift=="VEIN")
                ## GF Shift the AIF
                cp <- aif.shift.func(t, cp, time_shift)

            f <- Ktrans*exp(kep*t)*cp
            int <- c(0, cumsum(dt*(tail(f, -1)+head(f, -1)))/2)
            ct <- exp(-kep*t)*int

            return(ct)
        }
        if(zing==1)
            return("modelT")
    }



    ## ######################################################################################
    ## #####SPECIFY THE ONE COMPARTMENT EXTENDED TOFTS MODEL WITH ESTIMATED TIME DELAY#######
    ## ######################################################################################
    roi.modelxT <- function(p, t, dt, cp, zing=0) {
        if(zing==0){

            Ktrans <- p[1]
            kep <- p[2]
            vb <- p[3]
            ## GF Add a time shift parameter
            ## GF The fix.tlag argument allows tlag to be estimated based on median roi values and fixed to that estimated value for subsequest per-voxel fits
            if(fix.tlag != TRUE & AIF.shift != "NONE")
                time_shift <- p[4]

            if(fix.tlag == TRUE & AIF.shift != "NONE")
                time_shift <- param.est.whole.roi.xTofts$tlag


            ## GP ODE: ct'(t) = cp(t)*ktrans - ct(t)*kep ; c(0) = 0
            ## GP closed form solution: ct(t) = exp(-kep*t) * int(ktrans*exp(kep*tau)*cp(tau), tau=1..t)
            ## GP we are here approximating int(f(t)dt) by the the non-uniform trapezoidal rule:
            ## GP int(f(t)dt, t=a..b) = sum((t_{k+1} - t_{k})*(f(x_{k+1})+f(x_{k})), k=1..n)


            if(AIF.shift=="ARTERY" | AIF.shift=="VEIN")
                ## GF Shift the AIF
                cp <- aif.shift.func(t, cp, time_shift)

            f <- Ktrans*exp(kep*t)*cp
            int <- c(0, cumsum(dt*(tail(f, -1)+head(f, -1)))/2)
            ct <- exp(-kep*t)*int

            ## GP last transformation
            ct <- ct + vb*cp

            return(ct)
        }
        if(zing==1)
            return("modelxT")
    }

    ## #############################################################
    ## ##BASIC DECONVOLUTION FUNCTION###############################
    ## #############################################################
    calch <- function(u, y, TIME_trunc){
        ## TRY ALTERNATE PREPLOT PARAMETER HERE
        locfit_y <- preplot(y, newdata=0:max(TIME_trunc))
        ## locfit_y <- preplot(y, newdata=TIME_trunc)
        y_smooth <- locfit_y$fit

        ## Truncate pre-peak VIF values
        u <- u[match(u[u==max(u)], u):length(u)]
        y_smooth <- y_smooth[match(u[u==max(u)], u):length(u)]

        ## Generate 'A' matrix
        n<-length(u)
        A<-matrix(0,nrow=n,ncol=n)
        ind<-row(A)-col(A)
        ind[ind<0] <- (-1)
        ind<-ind+2
        A <- matrix(c(0,u)[ind],nrow=n,ncol=n)

        ## Solve for Impulse Response Function
        h <- solve(A,y_smooth)

        h.time.vector <- (1:length(h))-1
        out <- list(h.time.vector, h)
        names(out) <-c("t", "IRF")
        return(out)
    }


    ## #############################################################
    ## ##CALCULATE IRF AND AUC AND AUC/MRT OF IRF###################
    ## #############################################################
    calchFUNC <- function(vector.times, AIF, map_cc_slice, correct.vp = TRUE, alpha.AIF = c(0.1, 0.5), vp.nom=0.1, kep.nom=0.5){
        AUMC <- function(AUMC.median, h.median, irf_time_vec, r){
            AUMC.median <- AUMC.median + 0.5*(h.median[r]*irf_time_vec[r] + h.median[r+1]*irf_time_vec[r+1])
        }

        artery_data <- data.frame(vector.times*60, AIF)
        names(artery_data) <- c("TIME", "ARTERY")
        data_artery_peak <- subset(artery_data, artery_data$ARTERY == max(artery_data$ARTERY))
        data_remove_artery_prepeak <- subset(artery_data, artery_data$TIME >= data_artery_peak$TIME)
        frames_to_peak <- length(artery_data[,1]) - length(data_remove_artery_prepeak[,1]) + 1
        TIME   <- data_remove_artery_prepeak$TIME
        ARTERY <- data_remove_artery_prepeak$ARTERY
        TIME_trunc <- TIME[seq(1,length(TIME)-1,by=1)]
        TIME_trunc <- TIME_trunc - TIME_trunc[1]
        ARTERY_trunc <- ARTERY[seq(1,length(ARTERY)-1,by=1)]

        ARTERY_smooth <- locfit.robust(ARTERY_trunc~TIME_trunc, acri="cp", alpha=alpha.AIF)
        AIF_smooth <- ARTERY_smooth

        locfit_u <- preplot(AIF_smooth, newdata=0:max(TIME_trunc))
        u_smooth <- locfit_u$fit

        Tmax <- max(TIME_trunc)

        ## #########Perform deconvolution analysis on median voxel-wise curves#

        TUMOR.median <- map_cc_slice
        TUMOR.median <- TUMOR.median[seq(frames_to_peak,length(TUMOR.median),by=1)]

        if(vp.nom > 0)
            TUMOR.median_corr <- TUMOR.median - vp.nom*ARTERY

        if(vp.nom <= 0)
            TUMOR.median_corr <- TUMOR.median

        TUMOR.median_corr_shifted <- TUMOR.median_corr[seq(2,length(TUMOR.median_corr),by=1)]
        TUMOR.median_smooth <- locfit.robust(TUMOR.median_corr_shifted~TIME_trunc, acri="cp")

        calch.out <- calch(u_smooth, TUMOR.median_smooth, TIME_trunc)

        h.median <- calch.out$IRF
        irf_time_vec <- calch.out$t
        n <- length(h.median)

        AUC.median <- 0
        AUMC.median <- 0
        for(r in 1:(n-1)){
            h_sum <- h.median[r] + h.median[r+1]
            t_sum <- irf_time_vec[r] + irf_time_vec[r+1]
            AUC.median <- AUC.median + 0.5*h_sum
            AUMC.median <-  AUMC(AUMC.median, h.median, irf_time_vec, r)
        }

        AUCMRT.median <- AUC.median/(AUMC.median/AUC.median)*60

        if(kep.nom > 0){
            t_scan <- max(TIME_trunc)/60
            ve_trunc_error <- 1-exp(-kep.nom*t_scan)
            Ktrans_trunc_error <- (1-exp(-kep.nom*t_scan))^2/(1-(1+kep.nom*t_scan)*exp(-kep.nom*t_scan))
            AUC.median <- AUC.median / ve_trunc_error
            AUCMRT.median <- AUCMRT.median / Ktrans_trunc_error
        }

        irf_time_vec <- irf_time_vec/60
        out <- list(AUC.median, AUCMRT.median, h.median, irf_time_vec)

        names(out) <- c("AUCh", "AUChMRTh", "IRF", "t")
        return(out)
    }


    ## ###########################################
    ## ##SPECIFY THE OBJECTIVE FUNCTION ##########
    ## ###########################################
    Obj_roi <- function(p, model, t, dt, cp, roi) {
        sum((model(p, t, dt, cp) - roi)^2)
    }

    ## ######## ##    ## ########
    ## ##       ###   ## ##     ##
    ## ##       ####  ## ##     ##
    ## ######   ## ## ## ##     ##
    ## ##       ##  #### ##     ##
    ## ##       ##   ### ##     ##
    ## ######## ##    ## ########

    ## ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
    ## ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
    ## ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
    ## ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
    ## ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
    ## ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
    ## ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######


    map_cc_slice <- NULL
    map_cc_roi <- NULL
    aif <- NULL
    aif.shifted <- NULL
    map.times <- NULL
    map.KtransxT <- NULL
    map.kepxT <- NULL
    map.vexT <- NULL
    map.vbxT <- NULL
    map.fitfailuresxT <- NULL
    map.KtransT <- NULL
    map.kepT <- NULL
    map.veT <- NULL
    map.fitfailuresT <- NULL
    mask.roi <- NULL
    nx <- NULL
    ny <- NULL
    nt <- NULL
    map.KtransxT <- NULL
    map.kepxT <- NULL
    map.vexT <- NULL
    map.vbxT <- NULL
    map.AIC.compare <- NULL
    map.AIC.T <- NULL
    map.EF <- NULL
    map.AIC.xT <- NULL
    roi.median.fitted.Tofts <- NULL
    param.est.whole.roi.Tofts <- NULL
    roi.median.fitted.xTofts <- NULL
    param.est.whole.roi.xTofts <- NULL
    cv.whole.roi.xTofts <- NULL

    ## #################################################
    ## #PACKAGE INFORMATION#############################
    ## #################################################
    cat("\n")
    cat("#########################################################################", "\n")
    cat("##-------------- KINETIC ANALYSIS TOOL (KAT) FOR DCE-MRI --------------##", "\n")
    cat("#########################################################################", "\n")
    cat("##---------------------- R package version",KAT.version,"----------------------##", "\n")
    cat("#########################################################################", "\n")
    cat("\n")

    ## #IMPORT MATLAB FILE INTO R#######################
    filea <- strsplit(file, split="/")[[1]]
    fileb <- filea[length(filea)]

    if(file!="concatenate.KAT.with.KAT.checkData.RData")
        cat("loading", fileb, "into R...", "\n")
    ptm <- proc.time()[3]

    if(file.format == "matlab"){
        mat_data <- readMat(file)

        if(length(names(mat_data)) < 10)
            file.original <- TRUE

        if(length(names(mat_data)) >= 10)
            file.original <- FALSE
    }

    if(file.format == "RData"){
        if(file!="concatenate.KAT.with.KAT.checkData.RData")
            load(file)

        if(length(names(dcemri.data)) < 10){
            file.original <- TRUE
            ROIcounter <- apply(dcemri.data$maskROI, 3, max)
            results_file_temp <- results_file
        }

        if(length(names(dcemri.data)) >= 10){
            file.original <- FALSE
            ROIcounter <- 1
        }
    }

    if(file!="concatenate.KAT.with.KAT.checkData.RData"){
        cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
        cat("--------", "\n")
    }

    for(slicenumber in 1:(length(ROIcounter))){

        if(ROIcounter[slicenumber]==1){
            if(file.original == TRUE){
                slice <- slicenumber
                cat("***** ROI DETECTED IN SLICE", slice, " *****", "\n")
                cat("--------", "\n")
            }

            if(file.original == TRUE){

                ## ADD INITIAL VALUE OF TIME SHIFT INTO INNITIAL PARAMS FOR TOFTS AND XTOFTS MODELS

                if(AIF.shift!="VEIN" & AIF.shift!="ARTERY" & AIF.shift!="NONE")
                    stop("You must specify the argument AIF.shift argument as VEIN, ARTERY or NONE, indicating that the AIF you are using is based on data from a vein or artery or NONE if tlag should be set to 0.  This will ensure that the time lag parameter in the Tofts and xTofts models has the appropriate inital value and is bounded correctly; either -Inf to 0 (for VEIN) or 0 to Inf (for ARTERY).")

                ## APPEND SLICE NUMBER TO FILE NAME
                filenameTag <- paste("_slice", slice, sep="")
                results_file <- paste(results_file_temp, filenameTag, sep="")

                roi.model <- roi.modelxT

                if(slice=="" || slice < 0)
                    stop("The slice argument has not been properly specified; slice=``slice number''")

                cat("extracting slice", slice, "for analysis...", "\n")
                ptm <- proc.time()[3]

                if(file.format == "matlab"){
                    ## #EXTRACT THE TIME VECTOR AND CONVERT TO MINUTES########################
                    map.times <- as.vector(mat_data$map[[4]]/60)

                    ## #EXTRACT CONTRAST AGENT CONCENTRATIONS FOR SLICE OF INTEREST###########
                    map_cc <- mat_data$map[[3]]
                    map_cc_slice <- map_cc[,,slice,]

                    if(length(map.times)!=length(map_cc_slice[1,1,]))
                        stop("The length of the time vector does not match the length of the contrast agent concentration vector")
                }

                if(file.format == "RData"){
                    ## #EXTRACT THE TIME VECTOR AND CONVERT TO MINUTES########################
                    map.times <- as.vector(dcemri.data$vectorTimes/60)


                    ## #EXTRACT CONTRAST AGENT CONCENTRATIONS FOR SLICE OF INTEREST###########
                    map_cc <- dcemri.data$mapCC

                    map_cc_slice <- map_cc[,,slice,]

                    if(length(map.times)!=length(map_cc_slice[1,1,]))
                        stop("The length of the time vector does not match the length of the contrast agent concentration vector")
                }

                nt <- length(map_cc_slice[1,1,])
                ny <- length(map_cc_slice[,1,1])
                nx <- length(map_cc_slice[1,,1])

                ccTEMP <- rep(0, prod(dim(map_cc_slice)))
                dim(ccTEMP) <- dim(map_cc_slice)[c(2,1,3)]

                for(i in 1:nt)
                    ccTEMP[,,i] <- rot90(map_cc_slice[,,i],3)
                map_cc_slice <- ccTEMP


                ## #EXTRACT ROI MASK FOR SLICE OF INTEREST################################
                if(file.format == "matlab")
                    mask.roi <- mat_data$mask[[1]][,,slice]

                if(file.format == "RData")
                    mask.roi <- dcemri.data$mask[,,slice]

                mask.roi <- rot90(mask.roi,3)

                if(max(mask.roi)!=1)
                    stop("Your ROI mask is either composed entirely of zeroes or contains nonnumeric elements; voxels within the ROI should have a value of ``1'' and all other voxels should have a value of ``0''.")


                ## #EXTRACT THE AIF VECTOR FROM THE DATA FILE#########################
                if(file.format == "matlab")
                    aif <- as.vector(mat_data$aif)
                if(file.format == "RData")
                    aif <- as.vector(dcemri.data$vectorAIF)

                if(file.format == "matlab")
                    rm(mat_data)

                cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
                cat("--------", "\n")

                ## #########################################
                ## ##APPLY ROI MASKS TO map_cc_slice########
                ## #########################################
                cat("applying ROI mask to cc matrix...", "\n")
                ptm <- proc.time()[3]
                map_cc_roi <- map_cc_slice

                mask.roi.temp <- mask.roi
                mask.roi.temp[mask.roi.temp==0] <- NA

                for(z in 1:nt)
                    map_cc_roi[,,z] <- map_cc_roi[,,z]*mask.roi.temp

                cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
                cat("--------", "\n")

                ## ##################################################
                ## ##EXTRACT THE WHOLE TUMOR PROFILE FROM ROI########
                ## ##################################################
                ptm <- proc.time()[3]
                cc.median <- seq(1:nt)
                for(t in 1:nt)
                    cc.median[t] <- median(map_cc_roi[,,t], na.rm=TRUE)

                cat("fitting", modeltype1, "and", modeltype2, "models to whole ROI data...", "\n")


                ## ESTIMATE INITIAL PARAMETER VALUES USING NUMERICAL DECONVOLUTION
                IRF.out <- calchFUNC(map.times, aif, cc.median)

                AUCcorrnom <- IRF.out$AUCh
                AUCMRTcorrnom <- IRF.out$AUChMRTh
                AUCMRTcorrnom.divby.AUCcorrnom <- IRF.out$AUChMRTh/IRF.out$AUCh

                ## INITIAL PARAM VALUES AND BOUNDS FOR TOFTS MODEL
                if(tlag.Tofts.on == TRUE & AIF.shift=="VEIN"){
                    p0.T.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, 0.12)
                    lower.wholeT <- c(lo, lo, 0)
                    upper.wholeT <- c(Inf,Inf,Inf)
                }
                if(tlag.Tofts.on == TRUE & AIF.shift=="ARTERY"){
                    p0.T.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, 0.12)
                    lower.wholeT <- c(lo, lo, 0)
                    upper.wholeT <- c(Inf,Inf,Inf)
                }
                if(tlag.Tofts.on == FALSE){
                    p0.T.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh)
                    lower.wholeT <- c(lo, lo)
                    upper.wholeT <- c(Inf,Inf)
                }


                ## INITIAL PARAM VALUES AND BOUNDS FOR XTOFTS MODEL
                if(AIF.shift=="VEIN"){
                    p0.xT.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, 0.05, 0.12)
                    lower.wholexT <- c(lo, lo, lo, 0)
                    upper.wholexT <- c(Inf,Inf,Inf,Inf)
                }
                if(AIF.shift=="ARTERY"){
                    p0.xT.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, 0.05, 0.12)
                    lower.wholexT <- c(lo, lo, lo, 0)
                    upper.wholexT <- c(Inf,Inf,Inf,Inf)
                }
                if(AIF.shift=="NONE"){
                    p0.xT.median <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, 0.05)
                    lower.wholexT <- c(lo, lo, lo)
                    upper.wholexT <- c(Inf,Inf,Inf)
                }

                ## GENERATE A FILE OF MEDIAN PROFILE OVER ROI FOR SAAM II
                SAAMII <- FALSE
                if(SAAMII == TRUE){
                    hw.x <- cbind(format(map.times, digits=3), format(aif, digits=3), format(cc.median, digits=3))

                    write.table("DATA", file=paste("slice", slice, "_forSAAMII", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
                    write.table("(SD 1)", file=paste("slice", slice, "_forSAAMII", sep=""), quote=FALSE, row.names=FALSE, append=TRUE, col.names=FALSE)
                    write.table(hw.x, file=paste("slice", slice, "_forSAAMII", sep=""), quote=FALSE, row.names=FALSE, append=TRUE, col.names=c("t", paste("AIF_s", slice, sep=""), paste("CC_s", slice, sep="")))
                    write.table("END", file=paste("slice", slice, "_forSAAMII", sep=""), quote=FALSE, row.names=FALSE, append=TRUE, col.names=FALSE)
                }


                ## ###############################################
                ## ##FIT SELECTED MODEL TO MEDIAN ROI DATA########
                ## ###############################################
                roi <- cc.median
                t <- map.times

                fix.tlag <- FALSE


                ## SWITCH TO THE EXTENDED TOFTS MODEL STRUCTURE############
                roi.model <- roi.modelxT

                if(method.optimization=="L-BFGS-B")
                    fit.roi.median.xTofts <- try(optim(p0.xT.median, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = "L-BFGS-B", lower=lower.wholexT, upper=upper.wholexT, hessian=TRUE), silent=try.silent)

                if(method.optimization!="L-BFGS-B")
                    fit.roi.median.xTofts <- try(optim(p0.xT.median, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = method.optimization, hessian=TRUE), silent=try.silent)

                if(class(fit.roi.median.xTofts) != "try-error"){
                    if(fit.roi.median.xTofts$convergence == 0){
                        roi.median.fitted.xTofts <- roi.model(p=fit.roi.median.xTofts$par, t=map.times, dt=diff(map.times), cp=aif)
                        params.xTofts <- fit.roi.median.xTofts$par

                        param.est.whole.roi.xTofts <- list(fit.roi.median.xTofts$par[1], fit.roi.median.xTofts$par[2], fit.roi.median.xTofts$par[3], fit.roi.median.xTofts$par[4])
                        names(param.est.whole.roi.xTofts) <- c("Ktrans", "kep", "vb", "tlag")
                    }
                }

                if(class(fit.roi.median.xTofts) == "try-error"){
                    cat("A problem occured when fitting the extended Tofts model to median intensity/concentration data across the ROI. Optimization algorithm did not converge. This may occur when a large fraction of voxels within the ROI are non-enhancing, as these are not excluded from analysis of the median intensity/concentration profile", "\n")
                }

                if(class(fit.roi.median.xTofts) != "try-error"){
                    if(fit.roi.median.xTofts$convergence != 0){
                        cat("A problem occured when fitting the Tofts model to median intensity/concentration data across the ROI (convergence code not equal to zero). Optimization algorithm did not converge. This may occur when a large fraction of voxels within the ROI are non-enhancing, as these are not excluded from analysis of the median intensity/concentration profile", "\n")
                    }
                }


                ## SWITCH TO TOFTS MODEL STRUCTURE#########################
                roi.model <- roi.modelT

                if(method.optimization=="L-BFGS-B")
                    fit.roi.median.Tofts <- try(optim(p0.T.median, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = "L-BFGS-B", lower=lower.wholeT, upper=upper.wholeT, hessian=TRUE), silent=try.silent)

                if(method.optimization!="L-BFGS-B")
                    fit.roi.median.Tofts <- try(optim(p0.T.median, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = method.optimization, hessian=TRUE), silent=try.silent)

                if(class(fit.roi.median.Tofts) != "try-error"){
                    if(fit.roi.median.Tofts$convergence == 0){
                        roi.median.fitted.Tofts <- roi.model(p=fit.roi.median.Tofts$par, t=map.times, dt=diff(map.times), cp=aif)

                        if(tlag.Tofts.on == FALSE){
                            param.est.whole.roi.Tofts <- list(fit.roi.median.Tofts$par[1], fit.roi.median.Tofts$par[2])
                            names(param.est.whole.roi.Tofts) <- c("Ktrans", "kep")
                        }

                        if(tlag.Tofts.on == TRUE){
                            param.est.whole.roi.Tofts <- list(fit.roi.median.Tofts$par[1], fit.roi.median.Tofts$par[2], fit.roi.median.Tofts$par[3])
                            names(param.est.whole.roi.Tofts) <- c("Ktrans", "kep", "tlag")
                        }
                    }
                }

                if(class(fit.roi.median.Tofts) == "try-error"){
                    cat("A problem occured when fitting the Tofts model to median intensity/concentration data across the ROI (a try-error occured). Optimization algorithm did not converge. This may occur when a large fraction of voxels within the ROI are non-enhancing, as these are not excluded from analysis of the median intensity/concentration profile", "\n")
                }

                if(class(fit.roi.median.Tofts) != "try-error"){
                    if(fit.roi.median.Tofts$convergence != 0){
                        cat("A problem occured when fitting the Tofts model to median intensity/concentration data across the ROI (convergence code not equal to zero). Optimization algorithm did not converge. This may occur when a large fraction of voxels within the ROI are non-enhancing, as these are not excluded from analysis of the median intensity/concentration profile", "\n")
                    }
                }

                roi.model <- roi.modelxT

                ## PERFORM TRUNCATION CORRECTION ON INITIAL PARAMETER VALUES ESTIMATED BY DECONVOLUTION
                ## AND USE AS INITIAL PARAMS FOR PER-VOXEL FITTING (BEGIN)
                IRF.out <- calchFUNC(vector.times=map.times, AIF=aif, map_cc_slice=cc.median, vp.nom=param.est.whole.roi.xTofts$vb, kep.nom=param.est.whole.roi.xTofts$kep)
                p0.T <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh)
                p0.xT <- c(IRF.out$AUChMRTh, IRF.out$AUChMRTh/IRF.out$AUCh, param.est.whole.roi.xTofts$vb)
                names(p0.T) <- c("Ktrans", "kep")
                names(p0.xT) <- c("Ktrans", "kep", "vb")
                ## PERFORM TRUNCATION CORRECTION ON INITIAL PARAMETER VALUES ESTIMATED BY DECONVOLUTION
                ## AND USE AS INITIAL PARAMS FOR PER-VOXEL FITTING (END)


                ## SAVE UNCORRECTED AND CORRECTED IRF RESULTS INTO A SINGLE R OBJECT
                AUCcorr <- IRF.out$AUCh
                AUCMRTcorr <- IRF.out$AUChMRTh
                AUCMRTcorr.divby.AUCcorr <- IRF.out$AUChMRTh/IRF.out$AUCh
                IRF.results <- c(AUCcorrnom, AUCMRTcorrnom, AUCMRTcorrnom.divby.AUCcorrnom, AUCcorr, AUCMRTcorr, AUCMRTcorr.divby.AUCcorr)
                names(IRF.results) <- c("AUCcorrnom(ve)", "AUCMRTcorrnom(Ktrans)", "AUCMRTcorrnom.divby.AUCcorrnom(kep)", "AUCcorr(ve)", "AUCMRTcorr(Ktrans)", "AUCMRTcorr.divby.AUCcorr(kep)")

                cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
                cat("--------", "\n")

                ## #CALCULATE %CVs FOR FITTED PARAMETERS (xTofts model only)########
                roi.model <- roi.modelxT

                if(class(fit.roi.median.xTofts) != "try-error"){
                    if(fit.roi.median.xTofts$convergence == 0){
                        RSS <- sum((roi.model(p=fit.roi.median.xTofts$par, t=map.times, dt=diff(map.times), cp=aif)- roi)^2)

                        param_est <- fit.roi.median.xTofts$par
                        nD <- nt
                        nP <- length(param_est)
                        hess <- fit.roi.median.xTofts$hessian
                        cov <- try((nP*RSS/(nD-nP))*solve(hess), silent=try.silent)
                        ## create placeholder cv object in case try-error occurs
                        cv <- param_est
                        cv[] <- NA

                        if(class(cov) != "try-error"){
                            sd <- try(sqrt(diag(cov)), silent=try.silent)
                            if(class(sd) != "try-error"){
                                cv <- sd/param_est*100
                                cv.whole.roi.xTofts <- as.numeric(format(cv, digits=1))
                            }
                        }

                        param.roi <- as.numeric(format(param_est, digits=3))
                    }
                }

                ## ####################################
                ## ##SAVE SHIFTED AIF/VIF##############
                ## ####################################
                if(AIF.shift=="ARTERY" | AIF.shift=="VEIN")
                    aif.shifted <- aif.shift.func(map.times, aif, param.est.whole.roi.xTofts$tlag)

                ## #########################################
                ## ##PLOT AIFs and MEDIAN CC PROFILE########
                ## #########################################
                if(show.rt.fits==TRUE){

                    if(AIF.shift=="ARTERY"){
                        dev.new(width=6, height=6, xpos=1500, ypos=0)
                        plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Arterial Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n")
                        lines(map.times, aif, col="red")
                        lines(map.times, aif.shifted)
                        legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw AIF", "shifted AIF"), c("red", "black"))
                    }
                    if(AIF.shift=="VEIN"){
                        dev.new(width=6, height=6, xpos=1500, ypos=0)
                        plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Venous Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n")
                        lines(map.times, aif, col="blue")
                        lines(map.times, aif.shifted)
                        legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw VIF", "shifted VIF"), c("blue", "black"))
                    }
                    if(AIF.shift=="NONE"){
                        dev.new(width=6, height=6, xpos=1500, ypos=0)
                        plot(map.times, aif, ylim=c(0, max(aif)), xlab="min", ylab="contrast agent", main="Vascular Input Function", type="n")
                        lines(map.times, aif)
                        legend(x=0.55*max(map.times), y=max(aif), c("raw VIF"), c("black"))
                    }
                }

                if(show.rt.fits==TRUE){
                    dev.new(width=6, height=6, xpos=1500, ypos=575)
                    plot(map.times, roi, xlab="min", ylab="contrast agent", main="median contrast agent conc; Tofts=blue, xTofts=red", cex=2)
                    text(x=0.65*max(map.times), y=0.40*max(roi, na.rm=TRUE), labels="EXTENDED TOFTS PARAMS", cex=1)

                    if(class(fit.roi.median.xTofts) != "try-error"){
                        if(fit.roi.median.xTofts$convergence == 0){
                            text(x=0.65*max(map.times), y=0.35*max(roi, na.rm=TRUE), labels=paste(expression(K^trans), " = ", format(param.est.whole.roi.xTofts$Ktrans, digits=3), " 1/min (", cv.whole.roi.xTofts[1], "%)", sep=""), cex=1)
                            text(x=0.65*max(map.times), y=0.30*max(roi, na.rm=TRUE), labels=paste(expression(k_ep), " = ", format(param.est.whole.roi.xTofts$kep, digits=3), " 1/min (", cv.whole.roi.xTofts[2], "%)", sep=""), cex=1)
                            text(x=0.65*max(map.times), y=0.25*max(roi, na.rm=TRUE), labels=paste(expression(v_b), " = ", format(param.est.whole.roi.xTofts$vb, digits=3), " dimensionless (", cv.whole.roi.xTofts[3], "%)", sep=""), cex=1)
                            if(AIF.shift != "NONE")
                                text(x=0.65*max(map.times), y=0.20*max(roi, na.rm=TRUE), labels=paste(expression(t_lag), " = ", format(60*param.est.whole.roi.xTofts$tlag, digits=3), " sec (", cv.whole.roi.xTofts[4], "%)", sep=""), cex=1)
                            lines(map.times, roi.median.fitted.xTofts, col="red", lwd=5)
                        }
                    }

                    if(class(fit.roi.median.Tofts) != "try-error"){
                        if(fit.roi.median.Tofts$convergence == 0){
                            lines(map.times, roi.median.fitted.Tofts, col="blue", lwd=2)
                        }
                    }
                }

            }

            ## ############################################################
            ## ##LOAD MATLAB FILE PREVIOUSLY GENERATED BY THIS SCRIPT######
            ## ############################################################
            if(file.original == FALSE){
                fix.tlag <- TRUE

                ## #IMPORT MATLAB FILE INTO R#######################
                if(file.format=="matlab"){
                    mat_data <- readMat(file)
                    args <- mat_data$args[,,1]
                    results_file <- args$resultsfile
                    ID.visit <- args$IDvisit
                    slice <- args$slice
                    method.optimization <- args$methodoptimization
                    show.rt.fits <- args$showrtfits
                    param.for.avdt <- param.for.avdt
                    range.map <- args$rangemap
                    cutoff.map <- args$cutoffmap
                    lo <- args$lo
                    p0.xT <- mat_data$p0xT
                    p0.T <- mat_data$p0T
                    AIF.shift <- args$AIFshift
                    nx <- mat_data$nx
                    ny <- mat_data$ny
                    nt <- mat_data$nt
                    map_cc_slice <- mat_data$cc
                    map_cc_roi <- mat_data$ccroi
                    map.times <- mat_data$maptimes
                    aif <- mat_data$aif
                    aif.shifted <- mat_data$aifshifted
                    map.KtransxT <- mat_data$mapKtransxT
                    map.kepxT <- mat_data$mapkepxT
                    map.vexT <- mat_data$mapvexT
                    map.vbxT <- mat_data$mapvbxT
                    map.KtransT.cv <- mat_data$mapKtransTcv
                    map.kepT.cv <- mat_data$mapkepTcv
                    map.KtransxT.cv <- mat_data$mapKtransxTcv
                    map.kepxT.cv <- mat_data$mapkepxTcv
                    map.vbxT.cv <- mat_data$mapvbxTcv
                    map.fitfailuresxT <- mat_data$mapfitfailuresxT
                    map.KtransT <- mat_data$mapKtransT
                    map.kepT <- mat_data$mapkepT
                    map.veT <- mat_data$mapveT
                    map.fitfailuresT <- mat_data$mapfitfailuresT
                    mask.roi <- mat_data$maskroi
                    param.est.medianT <- mat_data$paramestmedianT
                    param.est.medianxT <- mat_data$paramestmedianxT
                    cc.median <- mat_data$ccmedian
                    roi.median.fitted <- mat_data$roimedianfitted
                    param.est.whole.roi <- mat_data$paramestwholeroi
                    map.AIC.compare <- mat_data$mapAICcompare
                    map.AIC.T <- mat_data$mapAICT
                    map.EF <- mat_data$mapEF
                    map.AIC.xT <- mat_data$mapAICxT
                    roi.median.fitted.Tofts <- mat_data$roimedianfittedTofts
                    param.est.whole.roi.Tofts <- mat_data$paramestwholeroiTofts
                    roi.median.fitted.xTofts <- mat_data$roimedianfittedxTofts
                    param.est.whole.roi.xTofts <- mat_data$paramestwholeroixTofts
                    cv.whole.roi.xTofts <- mat_data$cvwholeroixTofts
                    params.xTofts <- mat_data$paramsxTofts

                    rm(mat_data)
                }

                if(file.format=="RData"){
                    load(file)
                    args <- dcemri.data$args
                    results_file <- args$resultsfile
                    ID.visit <- args$IDvisit
                    slice <- args$slice
                    method.optimization <- args$methodoptimization
                    show.rt.fits <- args$showrtfits
                    param.for.avdt <- param.for.avdt
                    range.map <- args$rangemap
                    cutoff.map <- args$cutoffmap
                    lo <- args$lo
                    p0.xT <- dcemri.data$p0xT
                    p0.T <- dcemri.data$p0T
                    AIF.shift <- args$AIFshift
                    nx <- dcemri.data$nx
                    ny <- dcemri.data$ny
                    nt <- dcemri.data$nt
                    export.matlab <- args$exportmatlab
                    map_cc_slice <- dcemri.data$cc
                    map_cc_roi <- dcemri.data$ccroi
                    map.times <- dcemri.data$maptimes
                    aif <- dcemri.data$aif
                    aif.shifted <- dcemri.data$aifshifted
                    map.KtransT.cv <- dcemri.data$mapKtransTcv
                    map.kepT.cv <- dcemri.data$mapkepTcv
                    map.KtransxT.cv <- dcemri.data$mapKtransxTcv
                    map.kepxT.cv <- dcemri.data$mapkepxTcv
                    map.vbxT.cv <- dcemri.data$mapvbxTcv
                    map.KtransxT <- dcemri.data$mapKtransxT
                    map.kepxT <- dcemri.data$mapkepxT
                    map.vexT <- dcemri.data$mapvexT
                    map.vbxT <- dcemri.data$mapvbxT
                    map.fitfailuresxT <- dcemri.data$mapfitfailuresxT
                    map.KtransT <- dcemri.data$mapKtransT
                    map.kepT <- dcemri.data$mapkepT
                    map.veT <- dcemri.data$mapveT
                    map.fitfailuresT <- dcemri.data$mapfitfailuresT
                    mask.roi <- dcemri.data$maskroi
                    param.est.medianT <- dcemri.data$paramestmedianT
                    param.est.medianxT <- dcemri.data$paramestmedianxT
                    cc.median <- dcemri.data$ccmedian
                    roi.median.fitted <- dcemri.data$roimedianfitted
                    param.est.whole.roi <- dcemri.data$paramestwholeroi
                    map.AIC.compare <- dcemri.data$mapAICcompare
                    map.AIC.T <- dcemri.data$mapAICT
                    map.EF <- dcemri.data$mapEF
                    map.AIC.xT <- dcemri.data$mapAICxT
                    roi.median.fitted.Tofts <- dcemri.data$roimedianfittedTofts
                    param.est.whole.roi.Tofts <- dcemri.data$paramestwholeroiTofts
                    roi.median.fitted.xTofts <- dcemri.data$roimedianfittedxTofts
                    param.est.whole.roi.xTofts <- dcemri.data$paramestwholeroixTofts
                    cv.whole.roi.xTofts <- dcemri.data$cvwholeroixTofts
                    params.xTofts <- dcemri.data$paramsxTofts
                    rm(dcemri.data)
                }
                modeltype1 <- "xTofts"
                modeltype2 <- "Tofts"
            }

            if(file.original == TRUE){
                ## ############################################################
                ## ##FIT THE EXTENDED TOFTS MODEL TO ROI VOXELS################
                ## ############################################################

                cat("fitting", modeltype1, "and", modeltype2, "models to ROI voxels...", "\n")

                ptm <- proc.time()[3]

                t <- map.times

                map.KtransxT <- matrix(NA, nrow=nx, ncol=ny)
                map.kepxT <- matrix(NA, nrow=nx, ncol=ny)
                map.vexT <- matrix(NA, nrow=nx, ncol=ny)
                map.vbxT <- matrix(NA, nrow=nx, ncol=ny)

                map.KtransxT.cv <- matrix(NA, nrow=nx, ncol=ny)
                map.kepxT.cv <- matrix(NA, nrow=nx, ncol=ny)
                map.vbxT.cv <- matrix(NA, nrow=nx, ncol=ny)

                map.fitfailuresxT <- matrix(NA, nrow=nx, ncol=ny)
                map.OptimValuexT <- matrix(NA, nrow=nx, ncol=ny)

                map.KtransT <- matrix(NA, nrow=nx, ncol=ny)
                map.kepT <- matrix(NA, nrow=nx, ncol=ny)
                map.veT <- matrix(NA, nrow=nx, ncol=ny)

                map.KtransT.cv <- matrix(NA, nrow=nx, ncol=ny)
                map.kepT.cv <- matrix(NA, nrow=nx, ncol=ny)

                map.fitfailuresT <- matrix(NA, nrow=nx, ncol=ny)
                map.OptimValueT <- matrix(NA, nrow=nx, ncol=ny)

                map.AIC.compare <- matrix(0, nrow=nx, ncol=ny)
                map.AIC.T <- matrix(NA, nrow=nx, ncol=ny)
                map.EF <- matrix(NA, nrow=nx, ncol=ny)
                map.AIC.xT <- matrix(NA, nrow=nx, ncol=ny)

                cc_fittedxT <- array(0, dim=c(nx,ny,nt))
                cc_fittedT <- array(0, dim=c(nx,ny,nt))

                nv <- 1
                nv1_q <- trunc(quantile(1:length(mask.roi[mask.roi==1]), probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)))

                if(show.rt.fits==TRUE)
                    dev.new(xpos=3500, ypos=0)

                ptm_slice <- proc.time()[3]

                ## #FUNCTION TO CALCULATE FRACTION OF VALUES IN A VECTOR THAT ARE GREATER THAN ZERO
                GTzero <- function(x){
                    length(as.vector(x>0)[as.vector(x>0)==TRUE])/length(x)
                }
                ## ###########

                for(x in 1:nx){
                    for(y in 1:ny){
                        if(mask.roi[x,y] == 1 & GTzero(map_cc_slice[x,y,])<=fracGTzero){
                            map.fitfailuresxT[x,y] <- -2
                            map.fitfailuresT[x,y] <- -2
                        }

                        if(mask.roi[x,y] == 1 & GTzero(map_cc_slice[x,y,])>fracGTzero)
                            map.EF[x,y] <- 1

                        if(mask.roi[x,y] == 1 & GTzero(map_cc_slice[x,y,])>fracGTzero){
                            nv <- nv+1

                            roi <- map_cc_slice[x,y,]

                            fix.tlag <- TRUE

                            roi.model <- roi.modelxT

                            if(verbose==TRUE){
                                cat("x =", x, "\n")
                                cat("y =", y, "\n")
                                cat("contrast agent curve =", roi, "\n")
                                cat("fitting xTofts model to voxel data...")
                            }

                            if(method.optimization=="L-BFGS-B")
                                fit_roi <- try(optim(p0.xT, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = "L-BFGS-B", lower=c(lo, lo, lo), upper=c(Inf,Inf,Inf), hessian=TRUE), silent=try.silent)

                            if(method.optimization!="L-BFGS-B")
                                fit_roi <- try(optim(p0.xT, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = method.optimization, hessian=TRUE), silent=try.silent)

                            if(verbose==TRUE)
                                cat("done", "\n")

                            if(class(fit_roi) != "try-error"){
                                if(fit_roi$convergence == 0){

                                    ## CALCULATE %CVs (begin)

                                    RSS <- sum((roi.model(p=fit_roi$par, t=map.times, dt=diff(map.times), cp=aif)- roi)^2)
                                    param_est <- fit_roi$par
                                    nD <- nt
                                    nP <- length(param_est)
                                    hess <- fit_roi$hessian
                                    cov <- try((nP*RSS/(nD-nP))*solve(hess), silent=try.silent)

                                    if(class(cov) != "try-error"){
                                        sd <- try(sqrt(diag(cov)), silent=try.silent)
                                        if(class(sd) != "try-error"){
                                            cv <- sd/param_est*100
                                            cv.roi.voxel <- as.numeric(format(cv, digits=1))
                                            map.KtransxT.cv[x,y] <- cv.roi.voxel[1]
                                            map.kepxT.cv[x,y] <- cv.roi.voxel[2]
                                            map.vbxT.cv[x,y] <- cv.roi.voxel[3]
                                        }
                                    }
                                    ## CALCULATE %CVs (end)

                                    map.KtransxT[x,y] <- fit_roi$par[1]
                                    map.kepxT[x,y] <- fit_roi$par[2]
                                    if(map.kepxT[x,y] < 0.00001)
                                        map.kepxT[x,y] <- 0.00001
                                    map.vexT[x,y] <- fit_roi$par[1]/map.kepxT[x,y]
                                    map.vbxT[x,y] <- fit_roi$par[3]
                                    map.OptimValuexT[x,y] <- fit_roi$value
                                }
                            }

                            if(class(fit_roi) == "try-error"){
                                if(verbose==TRUE)
                                    cat("...but a <try-error> occurred", "\n")
                                map.fitfailuresxT[x,y] <- 99
                            }

                            if(class(fit_roi) != "try-error")
                                map.fitfailuresxT[x,y] <- fit_roi$convergence

                            if(class(fit_roi) == "try-error"){
                                map.KtransxT[x,y] <- NA
                                map.kepxT[x,y] <- NA
                                map.vexT[x,y] <- NA
                                map.vbxT[x,y] <- NA
                            }

                            if(class(fit_roi) != "try-error"){
                                if(fit_roi$convergence > 0){
                                    map.KtransxT[x,y] <- NA
                                    map.kepxT[x,y] <- NA
                                    map.vexT[x,y] <- NA
                                    map.vbxT[x,y] <- NA
                                }
                            }


                            if(verbose==TRUE)
                                cat("simulating xTofts model at estimated parameter values...")


                            if(class(fit_roi) != "try-error"){
                                if(fit_roi$convergence == 0){
                                    simulation <- roi.model(p=fit_roi$par, t=map.times, dt=diff(map.times), cp=aif)
                                    cc_fittedxT[x,y,] <- simulation
                                }
                            }


                            if(verbose==TRUE)
                                cat("done", "\n")

                            ##
                            roi.model <- roi.modelT
                            ##

                            if(verbose==TRUE)
                                cat("fitting Tofts model to voxel data...")

                            if(method.optimization=="L-BFGS-B")
                                fit_roiT <- try(optim(p0.T, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = "L-BFGS-B", lower=lower.wholeT, upper=upper.wholeT, hessian=TRUE), silent=try.silent)

                            if(method.optimization!="L-BFGS-B")
                                fit_roiT <- try(optim(p0.T, Obj_roi, model=roi.model, t=map.times, dt=diff(map.times), cp=aif, roi=roi, method = method.optimization, hessian=TRUE), silent=try.silent)

                            if(verbose==TRUE)
                                cat("done", "\n")

                            if(class(fit_roiT) != "try-error") {
                                if(fit_roiT$convergence == 0){

                                    ## CALCULATE %CVs (begin)

                                    RSS <- sum((roi.model(p=fit_roiT$par, t=map.times, dt=diff(map.times), cp=aif)- roi)^2)

                                    nD <- nt
                                    nP <- length(param_est)
                                    param_est <- fit_roiT$par
                                    df <- nt - length(param_est)
                                    hess <- fit_roiT$hessian
                                    cov <- try((nP*RSS/(nD-nP))*solve(hess), silent=try.silent)

                                    if(class(cov) != "try-error"){
                                        sd <- try(sqrt(diag(cov)), silent=try.silent)
                                        if(class(sd) != "try-error"){
                                            cv <- sd/param_est*100
                                            cv.roi.voxel <- as.numeric(format(cv, digits=1))
                                            map.KtransT.cv[x,y] <- cv.roi.voxel[1]
                                            map.kepT.cv[x,y] <- cv.roi.voxel[2]
                                        }
                                    }
                                    ## CALCULATE %CVs (end)

                                    map.KtransT[x,y] <- fit_roiT$par[1]
                                    map.kepT[x,y] <- fit_roiT$par[2]
                                    if(map.kepT[x,y] < 0.00001)
                                        map.kepT[x,y] <- 0.00001
                                    map.veT[x,y] <- fit_roiT$par[1]/map.kepT[x,y]
                                    map.OptimValueT[x,y] <- fit_roiT$value
                                }
                            }

                            if(class(fit_roiT) == "try-error"){
                                map.fitfailuresT[x,y] <- 99
                                if(verbose==TRUE)
                                    cat("...but a <try-error> occurred", "\n")
                            }

                            if(class(fit_roiT) == "try-error"){
                                map.KtransT[x,y] <- NA
                                map.kepT[x,y] <- NA
                                map.veT[x,y] <- NA
                            }

                            if(class(fit_roiT) != "try-error"){
                                if(fit_roiT$convergence > 0){
                                    map.KtransT[x,y] <- NA
                                    map.kepT[x,y] <- NA
                                    map.veT[x,y] <- NA
                                }
                            }

                            if(class(fit_roiT) != "try-error")
                                map.fitfailuresT[x,y] <- fit_roiT$convergence

                            if(verbose==TRUE)
                                cat("simulating Tofts model at estimated parameter values...")

                            if(class(fit_roiT) != "try-error"){
                                if(fit_roiT$convergence==0){
                                    simulation_2 <- roi.model(p=fit_roiT$par, t=map.times, dt=diff(map.times), cp=aif)
                                    cc_fittedT[x,y,] <- simulation_2
                                }
                            }

                            if(verbose==TRUE)
                                cat("done", "\n")

                            roi.model <- roi.modelxT

                            if(class(fit_roi) != "try-error" & class(fit_roiT) != "try-error"){
                                if(fit_roi$convergence == 0 & fit_roiT$convergence == 0){

                                    if(verbose==TRUE)
                                        cat("calculating AICc values...")

                                    simulation_AIC <- simulation
                                    simulation_AIC_2 <- simulation_2


                                    ## TRY AIC INSTEAD OF LOG_LIKELIHOOD TEST (BEGIN)
                                    np_1 <- length(fit_roiT$par)
                                    np_2 <- length(fit_roi$par)  + 1 ## Add 1 to account for tlag which, which although not estimated on a per-voxel basis is estimated based on whole ROI data and is fixed in the per-voxel fits and should be penalized in the per-voxel fits

                                    ## Corrected AIC from Glatting 2007
                                    AIC1 <- nt*log(fit_roiT$value) + 2*(np_1+1) + (2*(np_1+1)*(np_1+2))/(nt-np_1-2)
                                    AIC2 <- nt*log(fit_roi$value) + 2*(np_2+1) + (2*(np_2+1)*(np_2+2))/(nt-np_2-2)
                                    AIC1 <- as.numeric(format(AIC1, digits=1))
                                    AIC2 <- as.numeric(format(AIC2, digits=1))
                                    if(AIC2 < AIC1)
                                        map.AIC.compare[x,y] <- 1
                                    if(AIC2 >= AIC1)
                                        map.AIC.compare[x,y] <- 2
                                    map.AIC.T[x,y] <- AIC1
                                    map.AIC.xT[x,y] <- AIC2
                                    ## TRY AIC INSTEAD OF LOG_LIKELIHOOD TEST (END)
                                }
                            }

                            if(class(fit_roi) != "try-error"){
                                if(fit_roi$convergence > 0){
                                    map.AIC.compare[x,y] <- NA
                                }
                            }

                            if(class(fit_roiT) != "try-error"){
                                if(fit_roiT$convergence > 0){
                                    map.AIC.compare[x,y] <- NA
                                }
                            }

                            if(class(fit_roi) == "try-error" | class(fit_roiT) == "try-error")
                                map.AIC.compare[x,y] <- NA


                            if(verbose==TRUE){
                                cat("done", "\n")
                                cat("======================", "\n")
                            }


                            ## ##########################

                            if(show.rt.fits==TRUE){
                                if(class(fit_roiT) != "try-error" & class(fit_roi) != "try-error"){
                                    if(fit_roiT$convergence==0 & fit_roi$convergence==0){
                                        plot(map.times, roi, ylab="contrast agent", xlab="min", main=paste(modeltype1, "(red) and", modeltype2, "(blue)"), cex=3)
                                        lines(map.times, simulation, col="red", lwd=5)
                                        lines(map.times, simulation_2, col="blue", lwd=2)
                                    }
                                }
                            }

                            if(nv==2)
                                cat("progress: ")
                            if(nv==nv1_q[[1]])
                                cat("10%..")
                            if(nv==nv1_q[[2]])
                                cat("20%..")
                            if(nv==nv1_q[[3]])
                                cat("30%..")
                            if(nv==nv1_q[[4]])
                                cat("40%..")
                            if(nv==nv1_q[[5]])
                                cat("50%..")
                            if(nv==nv1_q[[6]])
                                cat("60%..")
                            if(nv==nv1_q[[7]])
                                cat("70%..")
                            if(nv==nv1_q[[8]])
                                cat("80%..")
                            if(nv==nv1_q[[9]])
                                cat("90%..", "\n")
                            if(nv==nv1_q[[10]]-10)
                                cat("..10")
                            if(nv==nv1_q[[10]]-9)
                                cat("..9..")
                            if(nv==nv1_q[[10]]-8)
                                cat("8..")
                            if(nv==nv1_q[[10]]-7)
                                cat("7..")
                            if(nv==nv1_q[[10]]-6)
                                cat("6..")
                            if(nv==nv1_q[[10]]-5)
                                cat("5..")
                            if(nv==nv1_q[[10]]-4)
                                cat("4..")
                            if(nv==nv1_q[[10]]-3)
                                cat("3..")
                            if(nv==nv1_q[[10]]-2)
                                cat("2..")
                            if(nv==nv1_q[[10]]-1)
                                cat("1..", "\n")

                        }

                    }
                }

                cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
                cat("--------", "\n")

                graphics.off()

                KtransxT.median <- median(map.KtransxT[map.KtransxT>0 & map.fitfailuresxT==0], na.rm=TRUE)
                kepxT.median <- median(map.kepxT[map.kepxT>0 & map.fitfailuresxT==0], na.rm=TRUE)
                vexT.median <- median(map.vexT[map.vexT>0 & map.fitfailuresxT==0], na.rm=TRUE)
                vbxT.median <- median(map.vbxT[map.vbxT>=0 & map.fitfailuresxT==0], na.rm=TRUE)
                fitfailuresxT.total <- length(map.fitfailuresxT[map.fitfailuresxT >= 1]) / length(map.fitfailuresxT[map.fitfailuresxT >= 0]) * 100
                param.est.medianxT <- list(KtransxT.median, kepxT.median, vexT.median, vbxT.median, fitfailuresxT.total)
                names(param.est.medianxT) <- c("Ktrans.median", "kep.median", "ve.median", "vb.median", "percent.fitfailures")


                KtransT.median <- median(map.KtransT[map.KtransT>0 & map.fitfailuresT==0], na.rm=TRUE)
                kepT.median <- median(map.kepT[map.kepT>0 & map.fitfailuresT==0], na.rm=TRUE)
                veT.median <- median(map.veT[map.veT>0 & map.fitfailuresT==0], na.rm=TRUE)
                fitfailuresT.total <- length(map.fitfailuresT[map.fitfailuresT >= 1]) / length(map.fitfailuresT[map.fitfailuresT >= 0]) * 100
                param.est.medianT <- list(KtransT.median, kepT.median, veT.median, fitfailuresT.total)
                names(param.est.medianT) <- c("Ktrans.median", "kep.median", "ve.median", "percent.fitfailures")
            }


            ## ############################################################
            ## ##GENERATE A UNIQUE FILENAME BASED ON THE CURRENT DATE######
            ## ############################################################
            if(file.original==TRUE){
                if(file=="concatenate.KAT.with.KAT.checkData.RData")
                    ID.visit <- paste(strsplit(results_file, split=".mat")[[1]], "_s", slice, sep="")
                if(file!="concatenate.KAT.with.KAT.checkData.RData")
                    ID.visit <- paste(strsplit(file, split=".mat")[[1]], "_s", slice, sep="")
            }
            if(file.original==FALSE)
                ID.visit <- strsplit(file, split=".mat")[[1]]

            ID.visit <- strsplit(ID.visit, split="/")
            ID.visit <- ID.visit[[1]][length(ID.visit[[1]])]

            IDvp <- strsplit(ID.visit, split="_")
            ID.visit.forplot <- paste(IDvp[[1]][1], ".", IDvp[[1]][2], ".", IDvp[[1]][3], ".", IDvp[[1]][4], sep="")

            DATE <- date()
            if(length(strsplit(DATE,split="  ")[[1]])==2){
                full_date <- strsplit(DATE,split="  ")[[1]]
                full_date_1 <- full_date[1]
                full_date_1 <- strsplit(full_date_1, split=" ")[[1]]
                full_date_2 <- full_date[2]
                full_date_2 <- strsplit(full_date_2, split=" ")[[1]]
                month <- full_date_1[2]
                day   <- full_date_2[1]
                time  <- full_date_2[2]
                year  <- full_date_2[3]
            }
            if(length(strsplit(DATE,split="  ")[[1]])==1){
                full_date <- strsplit(DATE,split=" ")[[1]]
                month <- full_date[2]
                day <- full_date[3]
                time <- full_date[4]
                year <- full_date[5]
            }
            year <- strsplit(year, split="")[[1]]
            year <- paste(year[3],year[4],sep="")
            time_concat <- strsplit(time, split=":")[[1]]
            time_concat <- paste(paste(time_concat[1], time_concat[2], sep=""), time_concat[3], sep="")
            DATE <- paste(paste(paste(paste(day,month,sep=""), year, sep=""), "-", sep=""), time_concat, sep="")

            filename3 <- paste(ID.visit, "_KAT_", DATE, ".mat", sep="")
            filename3 <- sub(".RData", "", filename3)


            if(file.original==FALSE){
                roi.model <- roi.modelxT
            }

            ## ##NOTE THAT CALCULATION OF x_max, x_min, y_max, y_min IS PERFORMED WHEN file.original=TRUE AND file.original=FALSE##

            if(param.for.avdt == "Ktrans")
                MAP <- map.KtransxT
            if(param.for.avdt == "kep")
                MAP <- map.kepxT
            if(param.for.avdt == "ve")
                MAP <- map.vexT
            if(param.for.avdt == "vb")
                MAP <- map.vbxT

            ## ##ZOOM IN ON TUMOR ROI#####
            x_min <- 0
            x_max <- nx
            y_min <- 0
            y_max <- ny

            for(xx in 1:nx){
                if(sum(MAP[xx,], na.rm=TRUE)>0){
                    x_min <- xx - 3
                    break
                }
                if(sum(MAP[xx,], na.rm=TRUE)==0)
                    xx <- xx+1
            }

            for(y in 1:ny){
                if(sum(MAP[,y], na.rm=TRUE)>0){
                    y_min <- y - 3
                    break
                }
                if(sum(MAP[,y], na.rm=TRUE)==0)
                    y <- y+1
            }

            for(xx in nx:0){
                if(sum(MAP[xx,], na.rm=TRUE)>0){
                    x_max <- xx + 3
                    break
                }
                if(sum(MAP[xx,], na.rm=TRUE)==0)
                    xx <- xx-1
            }

            for(y in ny:0){
                if(sum(MAP[,y], na.rm=TRUE)>0){
                    y_max <- y + 3
                    break
                }
                if(sum(MAP[,y], na.rm=TRUE)==0)
                    y <- y-1
            }

            MAP_ul_s1 <- MAP[MAP>0]
            MAP_ul_s1 <- sort(MAP_ul_s1)
            MAP_ul <- range.map*(max(MAP_ul_s1[1:length(MAP_ul_s1)*cutoff.map]))

            if(file.original==FALSE){

                if(param.for.avdt == "Ktrans")
                    MAP <- map.KtransxT
                if(param.for.avdt == "kep")
                    MAP <- map.kepxT
                if(param.for.avdt == "ve")
                    MAP <- map.vexT
                if(param.for.avdt == "vb")
                    MAP <- map.vbxT


                if(file.format=="matlab"){
                    if(param.for.avdt == "Ktrans")
                        param.median <- param.est.medianxT[,,1]$Ktrans.median[1]

                    if(param.for.avdt == "kep")
                        param.median <- param.est.medianxT[,,1]$kep.median[1]

                    if(param.for.avdt == "ve")
                        param.median <- param.est.medianxT[,,1]$ve.median[1]

                    if(param.for.avdt == "vb")
                        param.median <- param.est.medianxT[,,1]$vb.median[1]

                    if(param.for.avdt == "tlag")
                        param.median <- param.est.medianxT[,,1]$tlag.median[1]
                }


                if(file.format=="RData"){
                    if(param.for.avdt == "Ktrans")
                        param.median <- param.est.medianxT$Ktrans.median

                    if(param.for.avdt == "kep")
                        param.median <- param.est.medianxT$kep.median

                    if(param.for.avdt == "ve")
                        param.median <- param.est.medianxT$ve.median

                    if(param.for.avdt == "vb")
                        param.median <- param.est.medianxT$vb.median

                    if(param.for.avdt == "tlag")
                        param.median <- param.est.medianxT$tlag.median
                }

                ## ##CALCULATE RANGE OF COLOR BAR (LEGEND)#####
                MAP_for_plot <- MAP

                ## ##REPLACE ANY NEGATIVE NUMBERS WITH ZEROS######
                MAP_for_plot[MAP_for_plot<0] <- 0

                ## ##TRUNCATE HIGH OUTLIER VALUES#####
                MAP_for_plot[MAP_for_plot>=MAP_ul*0.99] <- MAP_ul*0.99

                dev.new(width=6, height=6, xpos=1860, ypos=578)
                plot(map.times, cc.median, xlab = "min", ylab = "contrast agent", main = paste(modeltype1, "(red) and", modeltype2, "(blue) fitted to median whole ROI data"), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 2)

                lines(map.times, roi.median.fitted.xTofts, col = "red", lwd=5)
                lines(map.times, roi.median.fitted.Tofts, col = "blue", lwd=2)

                text(x=0.6*max(map.times), y=0.30*max(cc.median, na.rm=TRUE), labels="EXTENDED TOFTS PARAMS", cex=1)
                text(x=0.6*max(map.times), y=0.24*max(cc.median, na.rm=TRUE), labels=paste(expression(K^trans), " = ", format(param.est.whole.roi.xTofts$Ktrans, digits=3), " 1/min (", cv.whole.roi.xTofts[1], "%)", sep=""), cex=1)
                text(x=0.6*max(map.times), y=0.18*max(cc.median, na.rm=TRUE), labels=paste(expression(k_ep), " = ", format(param.est.whole.roi.xTofts$kep, digits=3), " 1/min (", cv.whole.roi.xTofts[2], "%)", sep=""), cex=1)
                text(x=0.6*max(map.times), y=0.12*max(cc.median, na.rm=TRUE), labels=paste(expression(v_b), " = ", format(param.est.whole.roi.xTofts$vb, digits=3), " dimensionless (", cv.whole.roi.xTofts[3], "%)", sep=""), cex=1)
                if(AIF.shift != "NONE")
                    text(x=0.6*max(map.times), y=0.06*max(cc.median, na.rm=TRUE), labels=paste(expression(t_lag), " = ", format(60*param.est.whole.roi.xTofts$tlag, digits=3), " sec (", cv.whole.roi.xTofts[4], "%)", sep=""), cex=1)

                dev.new(width=6, height=6, xpos=1860, ypos=0)
                image(map.AIC.compare, xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main=paste("best fit:", modeltype1, "(red) and", modeltype2, "(blue)"))
                image(map.AIC.compare, xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main=paste("best fit:", modeltype1, "(red) and", modeltype2, "(blue)"))


                dev.new(width=12.75, height=12.75, xpos=238, ypos=0)

                image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste("Model Type=", modeltype1, "     Median ", param.for.avdt,"=", format(param.median, digit=2), "     VisitID=", strsplit(args$IDvisit, split="_")[[1]][1], "     slice=", args$slice, sep=""))
                image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste("Model Type=", modeltype1, "     Median ", param.for.avdt,"=", format(param.median, digit=2), "     VisitID=", strsplit(args$IDvisit, split="_")[[1]][1], "     slice=", args$slice, sep=""))

                text(x_min/nx+(x_max/nx-x_min/nx)*0.97, y_min/ny+(y_max/ny-y_min/ny)*0.99, "close", col="green")
                text(x_min/nx+(x_max/nx-x_min/nx)*0.05, y_min/ny+(y_max/ny-y_min/ny)*0.99, "print to PDF", col="green")
                text(x_min/nx+(x_max/nx-x_min/nx)*0.84, y_min/ny+(y_max/ny-y_min/ny)*0.01, paste("KAT for DCEMRI v", KAT.version, ", Genentech PTPK", sep=""), col="darkgrey")
                legend <- seq(0,MAP_ul,by=0.001)
                dim(legend) <- c(1,length(legend))

                dev.new(width=2.5, height=12.75, xpos=0, ypos=0)

                image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAP_ul), xaxt="n", yaxt="n", xlab="0", main=format(MAP_ul, digits=3), ylab=expression(min^-1), cex.lab=1.25, cex.main=1.05)
                image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAP_ul), xaxt="n", yaxt="n", xlab="0", main=format(MAP_ul, digits=3), ylab=expression(min^-1), cex.lab=1.25, cex.main=1.05)


                if(AIF.shift=="ARTERY"){
                    dev.new(width=5.15, height=4.4, xpos=1500, ypos=0)
                    plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Arterial Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n")
                    lines(map.times, aif, col="red")
                    lines(map.times, aif.shifted)
                    legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw AIF", "shifted AIF"), c("red", "black"))
                }
                if(AIF.shift=="VEIN"){
                    dev.new(width=5.15, height=4.4, xpos=1500, ypos=0)
                    plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Venous Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n")
                    lines(map.times, aif, col="blue")
                    lines(map.times, aif.shifted)
                    legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw VIF", "shifted VIF"), c("blue", "black"))
                }
                if(AIF.shift=="NONE"){
                    dev.new(width=5.15, height=4.4, xpos=1500, ypos=0)
                    plot(map.times, aif, ylim=c(0, max(aif)), xlab="min", ylab="contrast agent", main="Vascular Input Function", type="n")
                    lines(map.times, aif)
                    legend(x=0.55*max(map.times), y=max(aif), c("raw VIF"), c("black"))
                }


                ## #######VOXEL LOCATOR####################
                dev.set(4)

                inf <- 1
                newplot <- 1
                newplot2 <- 1
                legend_count<-1
                legend_labels <- 1:1000
                legend_matrix <- matrix(0, ncol=ny, nrow=nx)


                cat("---", "\n")

                while(inf == 1){

                    z <- locator(1, type="o", col="green")

                    xx <- round(z$x*(nx-1)+1)
                    yy <- round(z$y*(ny-1)+1)

                    if(legend_matrix[xx,yy]==0){

                        legend(z$x, z$y, legend_labels[legend_count], col="green", text.col="green")
                        legend_matrix[xx,yy] <- 1

                        cat("Voxel Number/Coordinates: n=", legend_count, ", x=", xx, ", y=", yy, "\n", sep="")
                        cat("Parameter estimates (xTofts): Ktrans=", format(map.KtransxT[xx,yy],digits=3), ", ve=", format(map.KtransxT[xx,yy]/map.kepxT[xx,yy], digits=3), ", vb=", format(map.vbxT[xx,yy],digits=3), "\n", sep="")
                        cat("%CVs of Parameter estimates (xTofts): Ktrans=", format(map.KtransxT.cv[xx,yy],digits=2), ", kep (Ktrans/ve)=", format(map.kepxT.cv[xx,yy], digits=2), ", vb=", format(map.vbxT.cv[xx,yy],digits=2), "\n", sep="")
                        cat("Parameter estimates (Tofts): Ktrans=", format(map.KtransT[xx,yy],digits=3), ", ve=", format(map.KtransT[xx,yy]/map.kepT[xx,yy], digits=3), "\n", sep="")
                        cat("%CVs of Parameter estimates (Tofts): Ktrans=", format(map.KtransT.cv[xx,yy],digits=2), ", kep (Ktrans/ve)=", format(map.kepT.cv[xx,yy], digits=3), "\n", sep="")
                        cat("---", "\n")

                        legend_count <- legend_count+1
                    }

                    xdim <- x_max-x_min
                    ydim <- y_max-y_min

                    if(xx > (x_max-0.1*xdim) && yy > (y_max-0.02*ydim)){
                        graphics.off()
                        dev.off()
                    }

                    xx_old <- xx
                    yy_old <- yy

                    ## ####PLOT DATA AND SIMULATION############

                    conc <- 1:nt
                    for(i in 1:nt)
                        conc[i] <- map_cc_slice[xx,yy,i]

                    if(xx < (x_min + 0.12*xdim) && yy > (y_max-0.02*ydim)){

                        ## ##########PRINT SUMMARY FIGURES TO PDF FILE#################
                        pdf(file=paste(results_file, "-SUMMARY.pdf", sep=""), height=12, width=15)

                        layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow=TRUE), widths=c(1.5,5.5,5.5))
                        par(omi=c(0.15, 0.15, 0.15, 0.15))

                        image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAP_ul), xaxt="n", yaxt="n", xlab="0", main=format(MAP_ul, digits=3), ylab=expression(min^-1), cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
                        image(legend, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), zlim=c(0,MAP_ul), xaxt="n", yaxt="n", xlab="0", main=format(MAP_ul, digits=3), ylab=expression(min^-1), cex.lab=1.5, cex.main=1.5, cex.axis=1.5, add=TRUE)

                        image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste("Model Type=", modeltype1, "     Median ", param.for.avdt,"=", format(param.median, digit=2), "     VisitID=", strsplit(args$IDvisit, split="_")[[1]][1], "     slice=", args$slice, sep=""), cex.axis=1.5)

                        image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste("Model Type=", modeltype1, "     Median ", param.for.avdt,"=", format(param.median, digit=2), "     VisitID=", strsplit(args$IDvisit, split="_")[[1]][1], "     slice=", args$slice, sep=""), cex.axis=1.5, add=TRUE)

                        plot(map.times, cc.median, xlab="min", ylab="contrast agent", main="median contrast agent conc; Tofts=blue, xTofts=red", cex=3, cex.axis=1.5, cex.lab=1.5)
                        text(x=0.6*max(map.times), y=0.25*max(cc.median, na.rm=TRUE), labels="EXTENDED TOFTS PARAMS", cex=1.5)
                        text(x=0.6*max(map.times), y=0.20*max(cc.median, na.rm=TRUE), labels=paste(expression(K^trans), " = ", format(param.est.whole.roi.xTofts$Ktrans, digits=3), " 1/min (", cv.whole.roi.xTofts[1], "%)", sep=""), cex=1.5)
                        text(x=0.6*max(map.times), y=0.15*max(cc.median, na.rm=TRUE), labels=paste(expression(v_e), " = ", format(param.est.whole.roi.xTofts$ve, digits=3), " 1/min (", cv.whole.roi.xTofts[2], "%)", sep=""), cex=1.5)
                        text(x=0.6*max(map.times), y=0.10*max(cc.median, na.rm=TRUE), labels=paste(expression(v_b), " = ", format(param.est.whole.roi.xTofts$vb, digits=3), " dimensionless (", cv.whole.roi.xTofts[3], "%)", sep=""), cex=1.5)
                        if(AIF.shift != "NONE")
                            text(x=0.6*max(map.times), y=0.05*max(cc.median, na.rm=TRUE), labels=paste(expression(t_lag), " = ", format(60*param.est.whole.roi.xTofts$tlag, digits=3), " sec (", cv.whole.roi.xTofts[4], "%)", sep=""), cex=1.5)
                        lines(map.times, roi.median.fitted.xTofts, col="red", lwd=5)
                        lines(map.times, roi.median.fitted.Tofts, col="blue", lwd=2)
                        text(x=0.5*max(map.times), 1*max(cc.median, na.rm=TRUE),  paste("R package version =", KAT.version), cex=1.5)

                        image(array(1:2, dim=c(1,2)), col=palette(colorRampPalette(c("red","blue"))(2)), zlim=c(0,2), xaxt="n", yaxt="n", xlab=modeltype1, main=modeltype2, cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
                        image(array(1:2, dim=c(1,2)), col=palette(colorRampPalette(c("red","blue"))(2)), zlim=c(0,2), xaxt="n", yaxt="n", xlab=modeltype1, main=modeltype2, cex.lab=1.5, cex.main=1.5, cex.axis=1.5, add=TRUE)

                        image(map.AIC.compare, xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main="Summary of model discrimination analysis", cex.axis=1.5, cex.lab=1.5)
                        image(map.AIC.compare, xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), col=palette(colorRampPalette(c("black", "red","blue"))(3)), main="Summary of model discrimination analysis", cex.axis=1.5, cex.lab=1.5, add=TRUE)

                        if(AIF.shift=="ARTERY"){
                            plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Arterial Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n", cex=3, cex.axis=1.5, cex.lab=1.5)
                            lines(map.times, aif, col="red")
                            lines(map.times, aif.shifted)
                            legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw AIF", "shifted AIF"), c("red", "black"), cex=2)
                        }
                        if(AIF.shift=="VEIN"){
                            plot(map.times, aif, ylim=c(0, max(aif, aif.shifted)), xlab="min", ylab="contrast agent", main=paste("Shifted Venous Input Function (", format(param.est.whole.roi.xTofts$tlag*60, digits=3), " sec)", sep=""), type="n", cex=3, cex.axis=1.5, cex.lab=1.5)
                            lines(map.times, aif, col="blue")
                            lines(map.times, aif.shifted)
                            legend(x=0.55*max(map.times), y=max(aif, aif.shifted), c("raw VIF", "shifted VIF"), c("blue", "black"), cex=2)
                        }
                        if(AIF.shift=="NONE"){
                            plot(map.times, aif, ylim=c(0, max(aif)), xlab="min", ylab="contrast agent", main="Vascular Input Function", type="n", cex=3, cex.axis=1.5, cex.lab=1.5)
                            lines(map.times, aif)
                            legend(x=0.55*max(map.times), y=max(aif), c("raw VIF"), c("black"), cex=2)
                        }

                        dev.off()

                        cat("image printed to pdf.", "\n")
                        cat("---", "\n")
                    }


                    if(newplot==1){
                        dev.new(width=5.15, height=4.4, xpos=1500, ypos=432)
                        newplot <- 2
                    }

                    dev.set(7)


                    plot(map.times, conc, xlab="min", ylab="contrast agent", ylim=c(-max(conc, na.rm=TRUE)/5, 1.4*max(conc, na.rm=TRUE)), cex=1.5, main=paste(paste("red=", modeltype1, ", blue=", modeltype2, " (", sep=""), paste("x=", round(xx_old), ", y=", round(yy_old), ")", sep=""), sep=""))

                    value_xTofts <- MAP[xx,yy]

                    if(is.finite(value_xTofts)==FALSE)
                        value_xTofts <- 0

                    if(value_xTofts != 0){

                        paramsxT <- c(map.KtransxT[xx,yy], map.kepxT[xx,yy], map.vbxT[xx,yy], param.est.whole.roi.xTofts$tlag)
                        paramsT <- c(map.KtransT[xx,yy], map.kepT[xx,yy])

                        roi.model <- roi.modelxT

                        simulation <- roi.model(p=paramsxT, t=map.times, dt=diff(map.times), cp=aif)

                        roi.model <- roi.modelT

                        simulation_2<- roi.model(p=paramsT, t=map.times, dt=diff(map.times), cp=aif)

                        lines(map.times, simulation_2, col="blue", lwd=2)
                        lines(map.times, simulation, col="red", lwd=2, lty=2)

                        roi.model <- roi.modelxT
                    }


                    text(max(map.times)/4.5, 1.4*max(conc, na.rm=TRUE), "Fitted xTofts params")
                    text(max(map.times)/4.5, 1.4*max(conc, na.rm=TRUE)*0.93, paste(paste("Ktrans =", format(map.KtransxT[xx,yy], digits=3)),"min^-1"))
                    text(max(map.times)/4.5, 1.4*max(conc, na.rm=TRUE)*0.86, paste("ve =", format(map.vexT[xx,yy], digits=3)))
                    text(max(map.times)/4.5, 1.4*max(conc, na.rm=TRUE)*0.79, paste("vb =", format(map.vbxT[xx,yy], digits=3)))
                    if(newplot2==1){
                        dev.new(width=5.15, height=3, xpos=1500, ypos=864)
                        newplot2 <- 2
                    }
                    else
                        dev.set(8)

                    if(value_xTofts != 0){
                        simulation_AIC <- simulation
                        simulation_AIC_2 <- simulation_2

                        if(length(simulation)==length(conc)){
                            plot(map.times,conc-simulation, xlab="min", ylab="predicated - measured", cex=1.5, col="red", main=paste("AIC(Tofts)=", map.AIC.T[xx,yy], " AIC(xTofts)=", map.AIC.xT[xx,yy], sep=""))
                            lines(locfit(conc-simulation~map.times, acri="ici"), col="red", lwd=2)
                            abline(h=0, col="black", lwd=2)

                            if(is.finite(simulation_2[1])==TRUE){
                                points(map.times,conc-simulation_2, xlab="min", cex=1.5, col="blue")
                                lines(locfit(conc-simulation_2~map.times, acri="ici"), col="blue", lwd=2, lty=2)
                            }
                        }
                    }
                    dev.set(4)
                }
            }

            if(file.original==TRUE){
                ptm <- proc.time()[3]

                ## ######################
                ## ##EXPORT RESULTS######
                ## ######################
                proc.time.total <- format((proc.time()[3]-ptm_total)/60, digits=2)
                args <- list(as.character(file), as.character(results_file), as.character(method.optimization), show.rt.fits, as.character(param.for.avdt), range.map, cutoff.map, export.matlab, export.RData, verbose, show.errors, try.silent, fracGTzero, AIF.shift, slice, ID.visit)

                names(args) <- c("file", "resultsfile", "methodoptimization", "showrtfits", "paramforavdt", "rangemap", "cutoffmap", "exportmatlab", "exportRData", "verbose", "showerrors", "trysilent", "fracGTzero", "AIFshift", "slice", "IDvisit")

                roiplotparams <- list(x_min, x_max, y_min, y_max, MAP_ul)
                names(roiplotparams) <- c("xmin", "xmax", "ymin", "ymax", "MAPul")

                dummy_data <- dcemri.data

                dcemri.data <- list(args, map_cc_slice, map_cc_roi, cc.median, map.times, aif, aif.shifted, mask.roi, map.KtransxT, map.KtransxT.cv, map.kepxT, map.kepxT.cv, map.vbxT, map.vbxT.cv, map.vexT, map.OptimValuexT, map.fitfailuresxT, param.est.medianxT, roi.median.fitted.xTofts, param.est.whole.roi.xTofts, cv.whole.roi.xTofts, map.KtransT, map.KtransT.cv, map.kepT, map.kepT.cv, map.veT, map.OptimValueT, map.fitfailuresT, param.est.medianT, roi.median.fitted.Tofts, param.est.whole.roi.Tofts, proc.time.total, roiplotparams, KAT.version, map.AIC.xT, map.AIC.T, map.AIC.compare, nx, ny, nt, cc_fittedxT, cc_fittedT, p0.T, p0.xT, IRF.results, map.EF)

                names(dcemri.data) <- c("args", "cc", "ccroi", "ccmedian", "maptimes", "aif", "aifshifted", "maskroi", "mapKtransxT", "mapKtransxTcv", "mapkepxT", "mapkepxTcv", "mapvbxT", "mapvbxTcv", "mapvexT", "mapOptimValuexT", "mapfitfailuresxT", "paramestmedianxT", "roimedianfittedxTofts", "paramestwholeroixTofts", "cvwholeroixTofts", "mapKtransT", "mapKtransTcv", "mapkepT", "mapkepTcv", "mapveT", "mapOptimValueT", "mapfitfailuresT", "paramestmedianT", "roimedianfittedTofts", "paramestwholeroiTofts", "proctimetotal", "roiplotparams", "KATversion", "mapAICxT", "mapAICT", "mapAICcompare", "nx", "ny", "nt", "ccfittedxT", "ccfittedT", "p0T", "p0xT", "IRFresults", "mapEF")

                ## ##Use user specified absolute path; "results_file" argument ######
                if(export.RData==TRUE){
                    cat("writing results to ",  paste(results_file, ".RData", sep=""), "...", sep="", "\n")
                    save(dcemri.data, file=paste(results_file, ".RData", sep=""))
                }

                if(export.matlab==TRUE){
                    ## ##Use user specified absolute path; "results_file" argument ######
                    cat("writing results to ",  paste(results_file, ".mat", sep=""), "...", sep="", "\n")
                    writeMat(paste(results_file, ".mat", sep=""), args=args, mapccslice=map_cc_slice, mapccroi=map_cc_roi, ccmedian=cc.median, maptimes=map.times, aif=aif, aifshifted=aif.shifted, maskroi=mask.roi, mapKtransxT=map.KtransxT, mapKtransxTcv=map.KtransxT.cv, mapkepxT=map.kepxT, mapkepxTcv=map.kepxT.cv, mapvbxT=map.vbxT, mapvbxTcv=map.vbxT.cv, mapvexT=map.vexT, mapOptimValuexT=map.OptimValuexT, mapfitfailuresxT=map.fitfailuresxT, paramestmedianxT=param.est.medianxT, roimedianfittedxTofts=roi.median.fitted.xTofts, paramestwholeroixTofts=param.est.whole.roi.xTofts, cvwholeroixTofts=cv.whole.roi.xTofts, mapKtransT=map.KtransT, mapKtransTcv=map.KtransT.cv, mapkepT=map.kepT, mapkepTcv=map.kepT.cv, mapveT=map.veT, mapOptimValueT=map.OptimValueT, mapfitfailuresT=map.fitfailuresT, paramestmedianT=param.est.medianT, roimedianfittedTofts=roi.median.fitted.Tofts, paramestwholeroiTofts=param.est.whole.roi.Tofts, proctimetotal=proc.time.total, roiplotparams=roiplotparams, KATversion=KAT.version, mapAICxT=map.AIC.xT, mapAICT=map.AIC.T, mapAICcompare=map.AIC.compare, nx=nx, ny=ny, nt=nt, ccfittedxT=cc_fittedxT, ccfittedT=cc_fittedT, p0T=p0.T, p0xT=p0.xT, IRFresults=IRF.results, mapEF=map.EF)
                }

                dcemri.data <- dummy_data

                cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
                if(export.matlab==FALSE)
                    cat("Run KAT(filename.RData) to visualize results.", "\n")
                if(export.matlab==TRUE)
                    cat("Run KAT(filename.RData) or KAT(filename.mat) to visualize results.", "\n")
                cat("--------", "\n")
            }
        }
    }
}
