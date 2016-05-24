##' Detect rhythmic signals from time-series datasets with individual
##'   information
##'
##' This is a function that takes use of any one method from ARSER, JTK_CYCLE
##'   and Lomb-Scargle to detect rhythmic signals from time-series datasets
##'   containing individual information.
##'
##' This function is originally aimed to analyze large scale periodic data with
##'   individual information. Please pay attention to the data format of
##'   \code{datafile} and \code{designfile}(see \code{Examples} part).
##'   Time-series experimental values(missing values as \code{NA}) from
##'   all individuals should be stored in \code{datafile}, with the first row
##'   containing all library ID(unique identification number for each sample)
##'   and the first column containing all detected molecular names(eg.
##'   transcript or gene name). The \code{designfile} should at least have
##'   three columns-library ID, subject ID and sampling time column.
##'   Experimental group information of each subject ID may be in another
##'   column. In addition, sampling time information may be stored in multiple
##'   columns instead of one column. For example, sampling time-"36 hours" may
##'   be recorded as "day 2"(sampling day column, \code{design_dayColm}) plus
##'   "12 hours"(sampling hour column, \code{design_hrColm}). The library ID
##'   in \code{datafile} and \code{designfile} should be same. If there are
##'   different characters between library ID in these two files, try
##'   \code{design_libIDrename} to keep them same.
##'
##'   \code{ARS}, \code{JTK} or \code{LS} could be used to analyze time-series
##'   profiles individual by individual. \code{meta3d} requires that all
##'   individuals should be analyzed by the same method before integrating
##'   calculated p-value, period, phase, baseline value, amplitude and relative
##'   amplitude values group by group. However, the sampling pattern among
##'   individuals may be different and the requirement of sampling pattern for
##'   each method is not same(see more information about these methods and their
##'   limitations in \code{\link{meta2d}}). Please carefully select a proper
##'   method for the specific dataset. \code{meta3d} also help users select
##'   the suitable method through warning notes.
##'
##'   P-values from different individuals are integrated with Fisher's method
##'   (\code{"fisher"})(Fisher,1925; implementation code from \pkg{MADAM}).For
##'   short time-series profiles(eg. 10 time points or less), p-values given by
##'   Lomb-Scargle may be over conservative, which will also lead to
##'   conservative integrated p-values. The integrated period, baseline,
##'   amplitude and relative amplitude values are arithmetic mean of multiple
##'   individuals, respectively. The phase is
##'   \href{https://en.wikipedia.org/wiki/Mean_of_circular_quantities}{
##'   mean of circular quantities}(\code{adjustPhase = "predictedPer"})
##'   or a arithmetic mean (\code{adjustPhase = "notAdjusted"}) of multiple
##'   individual phases. For completely removing the potential problem of
##'   averaging phases with quite different period length(also mentioned
##'   in \code{\link{meta2d}}), setting \code{minper}, \code{maxper} and
##'   \code{ARSdefaultPer} to a same value may be the only known way. If
##'   \code{weightedMethod = TRUE} is selected, weighted scores(
##'   \code{-log10(p-values)}) will be taken into account in integrating
##'   period, phase, baseline, amplitude and relative amplitude.
##'
##' @param datafile a character string. The name of data file containing
##'   time-series experimental values of all individuals.
##' @param designfile a character string. The name of experimental design file,
##'   at least containing the library ID(column names of \code{datafile}),
##'   subject ID(the individual corresponding to each library ID), and
##'   sampling time information of each library ID.
##' @param outdir a character string. The name of directory used to store
##'   output files.
##' @param filestyle a character vector(length 1 or 3). The data format of
##'   input files, must be \code{"txt"}, or \code{"csv"}, or a character
##'   vector containing field separator character(\code{sep}), quoting
##'   character(\code{quote}), and the character used for decimal
##'   points(\code{dec}, for details see \code{\link[utils]{read.table}}).
##' @param design_libColm a numeric value. The order index(from left to right)
##'   of the column storing library ID in \code{designfile}.
##' @param design_subjectColm a numeric value. The order index(from left to
##'   right) of the column storing subject ID in \code{designfile}.
##' @param minper a numeric value. The minimum period length of interested
##'   rhythms. The default is \code{20} for circadian rhythms.
##' @param maxper a numeric value. The maximum period length of interested
##'   rhythms. The default is \code{28} for circadian rhythms.
##' @param cycMethodOne a character string. The selected method for analyzing
##'   time-series data of each individual, must be one of \code{"ARS"}(ARSER),
##'   \code{"JTK"}(JTK_CYCLE), or \code{"LS"}(Lomb-Scargle).
##' @param timeUnit a character string. The basic time-unit, must be one of
##'   \code{"day"}, \code{"hour"}(default for circadian study),
##'   \code{"minute"}, or \code{"second"} depending on specific experimental
##'   design.
##' @param design_hrColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling hour information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}.
##' @param design_dayColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling day information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_minColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling minute information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_secColm a numeric value. The order index(from left to right)
##'   of the column storing time point value-sampling second information in
##'   \code{designfile}. If there is no such column in \code{designfile},
##'   set it as \code{NULL}(default).
##' @param design_groupColm a numeric value. The order index(from left to
##'   right) of the column storing experimental group information of each
##'   individual in \code{designfile}. If there is no such column in
##'   \code{designfile}, set it as \code{NULL}(default) and take all
##'   individuals as one group.
##' @param design_libIDrename a character vector(length 2) containing  a
##'   matchable character string in each library ID of \code{designfile}, and
##'   a replacement character string. If it is not necessary to replace
##'   characters in library ID of \code{designfile}, set it as \code{NULL}(
##'   default).
##' @param adjustPhase a character string. The method used to adjust each
##'   calculated phase before getting integrated phase, must be one of
##'   \code{"predictedPer"}(adjust phase with predicted period length)
##'   or \code{"notAdjusted"}(not adjust phase).
##' @param combinePvalue a character string. The method used to integrate
##'   p-values of multiple individuals, currently only \code{"fisher"}(
##'   Fisher's method) could be selected.
##' @param weightedMethod logical. If \code{TRUE}(default), weighted score
##'   based on p-value of each individual will be used to integrate period,
##'   phase and amplitude values of multiple individuals.
##' @param outIntegration a character string. This parameter controls what
##'   kinds of analysis results will be outputted, must be one of
##'   \code{"both"}, \code{"onlyIntegration"}, or \code{"noIntegration"}.
##'   See \code{\link{meta2d}} for more information.
##' @param ARSmle a character string. The strategy of using MLE method in
##'   \code{"ARS"}, must be one of \code{"auto"}, \code{"mle"}, or
##'   \code{"nomle"}. See \code{\link{meta2d}} for more information.
##' @param ARSdefaultPer a numeric value. The expected period length of
##'   interested rhythm, which is a necessary parameter for \code{ARS}. See
##'   \code{\link{meta2d}} for more information.
##' @param dayZeroBased logical. If \code{TRUE}, the first sampling day is
##'   recorded as day zero in the \code{designfile}.
##' @param outSymbol a character string. A common prefix exists in the names of
##'   output files.
##' @return
##' \code{meta3d} will write analysis results to \code{outdir} instead of
##'   returning them as objects. Output files with "meta3dSubjectID" in
##'   the file name are analysis results for each individual. Files named with
##'   "meta3dGroupID" store integrated p-values, period, phase, baseline,
##'   amplitude and relative amplitude values from multiple individuals of
##'   each group and calculated FDR values based on integrated p-values.
##' @references
##' Glynn E. F., Chen J., and Mushegian  A. R. (2006). Detecting periodic
##'   patterns in unevenly spaced gene expression time series using
##'   Lomb-Scargle periodograms. \emph{Bioinformatics}, \bold{22(3)},
##'   310--316
##'
##' Fisher, R.A. (1925). \emph{Statistical methods for research workers}.
##'   Oliver and Boyd (Edinburgh).
##'
##' Kugler K. G., Mueller L.A., and Graber A. (2010). MADAM - an open source
##'   toolbox for meta-analysis. \emph{Source Code for Biology and Medicine},
##'   \bold{5}, 3.
##'
##' @examples
##' # write 'cycHumanBloodData' and 'cycHumanBloodDesign' into two 'csv' files
##' write.csv(cycHumanBloodData, file="cycHumanBloodData.csv",
##'   row.names=FALSE)
##' write.csv(cycHumanBloodDesign, file="cycHumanBloodDesign.csv",
##'   row.names=FALSE)
##'
##' # detect circadian transcripts with JTK in studied individuals
##' meta3d(datafile="cycHumanBloodData.csv", cycMethodOne="JTK",
##'   designfile="cycHumanBloodDesign.csv", outdir="example",
##'   filestyle="csv", design_libColm=1, design_subjectColm=2,
##'   design_hrColm=4, design_groupColm=3)
##' @export

meta3d <- function(datafile, designfile, outdir="metaout", filestyle,
                    design_libColm, design_subjectColm, minper=20, maxper=28,
                    cycMethodOne="JTK", timeUnit="hour", design_hrColm,
                    design_dayColm=NULL, design_minColm=NULL, design_secColm=NULL,
                    design_groupColm=NULL, design_libIDrename=NULL,
                    adjustPhase="predictedPer", combinePvalue="fisher",
                    weightedMethod=TRUE, outIntegration="both", ARSmle="auto",
                    ARSdefaultPer=24, dayZeroBased=FALSE, outSymbol="")
{

    ####start time
    run_start=proc.time();

    ####create the directory of storing output files, if it is not exist
    outdir2 <- unlist(strsplit(outdir,.Platform$file.sep));
    OUTDIR <- paste(outdir2,collapse=.Platform$file.sep);
    if (! file.exists(OUTDIR) )
    { dir.create(OUTDIR); }

    ####check the input parameters
    ###extract the field separator character(FILE_SEP), the set of quoting
    ###characters(FILE_QUOTE), the character used for decimal points(FILE_DEC).
    ###the 'getFileSignF()' is in 'metaSubF.R'
    file_sign <- getFileSignF(filestyle);
    FILE_SEP <- file_sign$sep;
    FILE_QUOTE <- file_sign$quote;
    FILE_DEC <- file_sign$dec;
    FILE_QUOTE2 <- file_sign$quote2;

    ###check 'minper' and 'maxper'
    if ( (length(minper) != 1) | (length(maxper) != 1) )
    {   stop("The 'minper' or 'maxper' contains more than one element.\n");  }
    if ( (!is.numeric(minper)) | (minper <= 0 ) )
    {	stop("The 'minper' should be a positive numeric value.\n");	}
    if ( (!is.numeric(maxper)) | (maxper <= 0 ) )
    {	stop("The 'maxper' should be a positive numeric value.\n");	}
    if (minper > maxper)
    {	stop("The 'minper' should not be larger than 'maxper'.\n");	}

    ###check whether output the integration results at the end of analysis
    if ( (outIntegration == "both") | (outIntegration == "onlyIntegration") | (outIntegration == "noIntegration") )  {
        INTEGRATION <- outIntegration;
    } else {
        stop(c("Please check the parameter of 'outIntegration', it should ",
               "be set as 'both', 'onlyIntegration' or 'noIntegration'.\n") );
    }

    ###check the expected period length used to adjust phase
    if (adjustPhase == "predictedPer") {
        ADPHAL <- "PERIOD";
    } else if (adjustPhase == "notAdjusted") {
        ADPHAL <- 0;
    } else {
        stop("The 'adjustPhase' should be set as 'predictedPer' or 'notAdjusted'.\n");
    }

    ###check 'combinePvalue' and 'weightedMethod'
    if ( (combinePvalue != "Fisher") & (combinePvalue != "fisher") )
    {
        stop("The 'combinePvalue' should be set as 'Fisher' (equal to 'fisher').\n");
    }
    if ( !is.logical(weightedMethod) )
    {
        stop("The 'weightedMethod' should be set as a logical value (TRUE or FALSE).\n");
    }

    ###check 'cycMethodOne'
    if (length(cycMethodOne) == 1)
    {
        if ( (cycMethodOne != "ARS") & (cycMethodOne != "JTK") & (cycMethodOne != "LS") )
        {
            stop(c("Please check the parameter of 'cycMethodOne', one method ",
                   "could only be selected from 'ARS', 'JTK' and 'LS'.\n") );
        }
    }  else {
        stop(c("One and only one method should be selected from 'ARS', 'JTK' or 'LS'.\n") );
    }

    ####extract the essential information in the designfile
    ###check the 'designfile'
    designD <- read.table(file=designfile, header=TRUE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
    IDS <- dimnames(designD)[[1]];
    if (ncol(designD) < 3)
    {
        stop(c("The data in 'designfile' is stored in a data frame with less ",
               "than three columns ('design_libColm', 'design_subjectColm', ",
               "'design_hrColm' are necessary), please check whether ",
               "'filestyle' is correctly set.\n") );
    }
    if (nrow(designD) == 0)
    {
        stop(c("The 'designfile' contains only one line which is taken as ",
               "column names, please check the input file.\n") );
    }

    ###extract time information
    timeT <- dayT <- hrT <- minT <- secT <- rep(0, length(IDS));
    ###hour information; 'extractTimeColF()' in 'meta3dSubF.R'
    if (length(design_hrColm) > 0)
    {  hrT <- extractTimeColF(design_timeColm=design_hrColm, notekeys=c("design_hrColm", "hour"), designD=designD, IDS=IDS);  }
    ###day information
    if (length(design_dayColm) > 0)
    {  dayT <- extractTimeColF(design_timeColm=design_dayColm, notekeys=c("design_dayColm", "day"), designD=designD, IDS=IDS);  }
    ###minute information
    if (length(design_minColm) > 0)
    {  minT <- extractTimeColF(design_timeColm=design_minColm, notekeys=c("design_minColm", "minute"), designD=designD, IDS=IDS);  }
    ###second information
    if (length(design_secColm) > 0)
    {  secT <- extractTimeColF(design_timeColm=design_secColm, notekeys=c("design_secColm", "second"), designD=designD, IDS=IDS);  }
    ###extract numeric values storing in time information
    ###the function of 'getNumber()' is in 'meta3dSubF.R'.
    if ( !all(is.numeric(dayT)) )
    {    dayT <- getNumber(dayT);   }
    if ( !all(is.numeric(hrT)) )
    {    hrT <- getNumber(hrT);   }
    if ( !all(is.numeric(minT)) )
    {    minT <- getNumber(minT);   }
    if ( !all(is.numeric(secT)) )
    {    secT <- getNumber(secT);   }

    ###transfer multiple columns into one time column,
    ###and the default time unit is 'hour'
    if ( (timeUnit != "day") & (timeUnit != "hour") & (timeUnit != "minute") & (timeUnit != "second") )
    {
        stop("Please set 'timeUnit' as 'day', 'hour', 'minute', or 'second'.\n");
    }
    if ( !is.logical(dayZeroBased) )
    {
        stop("The 'dayZeroBased' should be set as a logical value (TRUE or FALSE).\n");
    }
    unitT <- c(24, 1, 1/60, 1/3600);
    if (timeUnit == "day")
    {    unitT <- c( 1, 1/24, 1/(60*24), 1/(3600*24) );   }
    if (timeUnit == "minute")
    {   unitT <- c(24*60, 60, 1, 1/60);    }
    if (timeUnit == "second")
    {   unitT <- c(24*3600, 3600, 60, 1);    }
    if ( (!dayZeroBased) & (sum(dayT)) )
    { dayT <- dayT - 1;  }
    timeT <- cbind(dayT*unitT[1], hrT*unitT[2], minT*unitT[3], secT*unitT[4]);
    ###check whether negative value existed in timeT
    negti <- NULL;
    for (ti in 1:nrow(timeT))
    {
        teptm <- timeT[ti,];
        if (length(teptm[teptm < 0]) > 0)
        {    negti <- c(negti, ti);    }
    }
    if (length(negti) > 0)
    {
        stop(c("The below lines contains negative time information",
                paste(", NegLine", negti, sep=""), ", please modify time ",
                "information in these lines before running this function.\n") );
    }
    timeT <- apply(timeT, 1, sum);

    ###check 'design_libColm' and 'design_subjectColm'
    ###the 'checkColmF()' and 'checkColmF2()' in 'meta3dSubF.R'.
    if ( (checkColmF(design_Colm=design_libColm, notekeys="design_libColm")) &
         (checkColmF(design_Colm=design_subjectColm, notekeys="design_subjectColm")) )
    {
        checkColmF2(design_Colm = design_libColm, notekeys = "design_libColm", designD_ncol = ncol(designD) )
        checkColmF2(design_Colm = design_subjectColm, notekeys = "design_subjectColm", designD_ncol = ncol(designD) )
    }

    ###extract libraryID, subjectID
    libID <- designD[IDS, design_libColm];
    subjectID <- designD[IDS, design_subjectColm];
    subjectID <- gsub("\\s+", "", subjectID);

    ###if 'design_libIDrename' is not NULL, modify the libraryID
    if (length(design_libIDrename))
    {
        if (length(design_libIDrename) == 2) {
            libID <- gsub(design_libIDrename[1], design_libIDrename[2], libID);
        }  else  {
            stop(c("For trying to replace some characters in the library names ",
                   "in designfile, please set 'design_libIDrename' as a character ",
                   "vector containing two elements, the first element is a matchable ",
                   "string, and the second is replacement character string.\n") );
        }
    }

    ####divide subjectID into different groups if design_groupColm is not set as '0'.
    groupID <- NULL;
    UNI_GROUPID <- NULL;
    GROUP_SUBJECTL <- list();
    if (length(design_groupColm) > 0)
    {
        ###check 'design_groupColm'
        if (checkColmF(design_Colm = design_groupColm, notekeys="design_groupColm")) {
            checkColmF2(design_Colm = design_groupColm, notekeys = "design_groupColm", designD_ncol = ncol(designD) )
        }
        ###extract 'groupID' information, and renamed 'subjectID'
        groupID <- designD[IDS, design_groupColm];
        groupID <- gsub("\\s+", "", groupID);
        subjectID <- paste(subjectID, groupID, sep="");
        names(groupID) <- subjectID;
        ###store subjectID group by group
        UNI_GROUPID <- unique(groupID);
        for (k in 1:length(UNI_GROUPID))
        {
            group_subjectID <- groupID[groupID == UNI_GROUPID[k]];
            group_subjectID <- unique(names(group_subjectID));
            GROUP_SUBJECTL[[k]] <- group_subjectID;
        }
    }

    ####store extracted information into a dataframe-DESIGND
    DESIGND <- NULL;
    if (sum(timeT) > 0)
    {
        if ( length(timeT) == nrow(designD) )  {
            if ( length(groupID) )
            {
                DESIGND <- data.frame("libID"=libID, "subjectID"=subjectID, "groupID"=groupID, "timeT"=timeT, row.names=IDS, stringsAsFactors=FALSE);
            }  else  {
                DESIGND <- data.frame("libID"=libID, "subjectID"=subjectID, "timeT"=timeT, row.names=IDS, stringsAsFactors=FALSE);
            }
        }  else {
            stop("There is an unknown bug in this package, please contact the author.\n");
        }
    }  else if (sum(timeT) == 0) {
        stop(c("There is no time information in the design file ",
               "or has not set any time column before running this function.\n") );
    }  else  {
        stop("There is an unknown bug in this package, please contact the author.\n");
    }

    ####divide libraryID into different subgroups according to the subjectID
    UNI_SUBJECTID <- unique(subjectID);
    SUBJECT_DGL <- list();
    for (i in 1:length(UNI_SUBJECTID))
    {
        subject_index <- which(subjectID == UNI_SUBJECTID[i]);
        SUBJECT_DGL[[i]] <- DESIGND[subject_index,];
    }

	####check whether selected 'cycMethodOne' suitable for all subjectID
	###get sampling pattern of time-series data subject by subject
	###the 'checkDataTypeF()' is in 'meta3dSubF.R'.
	patternD <- checkDataTypeF(subject_dgl = SUBJECT_DGL);

	###minimal time points required for periodic analysis
	minimal_timepoints <- 3;
	methodFlag <- rep(TRUE, 3);
	names(methodFlag) <- c("ARS", "JTK", "LS");
	if (all(patternD$uniNum >= minimal_timepoints))  {
		dup_index <- which(patternD$dup == TRUE);
		miss_index <- which(patternD$miss == TRUE);
		uneven_index <- which(patternD$uneven == TRUE);
		nonInteger_index <- which(patternD$nonInteger == TRUE);

		if ( (length(dup_index) > 0) | (length(miss_index) > 0) |
		     (length(uneven_index) > 0) | (length(nonInteger_index) > 0) )
		{
			methodFlag["ARS"] <- FALSE;
		}
		if ( (length(uneven_index) > 0) | (length(nonInteger_index) > 0) )
		{
			methodFlag["JTK"] <- FALSE;
		}

		if (methodFlag[cycMethodOne] == FALSE)
		{
			stop(c(cycMethodOne, " could not be used to analyze all subjects. ",
			       "See introduction about method selection in 'Vignette' file ",
				   "for information. Please select one method from ",
				   paste(names(methodFlag[methodFlag == TRUE]), collapse=", "), ".\n" ))
		}
	}  else  {
		delete_index <- which(patternD$uniNum < minimal_timepoints);
		stop(c("Before running this function, please remove records in 'datafile' ",
		       "and 'designfile' associated with these subjects: ",
			    paste(UNI_SUBJECTID[delete_index], collapse=", "),
		       ". The unique time points in these subjects are smaller ",
			   "than required minimal time points.\n"));
	}

    ####analyze time-series dataset subject by subject
    ###read the 'datafile', and check whether correctly read the whole dataset
    EXPD <- read.table(file=datafile, header=TRUE, sep=FILE_SEP, quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
    EXP_IDS <- dimnames(EXPD)[[1]];

    ###check whether missing values exist in 'datafile' when 'cycMethodOne = "ARS" '
    if ( cycMethodOne == "ARS")
    {
        nomissFlag <- apply(EXPD[,-1], 2, function(z) {  all(!is.na(z))  } );
        if (!all(nomissFlag))
        {
            stop(c("'ARS' could not be used to analyze datasets with missing ",
			       "values. Please select one method from 'JTK' and 'LS'.\n"));
        }
    }

    ###check whether 'datafile' is correctly read-in
    if (ncol(EXPD) == 1)
    {
        stop(c("The 'datafile' is stored in a one-column dataframe, ",
               "please check whether 'filestyle' is correctly set.\n") );
    }
    if (nrow(EXPD) == 0)
    {
        stop(c("The 'datafile' contains only one line which is taken ",
               "as column names, please check the input file.\n") );
    }

    ###check whether number of libraryID in 'datafile' and 'designfile' is same
    if ( nrow(DESIGND) != ( ncol(EXPD) - 1) )
    {
        stop(c("Please check the 'designfile' and 'datafile', the column number ",
               "of 'datafile' should be one more than the number of total samples ",
               "in 'designfile'. See help files or example data for more information.\n") );
    }

    ###check whether only one row exists in 'datafile'
    SIG_row <- FALSE;
    if (nrow(EXPD) == 1)
    {
        EXPD <- rbind(EXPD, EXPD);
        dimnames(EXPD)[[1]] <- c("1","2");
        EXP_IDS <- c("1","2");
        SIG_row <- TRUE;
    }

    ###store the name list in 'EXP_IDNAME'
    EXP_IDNAME <- EXPD[,1];
    names(EXP_IDNAME) <- EXP_IDS;

    ###run 'meta2d()' function subject by subject
    SUBJECT_LSL <- list();
    for (i in 1:length(UNI_SUBJECTID))
    {
        subject_designD <- SUBJECT_DGL[[i]];
        subject_libID <- subject_designD$libID;
        subject_timeT <- subject_designD$timeT;
        ###sort time points in order
        names(subject_timeT) <- 1:length(subject_timeT);
        order_timeT <- names(sort(subject_timeT));
        order_timeT <- as.numeric(order_timeT);
        ###arrange libraryID of each subject according sorted time points
        expdataD <- EXPD[,subject_libID[order_timeT]];
        dimnames(expdataD)[[1]] <- EXP_IDS;
        dimnames(expdataD)[[2]] <- subject_libID[order_timeT];
        cat(paste("The 'meta3d' is processing ", UNI_SUBJECTID[i],
                  ", the ", i, " in total ", length(UNI_SUBJECTID),
                  " subjects.\n", sep="") );

        meta2d_outD <- runmeta2dF(indata=expdataD, intime=subject_timeT[order_timeT],
                            datatype=patternD[i,], rundir=OUTDIR, runMethod=cycMethodOne,
                            minper=minper, maxper=maxper, ARSmle=ARSmle, ARSdefaultPer=ARSdefaultPer);
        SUBJECT_LSL[[i]] <- meta2d_outD;
    }

    ###check whether there are analysis results
    if (!length(SUBJECT_LSL))
    {
        stop(c("No output result is from the analysis about '",
               datafile, "', please carefully check the 'datafile' ",
               "'designfile' and self-setted parameters in this function.\n") );
    }

	####output separate analysis result for each subject-'SUBJECT_LSL'
	###get common part of output files
	filename <- unlist(strsplit(designfile,.Platform$file.sep,fixed=TRUE));
	filename <- filename[length(filename)];
	if ( grepl("\\", filename, fixed=TRUE) )
	{
		filename <- unlist(strsplit(designfile,"\\",fixed=TRUE));
		filename <- filename[length(filename)];
	}  else if ( grepl("//", filename, fixed=TRUE) )  {
		filename <- unlist(strsplit(designfile,"//",fixed=TRUE));
		filename <- filename[length(filename)];
	}

	###output analysis results for each subject
	if ( INTEGRATION != "onlyIntegration" )
	{
		for ( j in 1:length(SUBJECT_LSL) )
		{
			subject_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "meta3dSubjectID_", UNI_SUBJECTID[j], "_", filename, sep="");
			subject_outputD <- SUBJECT_LSL[[j]];
			subject_nameorder <- as.character(subject_outputD[,"CycID"]);
			subject_outputD[,"CycID"] <- EXP_IDNAME[subject_nameorder];
			if (!SIG_row)
			{
				write.table(subject_outputD, file=subject_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
			}  else  {
				sigsubject_outputD <- subject_outputD;
				###the second and third column is pvalue and fdr, respectively
				sigsubject_outputD[,3] <- sigsubject_outputD[,2];
				write.table(sigsubject_outputD[1,], file=subject_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
			}
		}
	}

    ####integrate and output results
    INTEG_LSL <- list();
    if ( INTEGRATION != "noIntegration" )
    {
        ###group the subjects according to groupID and then integrate the results
        ###with the function named as 'integrateLSOUT' in 'meta3dSubF.R'.
        subject_index <- 1:length(UNI_SUBJECTID);
        names(subject_index) <- UNI_SUBJECTID;
        if (length(GROUP_SUBJECTL))
        {
            for (k in 1:length(UNI_GROUPID))
            {
                ###get index of subjects under one group
                subject_lsID <- as.character(GROUP_SUBJECTL[[k]]);
                subject_lsindex <- as.numeric(subject_index[subject_lsID]);
                ###store analysis results from all subjects of the same group in one list
                i <- 1;
                subjects_integL <- list();
                for (n in subject_lsindex)
                {
                    subjects_integL[[i]] <- SUBJECT_LSL[[n]];
                    i <- i + 1;
                }
                ###integrate analysis results under the same group
                INTEG_LSL[[k]] <- integrateLSOUT(lsoutL=subjects_integL, adjustV=ADPHAL, weightS=weightedMethod);
            }
        }  else  {
            ###take all subjects as one group and integrate all subjects
            INTEG_LSL[[1]] <- integrateLSOUT(lsoutL=SUBJECT_LSL, adjustV=ADPHAL, weightS=weightedMethod);
        }

        ###output integrated results group by group-'INTEG_LSL'
        for ( m in 1:length(INTEG_LSL) )
        {
            ###output file name
            if (length(GROUP_SUBJECTL))  {
                group_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "meta3dGroupID_", UNI_GROUPID[m], "_", filename, sep="");
            }  else  {
                group_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "meta3dGroupID_AllSubjects_", filename, sep="");
            }
            ###get integrated results
            group_outputD <- INTEG_LSL[[m]];
            ###replace ID order with ID name
            group_nameorder <- as.character(group_outputD[,"CycID"]);
            group_outputD[,"CycID"] <- EXP_IDNAME[group_nameorder];
            ###output results
            if (!SIG_row)  {
                write.table(group_outputD, file=group_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
            }  else  {
                group_outputD[, "meta3d_BH.Q"] <- group_outputD[, "meta3d_Pvalue"];
                write.table(group_outputD[1,], file=group_outname, row.names=FALSE, sep=FILE_SEP, quote=FILE_QUOTE2, dec=FILE_DEC);
            }
        }
    }

    ####output the analysis note
    cat(paste("DONE! The analysis about '", datafile, "' and '",
            designfile, "'" , " has been finished.\n", sep=""));
    print( c("Time used:",(proc.time()-run_start)) );
    cat("\n\n");
}
