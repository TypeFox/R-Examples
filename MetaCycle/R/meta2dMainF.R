##' Detect rhythmic signals from time-series datasets with multiple methods
##'
##' This is a function that incorporates ARSER, JTK_CYCLE and Lomb-Scargle to
##'   detect rhythmic signals from time-series datasets.
##'
##' \href{https://github.com/cauyrd/ARSER}{ARSER}(Yang, 2010),
##'   \href{http://openwetware.org/wiki/HughesLab:JTK_Cycle}{JTK_CYCLE}(
##'   Hughes, 2010), and
##'   \href{http://research.stowers-institute.org/efg/2005/LombScargle/}{
##'   Lomb-Scargle}(Glynn, 2006) are three popular methods of detecting
##'   rhythmic signals. \code{ARS} can not analyze unevenly sampled datasets,
##'   or evenly sampled datasets but with missing values, or with replicate
##'   samples, or with non-integer sampling interval. \code{JTK} is not
##'   suitable to analyze unevenly sampled datasets or evenly sampled datasets
##'   but with non-integer sampling interval. If set \code{analysisStrategy}
##'   as \code{"auto"}(default), \code{meta2d} will automatically select
##'   proper method from \code{cycMethod} for each input dataset. If the user
##'   clearly know that the dataset could be analyzed by each method defined
##'   by \code{cycMethod} and do not hope to output integrated values,
##'   \code{analysisStrategy} can be set as \code{"selfUSE"}.
##'
##'   \code{ARS} used here is translated from its python version which always
##'   uses \code{"yule-walker"}, \code{"burg"}, and \code{"mle"} methods(see
##'   \code{\link[stats]{ar}}) to fit autoregressive models to time-series
##'   data. Fitting by \code{"mle"} will be very slow for datasets
##'   with many time points. If \code{ARSmle = "auto"} is used,
##'   \code{meta2d} will only include \code{"mle"} when number of time points
##'   is smaller than 24. In addition, one evaluation work(Wu, 2014) indicates
##'   that \code{ARS} shows relative high false positive rate in analyzing
##'   high-resolution datasets (1h/2days and 2h/2days). \code{JTK}(version 3)
##'   used here is the latest version, which improves its p-value calculation
##'   in analyzing datasets with missing values.
##'
##'   The power of detecting rhythmic signals for an algorithm is associated
##'   with the nature of data and interested periodic pattern(Deckard, 2013),
##'   which indicates that integrating analysis results from multiple methods
##'   may be helpful to rhythmic detection. For integrating p-values,
##'   Bonferroni correction(\code{"bonferroni"}) and Fisher's method(
##'   \code{"fisher"}) (Fisher, 1925; implementation code from \pkg{MADAM})
##'   could be selected, and \code{"bonferroni"} is usually more conservative
##'   than \code{"fisher"}. The integrated period is arithmetic mean of
##'   multiple periods. For integrating phase, \code{meta2d} takes use of
##'   \href{https://en.wikipedia.org/wiki/Mean_of_circular_quantities}{
##'   mean of circular quantities}. Integrated period and phase is further
##'   used to calculate the baseline value and amplitude through fitting a
##'   constructed periodic model.
##'
##'   Phases given by \code{JTK} and \code{LS} need to be adjusted with their
##'   predicted period (\code{adjustedPhase = "predictedPer"}) before
##'   integration. If \code{adjustedPhas = "notAdjusted"} is selected, no
##'   integrated phase will be calculated. If set \code{weightedPerPha} as
##'   \code{TRUE}, weighted scores will be used in averaging periods and
##'   phases. Weighted scores for one method are based on all its reported
##'   p-values, which means a weighted score assigned to any one profile will
##'   be affected by all other profiles. It is always a problem of averaging
##'   phases with quite different period lengths(eg. averaging two phases
##'   with 16-hours' and 30-hours' period length). Currently, setting
##'   \code{minper}, \code{maxper} and \code{ARSdefaultPer} to a same value
##'   may be the only way of completely eliminating such problem.
##'
##'   This function is originally aimed to analyze large scale periodic data(
##'   eg. circadian transcriptome data) without individual information.
##'   Please pay attention to data format of input file(see \code{Examples}
##'   part). Except the first column and first row, others are time-series
##'   experimental values(setting missing values as \code{NA}).
##'
##' @param infile a character string. The name of input file containing
##'   time-series data.
##' @param outdir a character string. The name of directory used to store
##'   output files.
##' @param filestyle a character vector(length 1 or 3). The data format of
##'   input file, must be \code{"txt"}, or \code{"csv"}, or a character vector
##'   containing field separator character(\code{sep}), quoting character
##'   (\code{quote}), and the character used for decimal points(\code{dec},
##'   for details see \code{\link[utils]{read.table}}).
##' @param timepoints a numeric vector corresponding to sampling time points
##'   of input time-series data; if sampling time points are in the first line
##'   of input file, it could be set as a character sting-"Line1" or "line1".
##' @param minper a numeric value. The minimum period length of interested
##'   rhythms. The default is \code{20} for circadian rhythms.
##' @param maxper a numeric value. The maximum period length of interested
##'   rhythms. The default is \code{28} for circadian rhythms.
##' @param cycMethod a character vector(length 1 or 2 or 3). User-defined
##'   methods for detecting rhythmic signals, must be selected as any one, any
##'   two or all three methods(default) from \code{"ARS"}(ARSER),
##'   \code{"JTK"}(JTK_CYCLE) and \code{"LS"}(Lomb-Scargle).
##' @param analysisStrategy a character string. The strategy used to select
##'   proper methods from \code{cycMethod} for analyzing input time-series
##'   data, must be \code{"auto"}(default), or \code{"selfUSE"}. See
##'   \code{Details} part for more information.
##' @param outputFile logical. If \code{TRUE}, analysis results will be wrote
##'   in the output files. If \code{FALSE}, analysis results will be returned
##'   as an R list.
##' @param outIntegration a character string. This parameter controls what
##'   kinds of analysis results will be outputted, must be one of \code{"both"}
##'   (default), \code{"onlyIntegration"}(only output integration file), or
##'   \code{"noIntegration"}(not output integration file).
##' @param adjustPhase a character string. The method used to adjust original
##'   phase calculated by each method in integration file, must be one of
##'   \code{"predictedPer"}(adjust phase with predicted period length) or
##'   \code{"notAdjusted"}(not adjust phase).
##' @param combinePvalue a character string. The method used to integrate
##'   multiple p-values, must be one of \code{"bonferroni"}(Bonferroni
##'   correction), or \code{"fisher"}(Fisher's method).
##' @param weightedPerPha logical. If \code{TRUE}, weighted scores based on
##'   p-value given by each method will be used to calculate the integrated
##'   period length and phase.
##' @param ARSmle a character string. The strategy of using MLE method in
##'   \code{\link[stats]{ar}} fit of \code{"ARS"}, must be one of
##'   \code{"auto"}(use MLE depending the number of time points), \code{"mle"}
##'   (always use MLE), or \code{"nomle"}(never use MLE).
##' @param ARSdefaultPer a numeric value. The expected period length of
##'   interested rhythm, which is a necessary parameter for \code{ARS}. The
##'   default is \code{24}(for circadian rhythms). Set it to another proper
##'   numeric value for other rhythms.
##' @param outRawData logical. If \code{TRUE}, raw time-series data will be
##'   added in the output files.
##' @param releaseNote logical. If \code{TRUE}, reminding or warning notes
##'   during the analysis will be released on the screen.
##' @param outSymbol a character string. A common prefix exists in the names of
##'   output files.
##' @return
##' \code{meta2d} will write analysis results in different files under
##'   \code{outdir} if set \code{outputFile = TRUE}. Files named with
##'   "ARSresult", "JTKresult" and "LSreult" store analysis results from
##'   \code{ARS}, \code{JTK} and \code{LS} respectively. The file named with
##'   "meta2d" is the integration file, and it stores integrated values in
##'   columns with a common name tag-"meta2d". The integration file also
##'   contains p-value, FDR value, period, phase(adjusted phase if
##'   \code{adjustedPhase = "predictedPer"}) and amplitude values calculated
##'   by each method.
##'   If \code{outputFile = FALSE} is selected, \code{meta2d} will return a
##'   list containing the following components:
##'   \tabular{rl}{
##'        ARS   \tab  analysis results from \code{ARS} method\cr
##'        JTK   \tab  analysis results from \code{JTK} method\cr
##'        LS    \tab  analysis results from \code{LS} method\cr
##'        meta  \tab  the integrated analysis results as mentioned above
##'   }
##' @references
##' Yang R. and  Su Z. (2010). Analyzing circadian expression data by
##'   harmonic regression based on autoregressive spectral estimation.
##'   \emph{Bioinformatics}, \bold{26(12)}, i168--i174.
##'
##' Hughes M. E., Hogenesch J. B. and Kornacker K. (2010). JTK_CYCLE: an
##'   efficient nonparametric algorithm for detecting rhythmic components in
##'   genome-scale data sets. \emph{Journal of Biological Rhythms},
##'   \bold{25(5)}, 372--380.
##'
##' Glynn E. F., Chen J. and Mushegian A. R. (2006). Detecting periodic
##'   patterns in unevenly spaced gene expression time series using
##'   Lomb-Scargle periodograms. \emph{Bioinformatics}, \bold{22(3)},
##'   310--316.
##'
##' Wu G., Zhu J., Yu J., Zhou L., Huang J. Z. and  Zhang Z. (2014). Evaluation
##'   of five methods for genome-wide circadian gene identification.
##'   \emph{Journal of Biological Rhythms}, \bold{29(4)}, 231--242.
##'
##' Deckard A., Anafi R. C., Hogenesch J. B., Haase S.B. and Harer J. (2013).
##'   Design and analysis of large-scale biological rhythm studies:
##'   a comparison of algorithms for detecting periodic signals in biological
##'   data. \emph{Bioinformatics}, \bold{29(24)}, 3174--3180.
##'
##' Fisher, R.A. (1925). \emph{Statistical methods for research workers}.
##'   Oliver and Boyd (Edinburgh).
##'
##' Kugler K. G., Mueller L.A. and Graber A. (2010). MADAM - an open source
##'   toolbox for meta-analysis. \emph{Source Code for Biology and Medicine},
##'   \bold{5}, 3.
##' @examples
##' # write 'cycMouseLiverProtein' into a 'txt' file
##' write.table(cycMouseLiverProtein, file="cycMouseLiverProtein.txt",
##'   sep="\t", quote=FALSE, row.names=FALSE)
##' # write 'cycSimu4h2d' and 'cycYeastCycle' into two 'csv' files
##' write.csv(cycSimu4h2d, file="cycSimu4h2d.csv", row.names=FALSE)
##' write.csv(cycYeastCycle, file="cycYeastCycle.csv", row.names=FALSE)
##'
##' # analyze 'cycMouseLiverProtein.txt' with JTK_CYCLE and Lomb-Scargle
##' meta2d(infile="cycMouseLiverProtein.txt", filestyle="txt",
##'   outdir="example", timepoints=rep(seq(0, 45, by=3), each=3),
##'   cycMethod=c("JTK","LS"), outIntegration="noIntegration")
##'
##' # analyze 'cycSimu4h2d.csv' with ARSER, JTK_CYCLE and Lomb-Scargle and
##' # output integration file with analysis results from each method
##' meta2d(infile="cycSimu4h2d.csv", filestyle="csv", outdir="example",
##'   timepoints="Line1")
##'
##' # analyze 'cycYeastCycle.csv' with ARSER, JTK_CYCLE and Lomb-Scargle to
##' # detect transcripts associated with cell cycle, and return analysis
##' # results instead of output them into files
##' cyc <- meta2d(infile="cycYeastCycle.csv",filestyle="csv",
##'   minper=80, maxper=96, timepoints=seq(2, 162, by=16),
##'   outputFile=FALSE, ARSdefaultPer=85, outRawData=TRUE)
##' head(cyc$ARS)
##' head(cyc$JTK)
##' head(cyc$LS)
##' head(cyc$meta)
##' @export
##' @importFrom gnm MPinv
##' @importFrom utils flush.console read.table write.table
##' @importFrom stats cov lm loess median optimize p.adjust pchisq pf pnorm
##'                   predict pt runif sd spec.ar var wilcox.test

meta2d <- function(infile, outdir="metaout", filestyle, timepoints, minper=20,
            maxper=28, cycMethod=c("ARS","JTK","LS"), analysisStrategy="auto",
            outputFile=TRUE, outIntegration="both", adjustPhase="predictedPer",
            combinePvalue="fisher",  weightedPerPha=FALSE,  ARSmle="auto",
            ARSdefaultPer=24,outRawData=FALSE,releaseNote=TRUE,outSymbol="")
{

    ####start time
    run_start=proc.time();

    ####extract 'infile', 'outdir', 'minper', 'maxper',
    ####set by users and store them as global variables
    INFILE <- infile;
    outdir2 <- unlist(strsplit(outdir,.Platform$file.sep));
    OUTDIR <- paste(outdir2,collapse=.Platform$file.sep);
    ###create the directory if it is not exist
    if (! file.exists(OUTDIR) )
    { dir.create(OUTDIR); }

    ####check whether input parameters are correct
    ###check 'minper' and 'maxper'
    ###the maximum and minimum period length should be positive number,
    ###and maximum period should not be smaller than minimum period length
    if ( (!is.numeric(minper)) | (minper <= 0 ) )
    {	stop("The 'minper' should be a positive numeric value.\n");	}
    if ( (!is.numeric(maxper)) | (maxper <= 0 ) )
    {	stop("The 'maxper' should be a positive numeric value.\n");	}
    if (minper > maxper)
    {	stop("The 'minper' should not be larger than 'maxper'.\n");	}
    MINPER <- minper;
    MAXPER <- maxper;
    ###check whether the 'timepoints' is a effective numeric vector, or
    ###the 'timepoints' information is stored in the first line of 'infile'
    ###the 'checkTimeF' is in 'meta2dSubF.R'
    ###NUMT is TRUE if 'timepoints' is a numeric vector
    NUMT <- checkTimeF(timepoints);

    ###check whether output integration results at the end of analysis
    if ( (outIntegration == "both") | (outIntegration == "onlyIntegration") | (outIntegration == "noIntegration") )  {
        INTEGRATION <- outIntegration;
    } else {
        stop(c("Please check the parameter of 'outIntegration', it should ",
               "be set as 'both', 'onlyIntegration' or 'noIntegration'.\n") );
    }
    ###check defined period length used to adjust phase
    if (adjustPhase == "predictedPer") {
        ADPHA <- "PERIOD";
    } else if (adjustPhase == "notAdjusted") {
        ADPHA <- 0;
    } else {
        stop("The 'adjustPhase' should be set as 'predictedPer', or 'notAdjusted'.\n");
    }
    ###check which method is used to integrate multiple P-values
    if ( (combinePvalue != "bonferroni")  &  (combinePvalue != "fisher") )
    {
        stop(c("The 'combinePvalue' should be set as 'bonferroni', or 'fisher'.\n") );
    }
    ###check 'weightedPerPha'
    if ( !is.logical(weightedPerPha) )
    {
        stop("The 'weightedPerPha' should be set as a logical value (TRUE or FALSE).\n");
    }
    ###check 'outRawData'
    if ( !is.logical(outRawData) )
    {
        stop("The 'outRawData' should be set as a logical value (TRUE or FALSE).\n");
    }

    ####extract field separator character(FILE_SEP), set of quoting
    ####characters(FILE_QUOTE), the character used for decimal points(FILE_DEC)
    ###the 'getFileSignF()' is in 'metaSubF.R'
    file_sign <- getFileSignF(filestyle);
    FILE_SEP <- file_sign$sep;
    FILE_QUOTE <- file_sign$quote;
    FILE_DEC <- file_sign$dec;
    FILE_QUOTE2 <- file_sign$quote2;

    ####extract user-defined method for analyzing expression profiles
    CIRM <- rep(FALSE,3);
    names(CIRM) <- c("ARS","JTK","LS");
    cycMethod <- unique(cycMethod);
    for (i in 1:length(cycMethod))
    {
        if ( (cycMethod[i] != "ARS") & (cycMethod[i] != "JTK") & (cycMethod[i] != "LS") )
        {
            stop(c("Please check 'cycMethod', one or multiple methods ",
               "could only be selected from 'ARS', 'JTK' and 'LS'.\n") );
        }
    }
    CIRM[cycMethod] <- TRUE;

    ####read time-series values and extract row and column's name
    EXP_dataframe <- ID_ORDER <- LIB_ID <- NULL;
    if (NUMT)
    {
        EXP_dataframe <- read.table(INFILE, header=TRUE, sep=FILE_SEP,
		                            quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
        ID_ORDER <- dimnames(EXP_dataframe)[[1]];
        LIB_ID <- dimnames(EXP_dataframe)[[2]];
    }  else  {
        EXP_dataframe <- read.table(INFILE, header=FALSE, sep=FILE_SEP,
		                            quote=FILE_QUOTE, dec=FILE_DEC, stringsAsFactors=FALSE);
        dimnames(EXP_dataframe)[[2]] <- as.character(EXP_dataframe[1,]);
        timepoints <- as.numeric(EXP_dataframe[1, -1]);
        if ( any(is.na(timepoints)) | any(is.nan(timepoints)) )
        {
            stop(c("Please replace non-numeric character or 'NA' with numeric ",
               "value (time points) in the first line of 'infile'.\n") );
        }
        if (min(timepoints) < 0 )
        {
            stop("The minimal time points value should not be a negative value.\n");
        }
        EXP_dataframe <- EXP_dataframe[-1,];
        dimnames(EXP_dataframe)[[1]] <- 1:nrow(EXP_dataframe);
        ID_ORDER <- dimnames(EXP_dataframe)[[1]];
        LIB_ID <- dimnames(EXP_dataframe)[[2]];
    }
    ###check whether correctly read-in the data
    if (ncol(EXP_dataframe) == 1)
    {
        stop(c("The 'infile' is stored in a one-column dataframe, ",
           "please check whether 'filestyle' is correctly set.\n") );
    }
    if (nrow(EXP_dataframe) == 0)
    {
        stop(c("The 'infile' contains only one line which is taken ",
           "as column names, please check the input file.\n") );
    }
    if ( length(timepoints) != ( ncol(EXP_dataframe) - 1) )
    {
        stop(c("Please check 'infile' or 'timepoints', the column number ",
           "of 'infile' should be one more than the length of time points. ",
           "See the example data in this package.\n") );
    }
    ###check whether contain only one row data
    SIG_row <- FALSE;
    if (nrow(EXP_dataframe) == 1)
    {
        EXP_dataframe <- rbind(EXP_dataframe, EXP_dataframe);
        dimnames(EXP_dataframe)[[1]] <- c("1","2");
        ID_ORDER <- c("1","2");
        SIG_row <- TRUE;
    }
    ###store names of each row in 'OUT_ID'
    OUT_ID <- as.character(EXP_dataframe[,1]);
    names(OUT_ID) <- ID_ORDER;

    ####sort and check time points
    ####re-arrange column order of input data according time order
    ###sort 'timepoints' value with increasing order
    names(timepoints) <- 1:length(timepoints);
    timepoints <- sort(timepoints);
    time_orderIndex <- as.numeric(names(timepoints));
    ###re-arrange data according to ordered time points
    EXPM <- data.matrix(EXP_dataframe[,-1]);
    dimnames(EXPM)[[1]] <- ID_ORDER;
    EXPM <- EXPM[,time_orderIndex];
    expm_ID <- LIB_ID[-1];
    expm_ID <- expm_ID[time_orderIndex];
    dimnames(EXPM)[[2]] <- expm_ID;
    ###transfer ordered data to 'EXP_dataframe'
    EXP_dataframe[ID_ORDER, 2:ncol(EXP_dataframe)] <- EXPM[ID_ORDER,];
    LIB_ID <- c(LIB_ID[1], expm_ID);
    dimnames(EXP_dataframe)[[2]] <- LIB_ID;

    ###check number of unique time points
    ###uni_timepoints should be ordered
    uni_timepoints <- unique(timepoints);
    if ( length(uni_timepoints) < 4 )
    {
        if (releaseNote) {
			    cat(c("Warning: the number of time points is too small, ",
			    "it may be not a good choice to analyze this expression file ",
			    "with 'ARS','JTK' or 'LS'. If possible, please carefully ",
			    "check the output results and try other useful methods.\n") );
		    }
    }
    ###get the initial time points for JTK
    START_TIME <- uni_timepoints[1];

    ####check 'ARSmle' and 'ARSdefaultPer' if 'ARS' is selected
    if (CIRM["ARS"])
    {
		    ###the 'checkARSdefaultPerF()' is in 'meta2dSubF.R'
		    ARS_PER <- checkARSdefaulterPerF(ARSdefaultPer, minper=MINPER, maxper=MAXPER);
		    ###If ARSmle == "auto", when total number of time points is >= 24, 'mle'
		    ###will not be used for estimating AR coefficients('mle' is too slow)
		    if ( (ARSmle == "auto") & ( length(uni_timepoints) < 24) ) {
            ARS_MET <- c("yule-walker","mle","burg");
		    } else if ( (ARSmle == "auto") & (length(uni_timepoints) >= 24) ) {
			      ARS_MET <- c("yule-walker","burg");
		    } else if (ARSmle == "mle") {
			      ARS_MET <- c("yule-walker","mle","burg");
		    } else if (ARSmle == "nomle") {
			      ARS_MET <- c("yule-walker","burg");
		    } else {
			      stop(c("Please check the parameter of 'ARSmle', it ",
			       "should be set as 'auto', 'mle' or 'nomle'.\n") );
		    }
	  }

    ####extract key features of input data, including with/without non-integer interval,
    ####even/uneven sampling, with/without missing values, with/without replicates
    MISSING_VALUE <- WITH_REPLICATE <- non_integerInterval <- uneven_interval <- FALSE;
    if ( !all( round(diff(uni_timepoints)) == diff(uni_timepoints) ) )
    {	non_integerInterval <- TRUE;	}
    if ( length( unique(diff(uni_timepoints)) ) > 1 )
    {	uneven_interval <- TRUE;	}
    if ( (!all( !is.na(EXPM) )) | (!all( !is.nan(EXPM) )) )
    {	MISSING_VALUE <- TRUE;	}
    if ( length(timepoints) != length(uni_timepoints) )
    {	WITH_REPLICATE <- TRUE;	}

    ####set the output file name
    filename <- unlist(strsplit(INFILE,.Platform$file.sep,fixed=TRUE));
    filename <- filename[length(filename)];
    if ( grepl("\\", filename, fixed=TRUE) )
    {
        filename <- unlist(strsplit(INFILE,"\\",fixed=TRUE));
        filename <- filename[length(filename)];
    }  else if ( grepl("//", filename, fixed=TRUE) )  {
        filename <- unlist(strsplit(INFILE,"//",fixed=TRUE));
        filename <- filename[length(filename)];
    }
    if (length(outSymbol) > 1)
    {
        outSymbol <- as.character(outSymbol);
        outSymbol <- paste(outSymbol, collapse="");
    }
    ars_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "ARSresult_", filename,sep="");
    jtk_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "JTKresult_", filename,sep="");
    ls_outname <- paste(OUTDIR, .Platform$file.sep,  outSymbol, "LSresult_", filename,sep="");
    outfile_name <- c(ars_outname, jtk_outname, ls_outname);
    names(outfile_name) <- c("ARS","JTK","LS");
    outfile_tag <- rep(0,3);
    names(outfile_tag) <- c("ARS","JTK","LS");
    integ_outname <- paste(OUTDIR, .Platform$file.sep, outSymbol, "meta2d_", filename,sep="");

    ####select proper method to analyze profiles depending on sampling pattern
    ####ARS (even sampling & without non-integer intervals & without missing values & without replicates)
    ####JTK (even sampling & without non-integer intervals)
    ####LS is not restricted in analyzing profiles in current design
    ARS_OUTM <- JTK_OUTM <- LS_OUTM <- "None_output";
    if (analysisStrategy == "auto") {
        if ( (non_integerInterval) | (uneven_interval) ) {
            if (CIRM["LS"])
            {
                LS_OUTM <- runLS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, releaseNote);
                outfile_tag["LS"] <- 1;
            }
        } else if ( (MISSING_VALUE) | (WITH_REPLICATE) ) {
            if (CIRM["JTK"])
            {
                JTK_OUTM <- runJTK(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, releaseNote);
                outfile_tag["JTK"] <- 1;
            }
            if (CIRM["LS"])
            {
                LS_OUTM <- runLS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, releaseNote);
                outfile_tag["LS"] <- 1;
            }
        } else if ( (!non_integerInterval)&(!uneven_interval)&(!MISSING_VALUE)&(!WITH_REPLICATE) ) {
            if (CIRM["ARS"])
            {
                ARS_OUTM <- runARS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER,
				                   arsper=ARS_PER, arsmet=ARS_MET, releaseNote);
                outfile_tag["ARS"] <- 1;
            }
            if (CIRM["JTK"])
            {
                JTK_OUTM <- runJTK(EXP_dataframe, timepoints, minper=MINPER,maxper=MAXPER, releaseNote);
                outfile_tag["JTK"] <- 1;
            }
            if (CIRM["LS"])
            {
                LS_OUTM <- runLS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER, releaseNote);
                outfile_tag["LS"] <- 1;
            }
        } else {
            stop("Please contact the author of this package for this bug.\n");
        }
    } else if (analysisStrategy == "selfUSE") {
        if (releaseNote) {
			     cat(c("Warning: 'analysisStrategy' is set as 'selfUSE', please be sure ",
            "that the 'infile' could be analyzed by all selected methods in ",
            "'cycMethod'. If not, please set 'analysisStrategy' as 'auto'.\n") );
        }
        if (CIRM["ARS"])
        {
            ARS_OUTM <- runARS(EXP_dataframe, timepoints, minper=MINPER, maxper=MAXPER,
                               arsper=ARS_PER, arsmet=ARS_MET, releaseNote);
            outfile_tag["ARS"] <- 1;
        }
        if (CIRM["JTK"])
        {
            JTK_OUTM <- runJTK(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER,releaseNote);
            outfile_tag["JTK"] <- 1;
        }
        if (CIRM["LS"])
        {
            LS_OUTM <- runLS(EXP_dataframe,timepoints,minper=MINPER,maxper=MAXPER,releaseNote);
            outfile_tag["LS"] <- 1;
        }
    } else {
        stop(c("Please check the parameter of 'analysisStrategy', ",
           "it should be set as 'auto' or 'selfUSE'.\n") );
    }

    ####output analysis result from each selected method
    integration_header <- "CycID";
    tabletypeL <- list("sep"=FILE_SEP, "quote"=FILE_QUOTE2, "dec"=FILE_DEC);
    ###common column names of each method
    sig_header <- c("pvalue","BH.Q","period","adjphase","amplitude");
    if (INTEGRATION != "onlyIntegration")
    {
        ###output ARS analysis results
        if (outfile_tag["ARS"])
        {
            outSigResultF(outM=ARS_OUTM, SIG_row=SIG_row, outRawData=outRawData,
			   EXPM=EXPM,ID_ORDER=ID_ORDER, methodName="ARS", outfile_name=outfile_name,
			   features=c("fdr_BH", "pvalue"), tabletypeL=tabletypeL, outputFile);
            integration_header <- c(integration_header, paste("ARS", sig_header, sep="_"));
        }
        ###output JTK analysis results
        if (outfile_tag["JTK"])
        {
            outSigResultF(outM=JTK_OUTM, SIG_row=SIG_row, outRawData=outRawData,
               EXPM=EXPM, ID_ORDER=ID_ORDER, methodName="JTK", outfile_name=outfile_name,
               features=c("BH.Q", "ADJ.P"),tabletypeL=tabletypeL, outputFile);
            integration_header <- c(integration_header, paste("JTK", sig_header, sep="_"));
        }
        ###output LS analysis results
        if (outfile_tag["LS"])
        {
            outSigResultF(outM=LS_OUTM, SIG_row=SIG_row, outRawData=outRawData,
               EXPM=EXPM, ID_ORDER=ID_ORDER, methodName="LS", outfile_name=outfile_name,
               features=c("BH.Q", "p"), tabletypeL=tabletypeL, outputFile);
            integration_header <- c(integration_header, paste("LS", sig_header, sep="_"));
        }
    }

    ####integration step and output integrated results
    AMPTIM <- timepoints;
    ###do not set 'integration_outM' as NULL here
    integration_outM <- "None_output";
    outLIST <- list("ARS"=ARS_OUTM,"JTK"=JTK_OUTM,"LS"=LS_OUTM);
    if (  INTEGRATION != "noIntegration" )
    {
        ###integration step depending on the number of used methods
        integration_num <- length(outfile_tag[outfile_tag > 0]);
        if (integration_num > 1) {
            ###if 'infile' is analyzed by two or three methods
            ###put pvalue, phase and period from multiple methods into three matrix
            integration_outM <- pvaM <- phaM <- perM <- ID_ORDER;
            out_index <- as.numeric(which(outfile_tag == 1));
            out_index <- sort(out_index);
            for (index in out_index)
            {
                sep_outM <- outLIST[[index]];
                if (length(sep_outM) > 1) {
                    if (index == 1) {
                        ars_adjM <- adjPhaARS(sep_outM[ID_ORDER,],adjustV=ADPHA);
                        integration_outM <- cbind(integration_outM,sep_outM[ID_ORDER,c("pvalue","fdr_BH")],ars_adjM);
                        pvaM <- cbind(pvaM,sep_outM[ID_ORDER,"pvalue"]);
                        phaM <- cbind(phaM,ars_adjM[ID_ORDER,"phase"]);
                        perM <- cbind(perM,ars_adjM[ID_ORDER,"period"]);
                        if (INTEGRATION == "onlyIntegration")
                        {  integration_header <- c(integration_header, paste("ARS", sig_header, sep="_"));  }
                    } else if (index == 2) {
                        jtk_adjpha <- adjPhaJTK(sep_outM[ID_ORDER,], adjustV=ADPHA, START_TIME);
                        integration_outM <- cbind(integration_outM, sep_outM[ID_ORDER,c("ADJ.P","BH.Q","PER")],
						                          jtk_adjpha[ID_ORDER], sep_outM[ID_ORDER,"AMP"]);
                        pvaM <- cbind(pvaM, sep_outM[ID_ORDER,"ADJ.P"]);
                        phaM <- cbind(phaM, jtk_adjpha[ID_ORDER]);
                        perM <- cbind(perM, sep_outM[ID_ORDER,"PER"]);
                        if (INTEGRATION == "onlyIntegration")
                        {  integration_header <- c(integration_header, paste("JTK", sig_header, sep="_"));  }
                    } else if (index == 3) {
                        ls_adjpha <- adjPhaLS(sep_outM[ID_ORDER,],adjustV=ADPHA);
                        integration_outM <- cbind(integration_outM, sep_outM[ID_ORDER,c("p","BH.Q","Period")],
						                          ls_adjpha[ID_ORDER], sep_outM[ID_ORDER,"PhaseShiftHeight"]);
                        pvaM <- cbind(pvaM, sep_outM[ID_ORDER,"p"]);
                        phaM <- cbind(phaM, ls_adjpha[ID_ORDER]);
                        perM <- cbind(perM, sep_outM[ID_ORDER,"Period"]);
                        if (INTEGRATION == "onlyIntegration")
                        {  integration_header <- c(integration_header, paste("LS", sig_header, sep="_"));  }
                    }
                } else {
                    stop(c("No result from ", names(outfile_tag[index]),
                        " during integration. If 'analysisStrategy' was set as ",
                        "'selfUSE', please check whether ", names(outfile_tag[index]),
                        " is suitable for analyzing the input file. If it was set as ",
                        "'auto', please contact the author.\n") );
                }
            }
            ###get column names of integrated results
            if (adjustPhase == "notAdjusted")
            {
                ###if 'notAdjusted' is selected, only show integrated pvalues in the results
                integration_header <- c(integration_header, "meta2d_pvalue","meta2d_BH.Q");
            }  else  {
                integration_header <- c(integration_header, paste("meta2d", sig_header[1:3], sep="_"),
                                       "meta2d_phase", "meta2d_Base", "meta2d_AMP", "meta2d_rAMP");
            }
            ###integrate pvalue, period, phase, and re-calculate amplitude
            if (ncol(pvaM) > 1)
            {
                if (nrow(pvaM) >= 2) {
                    cirpvaM <- getCirOmicsPVA(pvaM, method=combinePvalue);
                    if (adjustPhase == "notAdjusted")
                    {
                        ###if 'notAdjusted' is selected, only integrate pvalues
                        integration_outM <- cbind(integration_outM, cirpvaM[ID_ORDER,]);
                    }  else  {
                        cirperphaM <- getCirOmicsPerPha(periodM=perM, phaseM=phaM,
						                                   pvalueM=pvaM, adjustV=ADPHA, weightedPerPha);
                        EXPMPERPHA <- cbind(EXPM[ID_ORDER,],cirperphaM[ID_ORDER,]);
                        ampM <- getCirOmicsAMP(exprperphaM=EXPMPERPHA, AMPTIM);
                        dimnames(ampM)[[1]] <- ID_ORDER;
                        integration_outM <- cbind(integration_outM, cirpvaM[ID_ORDER,],
						                                      cirperphaM[ID_ORDER,], ampM[ID_ORDER,]);
                    }
                } else {
                    stop(c("There is a bug in integration ",
                         "step, please contact the author.\n") );
                }
            }
        } else {
                ###the 'infile' is analyzed by only one method, ARS or JTK or LS
                ###the 'getSigIntegResultF()' is in 'meta2dSubF.R'
                sigIntegL <- getSigIntegResultF(sigoutL=outLIST, outfile_tag=outfile_tag, EXPM=EXPM,
                                   ID_ORDER=ID_ORDER, ADPHA=ADPHA, AMPTIM=AMPTIM, START_TIME=START_TIME,
                                   INTEGRATION=INTEGRATION, integration_outM=integration_outM, integration_header);
                integration_outM <- sigIntegL$integ;
                integration_header <- sigIntegL$header;
        }
        ###output integrated results
        if (!is.vector(integration_outM)) {
            if ( ncol(integration_outM) > 1 )
            {
                integration_outM <- cbind(OUT_ID[ID_ORDER],integration_outM[,-1]);
                dimnames(integration_outM)[[1]] <- ID_ORDER;
                dimnames(integration_outM)[[2]] <- integration_header;
                ###output integrated results if 'outputFile=TRUE'
                if (outputFile)
                {
                    ###the 'outIntegResultF()' is in 'meta2dSubF.R'
                    outIntegResultF(integration_outM=integration_outM, EXPM=EXPM,
                            ID_ORDER=ID_ORDER, integration_header=integration_header, SIG_row=SIG_row,
                            outRawData=outRawData, integ_outname=integ_outname, tabletypeL);
                }
            } else {
                cycMethod <- paste(cycMethod, collapse=",");
                stop(c("No integration results from ", cycMethod, ". There is ",
                     "a bug in the integration step, please contact the author.\n") );
            }
        }
    }

    ####output analysis note on the screen
    if (sum(outfile_tag) == 0)
    {
        cycMethod <- paste(cycMethod, collapse=",");
        stop(c("Warning: no output result from selected method ",
             cycMethod, ". Please check whether at least one is suitable ",
             "for analyzing the input file, and whether all parameters ",
             "are correctly set as suggest. If it is quite sure about ",
             "selected methods and parameters, please contact the author.\n"));
    } else {
        if (releaseNote)  {
            cat(c("DONE! The analysis about '",INFILE,"'"," has been finished.\n"));
            print( c("Time used:",(proc.time()-run_start)) );
            cat("\n\n");
        }
    }

    ####return analysis results if 'outputFile=FALSE'
    if (!outputFile)
    {
        meta2dL <- returnResultF(outLIST=outLIST, integration_outM=integration_outM, EXPM=EXPM,
		                 ID_ORDER=ID_ORDER, outfile_tag=outfile_tag, SIG_row=SIG_row, outRawData)
        return(meta2dL);
    }
}
