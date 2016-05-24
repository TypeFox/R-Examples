### Authors of original code of fisher's method: Karl Kugler <karl@eigenlab.net>, and Laurin AJ Mueller and Armin Graber
### Citation: MADAM - An Open Source Toolbox for Meta-Analysis. Source Code for Biology and Medicine 2010, 5:3
### The minor modification here is removing the code associated with 'multicore' package
###======================================================================================================================================
## function to calculate fisher sum
## p: vector of p-values
fisher.sum <- function(p, zero.sub=0.00001, na.rm=FALSE){
  if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
    stop("You provided bad p-values")
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  p[p==0] <- zero.sub
  if(na.rm)
    p<- p[!is.na(p)]
  S= -2*sum(log(p))
  res <- data.frame(S=S, num.p=length(p))
  return(res)
}
## main function of combining p-values by performing Fisher's method
fisher.method <- function(pvals, method=c("fisher"), p.corr=c("bonferroni","BH","none"), zero.sub=0.00001, na.rm=FALSE){
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum, zero.sub=zero.sub, na.rm=na.rm)))
  
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S, df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(fisher.sums$p.value, "bonferroni"),
                              BH = p.adjust(fisher.sums$p.value, "BH"),
                              none = fisher.sums$p.value)
  return(fisher.sums)
}
###======================================================================================================================================
###other functions shared by 'meta2d' and 'meta3d'
##adjusting phase with period length, eg. transfer phase from '26h' to '2h', if period length is 24h.
subAdjPha <- function(subpha, subper, adjustV)
{
    ###adjusting phase with expected period length
    if (adjustV == "PERIOD") {
        subper <- as.numeric(subper);
    } else if ( ( is.numeric(adjustV) ) & (adjustV > 0) ) {
        subper <- rep(adjustV, length(subpha));
    } else {
        return(as.numeric(subpha));
    }
    subpha <- as.numeric(subpha);
    pha_fold <- floor(subpha/subper);
    subadj_pha <- subpha - pha_fold*subper;
    return(subadj_pha);
}
##calculate the mean phase of multiple phases from different methods using 'mean of circular quantities'
circularMean <- function (z, subper, zweit, meanper, subadj)
{
	if ( length(z) > 0 )  {
		if (sum(zweit) > 0)
		{
			if (is.numeric(subadj))
			{
				if (subadj > 0)
				{   
					subper <- rep(subadj, length(subper)); 
				} else if (!subadj) {      
					##subadj = 0 (indicating no adjustment of phase), return linear mean values
					zmean <- sum(z*zweit)/sum(zweit);
					return(zmean);
				} else {
					stop( c("There is unknown bug associated with 'circularMean()'. ",
						  "Please contact the author. Thanks.\n") );
				}
			}
			
			##convert phase values to polar coordinates
			##for 'predictedPer', 'subper' here is period length from each method
			zpolar <- z/subper*2*pi;                                          
			##another strategy may use the same 'meanper' value for different methods (subper=meanper), which seems less reasonable
			##eg. if 'meanper'=25, the max per for JTK is 24 (one cycle sampling), it seems more reasonable to 
			##use predicted period by JTK for transferring its phase values to polar values than using the 'meanper'. 
			
			##convert polar coordinates to cartesian coordinates 
			siny <- sum(sin(zpolar)*zweit);                                   
			cosx <- sum(cos(zpolar)*zweit);
			##siny and cosx are not divided by 'sum(zweit)'; 'y/x' will get the same result if divided or not divided by the same value
			
			##get mean of circular quantities 
			meanpolar <- atan2(siny, cosx);                                   
			
			##convert to the mean phase value
			meanpha <- meanpolar/(2*pi)*meanper;
			if (meanpha < 0)
			{ meanpha <- meanpha +  meanper;  } 
			return(meanpha);
		}  else  {
			return(NaN);
		}
    }  else {
        return(NaN);
    }
}
##extract the field separator character(FILE_SEP), the set of quoting
##characters(FILE_QUOTE), the character used for decimal points(FILE_DEC).
getFileSignF <- function(filestyle)
{
    file_quote2 <- FALSE;
    if (length(filestyle) == 1) {
        if (filestyle=="csv") {
            file_sep <- ",";
            file_quote <- "\"";
            file_quote2 <- TRUE;
            file_dec <- ".";
        } else if (filestyle=="txt") {
            file_sep <- "\t";
            file_quote <- "";
            file_dec <- ".";
        } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' could be set as 'txt' or 'csv').\n") );
        }
    } else if (length(filestyle) > 1) {
        if (length(filestyle) == 3) {
            file_sep <- filestyle[1];
            file_quote <- filestyle[2];
            file_quote2 <- TRUE;
            file_dec <- filestyle[3];
        } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' should be assigned a vector containing ",
                   "three characters, which corresponding to symbols used to ",
                   "separate columns, quote values and used for decimal points).\n") );
        }
    } else {
            stop(c("Please set 'filestyle' before running this function ",
                   "(the 'filesyle' could be set as 'txt' or 'csv', or could be ",
                   "assigned a vector containing three characters, ",
                   "which corresponding to symbols used to separate columns, ",
                   "quote values and used for decimal points).\n") );
    }
	return(list("sep"=file_sep, "quote"=file_quote, "quote2"=file_quote2, "dec"=file_dec))
}
###======================================================================================================================================
