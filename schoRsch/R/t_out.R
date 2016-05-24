# -----------------------------------------------
# Function: Assemble summary table for t-tests.
# Requires a t.test output object as argument. 
# -----------------------------------------------
t_out <- function(toutput,
                  n.equal = TRUE,
                  welch.df.exact = TRUE,
                  welch.n = NA,
                  d.corr = TRUE,
                  print = TRUE) {
  
  # ---------------------------------------------
  # (1) Assemble summary table
  # ---------------------------------------------
  # Compute Cohen's d
  if (length(grep("Two",toutput$method))>0) {
    # Two samples, equal sample sizes
    if (identical(TRUE,n.equal)) {
      if (length(grep("Welch",toutput$method))>0) {
        d = toutput$statistic * sqrt(4/(welch.n));
      } else {
        d = toutput$statistic * sqrt(4/(toutput$parameter+2));
      }
    # Two samples, unequal sample sizes
    } else if (identical(FALSE,n.equal)) {
      stop("Please enter sample sizes, e.g., n.equal = c(12,8).",
           call. = TRUE)
    } else {
      if ((sum(n.equal)==toutput$parameter+2) && (length(n.equal==2))) {
        d = toutput$statistic * sqrt(1/n.equal[1] + 1/n.equal[2])
      } else {
        stop("Sample sizes are inconsistent with degrees of freedom.",
             call. = TRUE)       
      }
    }
  # Paired samples, including correction
  } else if (d.corr==TRUE) {
    d = toutput$statistic / sqrt(toutput$parameter+1) * sqrt(2);
  # Paired samples, uncorrected
  } else {
    d = toutput$statistic / sqrt(toutput$parameter+1);
  }
  # Round Welch-adjusted dfs
  if (length(grep("Welch",toutput$method))>0) {
  if (welch.df.exact==TRUE) {
    toutput$parameter <- format(round(toutput$parameter,2),nsmall=2)
  } else {
    toutput$parameter <- (welch.n - 2)
  }}
  # Put everything together
  outtable <- data.frame(
    t=format(round(toutput$statistic,2),nsmall=2),
    df = toutput$parameter,
    p=format(round(toutput$p.value,3),nsmall=3),
    d=format(round(d,2),nsmall=2));
  
  # ---------------------------------------------
  # (2) Format output table
  # ---------------------------------------------
  #
  # Correct p value display
  pcorr <- paste(", p = ", outtable$p, sep="")
  pcorr <- gsub("p = 1.000","p > .999", pcorr, fixed=TRUE)
  pcorr <- gsub("p = 0.000","p < .001", pcorr, fixed=TRUE)
  pcorr <- gsub("p = 0","p = ", pcorr, fixed=TRUE)
  
  outtext <- data.frame(
    Test=paste(toutput$method,":",sep=""),
    Results=paste("t(",outtable$df,") = " , outtable$t, pcorr,
             ", d = ",outtable$d,sep=""));
    
  if (print==TRUE) {
    print(outtext);
  } else {
    outtext; 
  }
}