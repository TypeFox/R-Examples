# -----------------------------------------------
# Function: Assemble summary table for ezANOVA.
# Requires an ezANOVA output object as argument. 
# -----------------------------------------------
anova_out <- function(ezout,
                      print   = TRUE,
                      sph.cor = "GG",
                      mau.p   = 0.05,
                      etasq   = "partial",
                      dfsep   = ", ") {
  
  # ---------------------------------------------
  # (1) Check input arguments
  # ---------------------------------------------
  x <- ezout;
  # Check for inconsistent sphericity mehod
  if (toupper(sph.cor)!="GG" &
        toupper(sph.cor)!="HF" &
        toupper(sph.cor)!="NO" ) {
    sph.cor="no"
    print(paste("Warning: Unknown correction method specified!",
                " Reporting uncorrected p-values instead.",sep=""),quote=FALSE)  
  }
  
  # Check for inconsistent sphericity mehod
  if (etasq!="partial" &
        etasq!="generalized" ) {
    etasq="partial"
    print(paste("Warning: Unknown effect size specified!",
                " Reporting partial eta squared instead.",sep=""),quote=FALSE)  
  }
  
  # ---------------------------------------------
  # (2) Assemble ANOVA table
  # ---------------------------------------------
  # Construct table
  doeswork <- 1;
  if ("ANOVA" %in% names(x)) {
    # Check whether SSn and SSd are present in input data
    if ("SSn" %in% colnames(x$ANOVA) && "SSd" %in% colnames(x$ANOVA)) {
      outtable <- data.frame(
        Effect=x$ANOVA$Effect,
        MSE=x$ANOVA$SSd/x$ANOVA$DFd,
        df1=x$ANOVA$DFn,
        df2=x$ANOVA$DFd,
        F=format(round(x$ANOVA$F,2),nsmall=2),
        p=format(round(x$ANOVA$p,3),nsmall=3),
        petasq=format(round(x$ANOVA$SSn/
                              (x$ANOVA$SSn+x$ANOVA$SSd),2),nsmall=2),
        getasq=format(round(x$ANOVA$ges,2),nsmall=2)
      );
    } else {
      outtable <- "Couldn't find Sum of Squares in ezANOVA output. Please use the 'detailed=TRUE' option!";
      doeswork <- 0;
    }
  } else {
    outtable <- "N/A ... possibly wrong input?";
    doeswork <- 0;
  }
  
  # Do remaining operations only if main output could be created
  if (doeswork < 1) {
    print(outtable);
  } else { 
    
    # ---------------------------------------------
    # (3) Sphericity tests
    # ---------------------------------------------  
    if ("Mauchly's Test for Sphericity" %in% names(x)) {
      outspher <- data.frame(
        Effect=x$"Mauchly's Test for Sphericity"$Effect,
        p_Mauchly=format(round(x$"Mauchly's Test for Sphericity"$p,3),nsmall=3),
        GGEpsilon=format(round(x$"Sphericity Corrections"$GGe,3),nsmall=3),
        p_GG=format(round(x$"Sphericity Corrections"$"p[GG]",3),nsmall=3),
        HFEpsilon=format(round(x$"Sphericity Corrections"$HFe,3),nsmall=3),
        p_HF=format(round(x$"Sphericity Corrections"$"p[HF]",3),nsmall=3)
      );
    } else {
      outspher <- "N/A"
      sph.cor = "no"
    }
    
    # ---------------------------------------------
    # (4) Prepare formatted output
    # ---------------------------------------------
    if ("ANOVA" %in% names(x)) {
      # Adjust p values when sphericity is violated
      txttable <- outtable;
      ajdffcts <- list();
      # Check all effects listed in "Sphericity Corrections"
      if (toupper(sph.cor)!="NO") {
        for (isph in 1:length(x$"Sphericity Corrections"$Effect)) {
          # Get effects of interest and check corresponding p_Mauchly
          if (x$"Mauchly's Test for Sphericity"$p[isph] <= mau.p) {
            eoi <- x$"Sphericity Corrections"$Effect[isph]
            # Cycle through ANOVA table and check for effects
            for (iaov in 1:length(x$ANOVA$Effect)) {
              # Adjust p-value
              if (x$ANOVA[iaov,1]==eoi) {
                if (toupper(sph.cor)=="GG") {
                  pmaucorr <- format(round(x$"Sphericity Corrections"$"p[GG]"[isph],3),nsmall=3);
                  levels(txttable$p) <- c(levels(txttable$p), pmaucorr)
                  txttable[iaov,6]=pmaucorr;
                } else if (toupper(sph.cor)=="HF") {
                  pmaucorr <- format(round(x$"Sphericity Corrections"$"p[HF]"[isph],3),nsmall=3);
                  levels(txttable$p) <- c(levels(txttable$p), pmaucorr)
                  txttable[iaov,6]=pmaucorr
                }
                ajdffcts <- c(ajdffcts,eoi);
              }
            }
          } # End: if p_Mauchly < p_crit
        }
        # Construct note
        if (length(ajdffcts) == 0) {
          note <- paste("No adjustments necessary (all p_Mauchly > ", mau.p,
                        ").",sep="")
        } else {
          note <- paste("p-values for the following effects were ", sph.cor,
                        "-adjusted (p_Mauchly <= ", mau.p, "): ",
                        paste(paste(ajdffcts,collapse="; ",sep=""),".",sep=""),sep="");
        }
        # Else/End: Check if sph.cor != "NO"
      } else {
        if (outspher!="N/A") {
          note <- "Reporting unadjusted p-values."
        }
      } 
      
      
      pcorr <- paste(", p = ", txttable$p, sep="")
      pcorr <- gsub("p = 1.000","p > .999", pcorr, fixed=TRUE)
      pcorr <- gsub("p = 0.000","p < .001", pcorr, fixed=TRUE)
      pcorr <- gsub("p = 0","p = ", pcorr, fixed=TRUE)
      
      # Get effect size
      if (etasq == "partial") {
        
        petasqcorr <- paste(", np2 = ", txttable$petasq, sep="")
        petasqcorr <- gsub("np2 = 1.00","np2 > .99", petasqcorr, fixed=TRUE)
        petasqcorr <- gsub("np2 = 0.00","np2 < .01", petasqcorr, fixed=TRUE)
        petasqcorr <- gsub("np2 = 0","np2 = ", petasqcorr, fixed=TRUE)
        
        outtext <- data.frame(
          Effect=x$ANOVA$Effect,
          Text=paste("F(", txttable$df1, "," ,txttable$df2, ") = ", txttable$F, 
                     pcorr, petasqcorr,sep=""));
        
      } else {
        
        getasqcorr <- paste(", ng2 = ", txttable$getasq, sep="")
        getasqcorr <- gsub("ng2 = 1.00","ng2 > .99", getasqcorr, fixed=TRUE)
        getasqcorr <- gsub("ng2 = 0.00","ng2 < .01", getasqcorr, fixed=TRUE)
        getasqcorr <- gsub("ng2 = 0","ng2 = ", getasqcorr, fixed=TRUE)
        
        outtext <- data.frame(
          Effect=x$ANOVA$Effect,
          Text=paste("F(", txttable$df1, dfsep, txttable$df2, ") = ", txttable$F, 
                     pcorr, getasqcorr,sep=""));
        
      }
      
    } else {
      outtext <- "N/A"
    }
    
    # ---------------------------------------------
    # (5) Combine and display ANOVA results
    # ---------------------------------------------  
    x <- list("--- ANOVA RESULTS     ------------------------------------" = outtable,
              "--- SPHERICITY TESTS  ------------------------------------" = outspher,
              "--- FORMATTED RESULTS ------------------------------------" = outtext);
    if (exists("note")) {
      x = c(x,"NOTE:"=note); 
    }
    if (print==TRUE) {
      print(x);
    } else {
      x;
    }
  } # Do only if doeswork > 0
}