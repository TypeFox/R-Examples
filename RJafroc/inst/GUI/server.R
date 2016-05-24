ReportTextForGUI <- function(dataset, method = "DBMH", fom = "wJAFROC", alpha = 0.05, covEstMethod = "Jackknife", nBoots = 200) {
  UNINITIALIZED <- -Inf
  if (method == "DBMH") {
    methodTxt <- "DBM-MRMC HILLIS SIGNIFICANCE TESTING"
    result <- DBMHAnalysis(dataset, fom, alpha)
  } else if (method == "ORH") {
    methodTxt <- "OBUCHOWSKI-ROCKETTE-HILLIS SIGNIFICANCE TESTING"
    result <- ORHAnalysis(dataset, fom, alpha, covEstMethod, nBoots)
  } 
  
  ciPercent <- 100 * (1 - alpha)
  reportTxt <- paste0("RJafroc SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR ",
                      "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, ",
                      "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE ",
                      "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER ",
                      "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, ",
                      "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS ",
                      "IN THE SOFTWARE.\n================================================================================\n")
  reportTxt <- paste(reportTxt, sprintf(paste("Package build stats:", packageDescription("RJafroc", fields = "Built"))), sep = "\n")
  dateTime <- paste0("Run date: ", base::format(Sys.time(), "%b %d %Y %a %X %Z"))
  reportTxt <- paste(reportTxt, sprintf(dateTime), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" FOM selected         :     %s", fom), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf("================================================================================\n"), sep = "\n")
  
  NL <- dataset$NL
  LL <- dataset$LL
  lesionNum <- dataset$lesionNum
  lesionID <- dataset$lesionID
  lesionWeight <- dataset$lesionWeight
  maxNL <- dim(NL)[4]
  dataType <- dataset$dataType
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  I <- length(modalityID)
  J <- length(readerID)
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2
  nLesionPerCase <- rowSums(lesionID != UNINITIALIZED)
  
  
  reportTxt <- paste(reportTxt, sprintf(" Significance testing method:  %s", methodTxt), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" Number of Readers          :  %d", J), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" Number of Treatments       :  %d", I), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" Number of Normal Cases     :  %d", K1), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" Number of Abnormal Cases   :  %d", K2), sep = "\n")
  reportTxt <- paste(reportTxt, sprintf(" Fraction of Normal Cases   :  %f", K1/K), sep = "\n")
  
  if (dataType == "FROC") {
    reportTxt <- paste(reportTxt, sprintf(" Min number of lesions per diseased case   :  %d", min(nLesionPerCase)), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Max number of lesions per diseased case   :  %d", max(nLesionPerCase)), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Mean number of lesions per diseased case  :  %f", mean(nLesionPerCase)), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Total number of lesions                   :  %d", sum(nLesionPerCase)), sep = "\n")
    
    nl <- NL[, , (K1 + 1):K, ]
    dim(nl) <- c(I, J, K2, maxNL)
    maxNLRating <- apply(nl, c(1, 2, 3), max)
    maxLLRating <- apply(LL, c(1, 2, 3), max)
    maxNLRating[which(maxNLRating == UNINITIALIZED)] <- -2000
    maxLLRating[which(maxLLRating == UNINITIALIZED)] <- -2000
    ILF <- sum(maxNLRating > maxLLRating) + 0.5 * sum(maxNLRating == maxLLRating)
    ILF <- ILF/I/J/K2
    reportTxt <- paste(reportTxt, sprintf(" Inc. Loc. Frac.          :  %f\n\n", ILF), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf("================================================================================\n"), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Avg. number of non-lesion localization marks per reader on non-diseased cases: %f", sum(NL[, , 1:K1, ] != UNINITIALIZED)/(I * J * K1)), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Avg. number of non-lesion localization marks per reader on diseased cases:  %f", sum(NL[, , (K1 + 1):K, ] != UNINITIALIZED)/(I * J * K2)), sep = "\n")
    reportTxt <- paste(reportTxt, sprintf(" Avg. number of lesion localization marks per reader :  %f\n", sum(LL != UNINITIALIZED)/(I * J * K2)), sep = "\n")
  }
  
  
  reportTxt <- paste(reportTxt, paste("================================================================================\n", 
                                      " ====================================================================", 
                                      " *****                        Overview                          *****", 
                                      " ====================================================================", 
                                      " Three analyses are presented: ",
                                      " (1) Analysis 1 treats both readers and cases as random samples",
                                      "     --results apply to the reader and case populations;", 
                                      " (2) Analysis 2 treats only cases as a random sample", 
                                      "     --results apply to the population of cases but only for the", 
                                      "     readers used in this study; and", 
                                      " (3) Analysis 3 treats only readers as a random sample", 
                                      "     --results apply to the population of readers but only for the",
                                      "     cases used in this study.\n", 
                                      " For all three analyses, the null hypothesis of equal treatments is", 
                                      sprintf(" tested in part (a), treatment difference %d%% confidence intervals",                                                                                                                                                                                                                                 ciPercent), sprintf(" are given in part (b), and treatment %d%% confidence intervals are", ciPercent), " given in part (c).  Parts (a) and (b) are based on the treatment x", " reader x case ANOVA while part (c) is based on the reader x case", 
                                      " ANOVA for the specified treatment; these ANOVA tables are displayed", 
                                      " before the analyses.  Different error terms are used as indicated", 
                                      " for parts (a), (b), and (c) according to whether readers and cases", 
                                      " are treated as fixed or random factors.  Note that the treatment", 
                                      " confidence intervals in part (c) are based only on the data for the", 
                                      " specified treatment, rather than the pooled data.  Treatment", 
                                      sprintf(" difference %d%% confidence intervals for each reader are presented", ciPercent), 
                                      " in part (d) of Analysis 2; each interval is based on the treatment", 
                                      " x case ANOVA table (not included) for the specified reader.\n", sep = "\n"), 
                     sep = "\n")
  
  reportTxt <- paste(reportTxt, paste(" ===========================================================================", 
                                      " *****                            Estimates                            *****", 
                                      " ===========================================================================\n", 
                                      "                        TREATMENT", sep = "\n"), sep = "\n")
  
  string <- "              "
  for (i in 1:I) {
    string <- paste0(string, "----------")
    if (i < I) {
      string <- paste0(string, "---")
    }
  }
  reportTxt <- paste(reportTxt, string, sep = "\n")
  string <- "  READER      "
  for (i in 1:I) {
    string <- paste0(string, sprintf("%-10.10s", dataset$modalityID[i]))
    if (i < I) {
      string <- paste0(string, "   ")
    }
  }
  reportTxt <- paste(reportTxt, string, sep = "\n")
  string <- "----------    "
  for (i in 1:I) {
    string <- paste0(string, "----------")
    if (i < I) {
      string <- paste0(string, "   ")
    }
  }
  reportTxt <- paste(reportTxt, string, sep = "\n")
  
  for (j in 1:J) {
    string <- sprintf("%-10.10s    ", dataset$readerID[j])
    for (i in 1:I) {
      string <- paste0(string, sprintf("%10.8f", result$fomArray[i, j]))
      if (i < I) {
        string <- paste0(string, "   ")
      }
    }
    reportTxt <- paste(reportTxt, string, sep = "\n")
  }
  reportTxt <- paste(reportTxt, "\n", sep = "\n")
  reportTxt <- paste(reportTxt, " TREATMENT MEANS (averaged across readers)", "----------    -----------------------------", sep = "\n")
  for (i in 1:I) {
    string <- paste0(sprintf("%-10.10s    %10.8f", dataset$modalityID[i], mean(result$fomArray[i, ])))
    reportTxt <- paste(reportTxt, string, sep = "\n")
  }
  reportTxt <- paste(reportTxt, "\n\n", sep = "\n")
  reportTxt <- paste(reportTxt, " TREATMENT MEAN DIFFERENCES", "----------   ----------    -----------", sep = "\n")
  for (i in 1:I) {
    if (i < I) {
      for (ip in (i + 1):I) {
        reportTxt <- paste(reportTxt, sprintf("%-10.10s - %-10.10s    %10.8f", dataset$modalityID[i], dataset$modalityID[ip], mean(result$fomArray[i, ]) - mean(result$fomArray[ip, ])), sep = "\n")
      }
    }
  }
  reportTxt <- paste(reportTxt, "\n\n\n", sep = "\n")
  if (method == "DBMH") {
    if (J > 1) {
      reportTxt <- paste(reportTxt, 
                         " ===========================================================================",
                         " *****                          ANOVA Tables                           *****", 
                         " ===========================================================================\n", 
                         " TREATMENT X READER X CASE ANOVA\n", 
                         "Source            SS               DF             MS        ", 
                         "------   --------------------    ------   ------------------", sep = "\n")
      for (l in 1:7) {
        reportTxt <- paste(reportTxt, sprintf(" %5s   %20.8f    %6d   %18.8f", result$anovaY[l, 1], result$anovaY[l, 2], result$anovaY[l, 3], result$anovaY[l, 4]), sep = "\n")
      }
      reportTxt <- paste(reportTxt, sprintf(" %5s   %20.8f    %6d", result$anovaY[8, 1], result$anovaY[8, 2], result$anovaY[8, 3]), sep = "\n")
      reportTxt <- paste(reportTxt, "\n\n", sep = "\n")
      reportTxt <- paste(reportTxt, " TREATMENT X READER X CASE ANOVA", sep = "\n")
      reportTxt <- paste(reportTxt, "\n\n", sep = "\n")
      reportTxt <- paste(reportTxt, "                        Mean Squares", sep = "\n")
      string <- " Source     df   "
      for (i in 1:I) {
        string <- paste0(string, sprintf("%-10.10s", dataset$modalityID[i]))
        if (i < I) {
          string <- paste0(string, "   ")
        }
      }
      reportTxt <- paste(reportTxt, string, sep = "\n")
      string <- " ------    ---   "
      for (i in 1:I) {
        string <- paste0(string, "----------   ")
      }
      reportTxt <- paste(reportTxt, string, sep = "\n")
      for (l in 1:3) {
        string <- sprintf("     %2s %6d   ", result$anovaYi[l, 1], result$anovaYi[l, 2])
        for (i in 1:I) {
          string <- paste0(string, sprintf("%10.8f", result$anovaYi[l, i + 2]))
          if (i < I) {
            string <- paste0(string, "   ")
          }
        }
        reportTxt <- paste(reportTxt, string, sep = "\n")
      }
    }
    
    reportTxt <- paste(reportTxt, 
                       " ===========================================================================", 
                       " *****                  Variance Components Estimates                  *****", 
                       " ===========================================================================\n", 
                       " DBM variance component and covariance estimates\n", 
                       "     DBM Component             Estimate    ", 
                       " -----------------------  ----------------", 
                       sprintf(" Var(R)                  %16.8f", result$varComp$varComp[1]), 
                       sprintf(" Var(C)                  %16.8f", result$varComp$varComp[2]), 
                       sprintf(" Var(T*R)                %16.8f", result$varComp$varComp[3]), 
                       sprintf(" Var(T*C)                %16.8f", result$varComp$varComp[4]), 
                       sprintf(" Var(R*C)                %16.8f", result$varComp$varComp[5]), 
                       sprintf(" Var(Error)              %16.8f", result$varComp$varComp[6]), sep = "\n")
    
  } else {
    reportTxt <- paste(reportTxt, 
                       " ===========================================================================", 
                       " *****                  Variance Components Estimates                  *****", 
                       " ===========================================================================\n", 
                       " Obuchowski-Rockette variance component and covariance estimates\n", 
                       "     OR Component             Estimate    ", 
                       " -----------------------  ----------------", 
                       sprintf(" Var(R)                  %16.8f", result$varComp$varCov[1]), 
                       sprintf(" Var(T*R)                %16.8f", result$varComp$varCov[2]), 
                       sprintf(" COV1                    %16.8f", result$varComp$varCov[3]), 
                       sprintf(" COV2                    %16.8f", result$varComp$varCov[4]), 
                       sprintf(" COV3                    %16.8f", result$varComp$varCov[5]), 
                       sprintf(" Var(Error)              %16.8f", result$varComp$varCov[6]), sep = "\n")
  }
  
  smallestDispalyedPval <- 1e-04
  reportTxt <- paste(reportTxt, "\n", sep = "\n")
  if (J > 1) {
    reportTxt <- paste(reportTxt, 
                       " ===========================================================================", 
                       " *****           Analysis 1: Random Readers and Random Cases           *****", 
                       " ===========================================================================\n\n", 
                       " (Results apply to the population of readers and cases)\n\n", sep = "\n")
    
    reportTxt <- paste(reportTxt, sprintf("    a) Test for H0: Treatments have the same %s figure of merit.\n\n", fom), sep = "\n")
    reportTxt <- paste(reportTxt, 
                       " Source        DF    Mean Square      F value  Pr > F ", 
                       " ----------  ------  ---------------  -------  -------", sep = "\n")
    if (method == "DBMH") {
      if (result$pRRRC >= smallestDispalyedPval) {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$anovaY[1, 4], result$fRRRC, result$pRRRC), sep = "\n")
      } else {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$anovaY[1, 4], result$fRRRC, smallestDispalyedPval), sep = "\n")
      }
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfRRRC, result$anovaY[4, 4] + max(result$anovaY[5, 4] - result$anovaY[7, 4])), sep = "\n")
      reportTxt <- paste(reportTxt, " Error term: MS(TR) + max[MS(TC) - MS(TRC), 0]\n", sep = "\n")
    } else {
      if (result$pRRRC >= smallestDispalyedPval) {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$msT, result$fRRRC, result$pRRRC), sep = "\n")
      } else {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$msT, result$fRRRC, smallestDispalyedPval), sep = "\n")
      }
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfRRRC, result$msTR + max(J * (result$varComp[3, 2] - result$varComp[4, 2]), 0)), sep = "\n")
      reportTxt <- paste(reportTxt, " Error term: MS(TR) + J * max[Cov2 - Cov3, 0]\n", sep = "\n")
    }
    
    if (result$pRRRC < alpha) {
      reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not equal,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfRRRC, result$fRRRC, result$pRRRC), 
                         sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not significantly different,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfRRRC, result$fRRRC, result$pRRRC), 
                         sep = "\n")
    }
    reportTxt <- paste(reportTxt, sprintf("    b) %d%% confidence intervals for treatment differences\n", ciPercent), sep = "\n")
    reportTxt <- paste(reportTxt, 
                       sprintf("       Treatment         Estimate   StdErr      DF      t     Pr > t          %d%% CI      ", ciPercent), 
                       "----------   ----------  --------  --------  -------  ------  -------  -------------------", 
                       sep = "\n")
    ii <- 1
    for (i in 1:I) {
      if (i < I) {
        for (ip in (i + 1):I) {
          reportTxt <- paste(reportTxt, sprintf("%-10.10s - %-10.10s  %8.5f  %8.5f  %7.2f  %6.2f  %7.4f  %8.5f , %8.5f\n", dataset$modalityID[i], dataset$modalityID[ip], result$ciDiffTrtRRRC[ii, 2], result$ciDiffTrtRRRC[ii, 3], 
                                                result$ciDiffTrtRRRC[ii, 4], result$ciDiffTrtRRRC[ii, 5], result$ciDiffTrtRRRC[ii, 6], result$ciDiffTrtRRRC[ii, 7], result$ciDiffTrtRRRC[ii, 8]), sep = "\n")
          ii <- ii + 1
        }
      }
    }
    reportTxt <- paste(reportTxt, "\n", sep = "\n")
    if (I == 2) {
      reportTxt <- paste(reportTxt, " H0: the two treatments are equal.", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, 
                         sprintf(" * H0: the %d treatments are equal.  To control the overall ", I), 
                         " type I error rate at .05, we conclude that treatment differences", 
                         " with p < .05 are significant only if the global test in ", 
                         " (a) is also significant (i.e, p < .05).", sep = "\n")
    }
    if (method == "DBMH") {
      reportTxt <- paste(reportTxt, " Error term: MS(TR) + max[MS(TC) - MS(TRC), 0]\n\n", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, " Error term: MS(TR) + J * max[Cov2 - Cov3, 0]\n\n", sep = "\n")
    }
    reportTxt <- paste(reportTxt, 
                       sprintf("    c) %d%% treatment confidence intervals based on reader x case ANOVAs", ciPercent), 
                       "       for each treatment (each analysis is based only on data for the", 
                       "       specified treatment\n", 
                       sep = "\n")
    reportTxt <- paste(reportTxt, 
                       sprintf("  Treatment     Area      Std Error     DF     %d%% Confidence Interval ", ciPercent),
                       "  ----------  ----------  ----------  -------  -------------------------", sep = "\n")
    for (i in 1:I) {
      reportTxt <- paste(reportTxt, sprintf("  %-10.10s  %10.8f  %10.8f  %7.2f  (%10.8f , %10.8f)", 
                                            result$ciAvgRdrEachTrtRRRC[i, 1], result$ciAvgRdrEachTrtRRRC[i, 2], result$ciAvgRdrEachTrtRRRC[i, 3], result$ciAvgRdrEachTrtRRRC[i, 4], result$ciAvgRdrEachTrtRRRC[i, 5], result$ciAvgRdrEachTrtRRRC[i, 6]), 
                         sep = "\n")
    }
    if (method == "DBMH") {
      reportTxt <- paste(reportTxt, " Error term: MS(R) + max[MS(C) - MS(RC), 0]\n\n\n", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, "\n\n\n", sep = "\n")
    }
  }
  
  
  
  reportTxt <- paste(reportTxt, 
                     " ===========================================================================", 
                     " *****           Analysis 2: Fixed Readers and Random Cases            *****", 
                     " ===========================================================================\n\n", 
                     " (Results apply to the population of cases but only for the readers", 
                     " used in this study)\n\n", sep = "\n")
  
  reportTxt <- paste(reportTxt, sprintf("    a) Test for H0: Treatments have the same %s figure of merit.\n\n", fom), sep = "\n")
  reportTxt <- paste(reportTxt, 
                     " Source        DF    Mean Square      F value  Pr > F ", 
                     " ----------  ------  ---------------  -------  -------", sep = "\n")
  if (method == "DBMH") {
    if (result$pFRRC >= smallestDispalyedPval) {
      reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$anovaY[1, 4], result$fFRRC, result$pFRRC), sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$anovaY[1, 4], result$fFRRC, smallestDispalyedPval), sep = "\n")
    }
    reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfFRRC, result$anovaY[5, 4]), sep = "\n")
    reportTxt <- paste(reportTxt, "  Error term: MS(TC)\n", sep = "\n")
  } else {
    if (result$pFRRC >= smallestDispalyedPval) {
      reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$msT, result$fFRRC, result$pFRRC), sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$msT, result$fFRRC, smallestDispalyedPval), sep = "\n")
    }
    if (J > 1) {
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfFRRC, (result$varComp[1, 2] - result$varComp[2, 2] + (J - 1) * (result$varComp[3, 2] - result$varComp[4, 2]))), sep = "\n")
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1 + (J - 1) * ( Cov2 - Cov3 )\n", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfFRRC, (result$varComp[1, 2] - result$varComp[2, 2])), sep = "\n")
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1\n", sep = "\n")
    }
  }
  
  if (result$pFRRC < alpha) {
    reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not equal,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfFRRC, result$fFRRC, result$pFRRC), sep = "\n")
  } else {
    reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not significantly different,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfFRRC, result$fFRRC, result$pFRRC), 
                       sep = "\n")
  }
  
  reportTxt <- paste(reportTxt, sprintf("    b) %d%% confidence intervals for treatment differences\n", ciPercent), sep = "\n")
  reportTxt <- paste(reportTxt, 
                     sprintf("       Treatment         Estimate   StdErr      DF      t     Pr > t          %d%% CI      ", ciPercent), 
                     "----------   ----------  --------  --------  -------  ------  -------  -------------------", 
                     sep = "\n")
  ii <- 1
  for (i in 1:I) {
    if (i < I) {
      for (ip in (i + 1):I) {
        reportTxt <- paste(reportTxt, sprintf("%-10.10s - %-10.10s  %8.5f  %8.5f  %7.2f  %6.2f  %7.4f  %8.5f , %8.5f\n", dataset$modalityID[i], dataset$modalityID[ip], result$ciDiffTrtFRRC[ii, 2], result$ciDiffTrtFRRC[ii, 3], result$ciDiffTrtFRRC[ii, 
                                                                                                                                                                                                                                                       4], result$ciDiffTrtFRRC[ii, 5], result$ciDiffTrtFRRC[ii, 6], result$ciDiffTrtFRRC[ii, 7], result$ciDiffTrtFRRC[ii, 8]), sep = "\n")
        ii <- ii + 1
      }
    }
  }
  reportTxt <- paste(reportTxt, "\n", sep = "\n")
  if (I == 2) {
    reportTxt <- paste(reportTxt, " H0: the two treatments are equal.", sep = "\n")
  } else {
    reportTxt <- paste(reportTxt, 
                       sprintf(" * H0: the %d treatments are equal.  To control the overall ", I), 
                       " type I error rate at .05, we conclude that treatment differences", 
                       " with p < .05 are significant only if the global test in ", 
                       " (a) is also significant (i.e, p < .05).", sep = "\n")
  }
  if (method == "DBMH") {
    reportTxt <- paste(reportTxt, " Error term: MS(TC) \n\n", sep = "\n")
  } else {
    if (J > 1) {
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1 + (J - 1) * ( Cov2 - Cov3 )\n", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1\n", sep = "\n")
    }
  }
  
  reportTxt <- paste(reportTxt, 
                     sprintf("    c) %d%% treatment confidence intervals based on reader x case ANOVAs", ciPercent), 
                     "       for each treatment (each analysis is based only on data for the", 
                     "       specified treatment\n", 
                     sep = "\n")
  reportTxt <- paste(reportTxt, 
                     sprintf("  Treatment     Area      Std Error     DF     %d%% Confidence Interval ", ciPercent), 
                     "  ----------  ----------  ----------  -------  -------------------------", sep = "\n")
  for (i in 1:I) {
    reportTxt <- paste(reportTxt, sprintf("  %-10.10s  %10.8f  %10.8f  %7.2f  (%10.8f , %10.8f)", 
                                          result$ciAvgRdrEachTrtFRRC[i, 1], result$ciAvgRdrEachTrtFRRC[i, 2], result$ciAvgRdrEachTrtFRRC[i, 3], result$ciAvgRdrEachTrtFRRC[i, 4], 
                                          result$ciAvgRdrEachTrtFRRC[i, 5], result$ciAvgRdrEachTrtFRRC[i, 6]), sep = "\n")
  }
  if (method == "DBMH") {
    reportTxt <- paste(reportTxt, " Error term: MS(C) \n\n\n", sep = "\n")
  } else {
    if (J > 1) {
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1 + (J - 1) * ( Cov2 - Cov3 )\n", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, " Error term: Var - Cov1\n", sep = "\n")
    }
  }
  if (method == "DBMH") {
    reportTxt <- paste(reportTxt, " TREATMENT X CASE ANOVAs for each reader\n\n", sep = "\n")
    reportTxt <- paste(reportTxt, "                        Sum of Squares", sep = "\n")
    string <- " Source     df   "
    for (j in 1:J) string <- paste0(string, sprintf("%-11.11s   ", dataset$readerID[j]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- " ------    ---   "
    for (j in 1:J) string <- paste0(string, sprintf("-----------   ", dataset$readerID[j]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("      T %6d   ", I - 1)
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$ssAnovaEachRdr[1, j + 2]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("      C %6d   ", K - 1)
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$ssAnovaEachRdr[2, j + 2]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("     TC %6d   ", (I - 1) * (K - 1))
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$ssAnovaEachRdr[3, j + 2]))
    reportTxt <- paste(reportTxt, string, "\n\n", sep = "\n")
    reportTxt <- paste(reportTxt, "                        Mean Squares", sep = "\n")
    string <- " Source     df   "
    for (j in 1:J) string <- paste0(string, sprintf("%-11.11s   ", dataset$readerID[j]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- " ------    ---   "
    for (j in 1:J) string <- paste0(string, sprintf("-----------   ", dataset$readerID[j]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("      T %6d   ", I - 1)
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$msAnovaEachRdr[1, j + 2]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("      C %6d   ", K - 1)
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$msAnovaEachRdr[2, j + 2]))
    reportTxt <- paste(reportTxt, string, sep = "\n")
    string <- sprintf("     TC %6d   ", (I - 1) * (K - 1))
    for (j in 1:J) string <- paste0(string, sprintf("%11.7f   ", result$msAnovaEachRdr[3, j + 2]))
    reportTxt <- paste(reportTxt, string, "\n\n\n\n", sep = "\n")
  }
  
  
  reportTxt <- paste(reportTxt, "    d) Treatment-by-case ANOVA CIs for each reader ", sep = "\n")
  reportTxt <- paste(reportTxt, "       (each analysis is based only on data for the specified reader)\n", sep = "\n")
  reportTxt <- paste(reportTxt, 
                     sprintf("  Reader         Treatment        Estimate  StdErr       DF      t     Pr > t          %d%% CI      ", ciPercent),
                     "---------- ---------- ----------  --------  --------  -------  ------  -------  -------------------", 
                     sep = "\n")
  l <- 1
  for (j in 1:J) {
    for (i in 1:I) {
      if (i < I) {
        for (ip in (i + 1):I) {
          reportTxt <- paste(reportTxt, sprintf("%-10.10s %-10.10s-%-10.10s  %8.5f  %8.5f  %7.2f  %6.2f  %7.4f  %8.5f , %8.5f", 
                                                dataset$readerID[j], dataset$modalityID[i], dataset$modalityID[ip], result$ciDiffTrtEachRdr[l, 3], 
                                                result$ciDiffTrtEachRdr[l, 4], result$ciDiffTrtEachRdr[l, 5], result$ciDiffTrtEachRdr[l, 6], 
                                                result$ciDiffTrtEachRdr[l, 7], result$ciDiffTrtEachRdr[l, 8], result$ciDiffTrtEachRdr[l, 9]), 
                             sep = "\n")
          l <- l + 1
        }
      }
    }
  }
  if (method == "ORH") {
    string <- "\nReader  Var(Error)     Cov1   \n------  ----------  ----------"
    reportTxt <- paste(reportTxt, string, sep = "\n")
    for (j in 1:J) {
      reportTxt <- paste(reportTxt, sprintf("%-6.6s  %10.8s  %10.8s", result$varCovEachRdr[j, 1], result$varCovEachRdr[j, 2], result$varCovEachRdr[j, 3]), sep = "\n")
    }
  }
  
  reportTxt <- paste(reportTxt, "\n\n", sep = "\n")
  
  if (J > 1) {
    reportTxt <- paste(reportTxt, 
                       " ===========================================================================", 
                       " *****           Analysis 3: Random Readers and Fixed Cases            *****", 
                       " ===========================================================================", 
                       " (Results apply to the population of readers but only for the cases used in this study)\n\n", sep = "\n")
    
    reportTxt <- paste(reportTxt, sprintf("    a) Test for H0: Treatments have the same %s figure of merit.\n\n", fom), sep = "\n")
    reportTxt <- paste(reportTxt, 
                       " Source        DF    Mean Square      F value  Pr > F ", 
                       " ----------  ------  ---------------  -------  -------", sep = "\n")
    if (method == "DBMH") {
      if (result$pRRFC >= smallestDispalyedPval) {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$anovaY[1, 4], result$fRRFC, result$pRRFC), sep = "\n")
      } else {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$anovaY[1, 4], result$fRRFC, smallestDispalyedPval), sep = "\n")
      }
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfRRFC, result$anovaY[4, 4]), sep = "\n")
    } else {
      if (result$pRRFC >= smallestDispalyedPval) {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  %7.4f", I - 1, result$msT, result$fRRFC, result$pRRFC), sep = "\n")
      } else {
        reportTxt <- paste(reportTxt, sprintf(" Treatment   %6d  %15.8f  %7.2f  <%6.4f", I - 1, result$msT, result$fRRFC, smallestDispalyedPval), sep = "\n")
      }
      reportTxt <- paste(reportTxt, sprintf(" Error       %6.2f  %15.8f", result$ddfRRFC, result$msTR), sep = "\n")
    }
    reportTxt <- paste(reportTxt, " Error term: MS(TR)\n", sep = "\n")
    
    if (result$pRRFC < alpha) {
      reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not equal,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfRRFC, result$fRRFC, result$pRRFC), 
                         sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, sprintf(" Conclusion: The %s FOMs of treatments are not significantly different,\n             F(%d,%3.2f) = %3.2f, p = %6.4f.\n\n", fom, I - 1, result$ddfRRFC, result$fRRFC, result$pRRFC), 
                         sep = "\n")
    }
    reportTxt <- paste(reportTxt, sprintf("    b) %d%% confidence intervals for treatment differences\n", ciPercent), sep = "\n")
    reportTxt <- paste(reportTxt, 
                       sprintf("       Treatment         Estimate   StdErr      DF      t     Pr > t          %d%% CI      ", ciPercent), 
                       "----------   ----------  --------  --------  -------  ------  -------  -------------------", 
                       sep = "\n")
    ii <- 1
    for (i in 1:I) {
      if (i < I) {
        for (ip in (i + 1):I) {
          reportTxt <- paste(reportTxt, sprintf("%-10.10s - %-10.10s  %8.5f  %8.5f  %7.2f  %6.2f  %7.4f  %8.5f , %8.5f\n", 
                                                dataset$modalityID[i], dataset$modalityID[ip], result$ciDiffTrtRRFC[ii, 2], 
                                                result$ciDiffTrtRRFC[ii, 3], result$ciDiffTrtRRFC[ii, 4], result$ciDiffTrtRRFC[ii, 5], 
                                                result$ciDiffTrtRRFC[ii, 6], result$ciDiffTrtRRFC[ii, 7], result$ciDiffTrtRRFC[ii, 8]), 
                             sep = "\n")
          ii <- ii + 1
        }
      }
    }
    
    reportTxt <- paste(reportTxt, "\n", sep = "\n")
    if (I == 2) {
      reportTxt <- paste(reportTxt, " H0: the two treatments are equal.", sep = "\n")
    } else {
      reportTxt <- paste(reportTxt, 
                         sprintf(" * H0: the %d treatments are equal.  To control the overall ", I), 
                         " type I error rate at .05, we conclude that treatment differences", 
                         " with p < .05 are significant only if the global test in ", 
                         " (a) is also significant (i.e, p < .05).", sep = "\n")
    }
    reportTxt <- paste(reportTxt, "\n\n", sep = "\n")
    
    reportTxt <- paste(reportTxt, 
                       "    c) Reader-by-case ANOVAs for each treatment (each analysis is based only on data for the", 
                       "       specified treatment\n", sep = "\n")
    reportTxt <- paste(reportTxt, 
                       sprintf("  Treatment     Area      Std Error     DF     %d%% Confidence Interval ", ciPercent), 
                       "  ----------  ----------  ----------  -------  -------------------------", sep = "\n")
    for (i in 1:I) {
      reportTxt <- paste(reportTxt, sprintf("  %-10.10s  %10.8f  %10.8f  %7.2f  (%10.8f , %10.8f)", 
                                            result$ciAvgRdrEachTrtRRFC[i, 1], result$ciAvgRdrEachTrtRRFC[i, 2], 
                                            result$ciAvgRdrEachTrtRRFC[i, 3], result$ciAvgRdrEachTrtRRFC[i, 4], 
                                            result$ciAvgRdrEachTrtRRFC[i, 5], result$ciAvgRdrEachTrtRRFC[i, 6]), 
                         sep = "\n")
    }
  }
  return (reportTxt)
} 


shinyServer(function(input, output) {  
  #source("system.file(\"GUI\", \"ReportTextForGUI.R\", package = \"RJafroc\")")
  values <- reactiveValues()
  output$ui <- renderUI({
    if (is.null(input$dataFile)){
      wellPanel(
        p(paste0("RJafroc SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR ",
                 "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, ",
                 "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE ",
                 "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER ",
                 "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, ",
                 "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS ",
                 "IN THE SOFTWARE."), style = "font-size:12pt")
      )
    }else{            
      fileName <- paste0(input$dataFile$datapath, ".", file_ext(input$dataFile$name))
      file.rename(input$dataFile$datapath, fileName)
      if (file_ext(input$dataFile$name) %in% c("xls", "xlsx")){
        dataset <- ReadDataFile(fileName)  
      }else if (file_ext(input$dataFile$name) %in% c("txt", "csv", "lrc")){
        dataset <- ReadDataFile(fileName, format = "MRMC")  
      }else if (file_ext(input$dataFile$name) == "imrmc"){
        dataset <- ReadDataFile(fileName, format = "iMRMC")  
      }else{
        stop("Invalid data format.")
      }
      
      values$dataset <- dataset
      fomBtn <- list()
      if (dataset$dataType == "ROC"){
        fomBtn[[1]] <- c("Wilcoxon" = "Wilcoxon",
                         "HrSe" = "HrSe",
                         "HrSp" = "HrSp",
                         "SongA1" = "SongA1",
                         "SongA2" = "SongA2",
                         "MaxLLF" = "MaxLLF",
                         "MaxNLF" = "MaxNLF",
                         "MaxNLFAllCases" = "MaxNLFAllCases",
                         "ExpTrnsfmSp" = "ExpTrnsfmSp",
                         "JAFROC" = "JAFROC",
                         "Weighted JAFROC" = "wJAFROC",
                         "JAFROC1" = "JAFROC1",
                         "Weighted JAFROC1" = "wJAFROC1")  
        fomBtn <- list(
          radioButtons("fom", "Figure of Merit",
                       fomBtn[[1]],
                       inline = TRUE,
                       selected = "wJAFROC"
          )
        )
      }else if (dataset$dataType == "FROC"){
        fomBtn[[1]] <- c("HrAuc" = "HrAuc",
                         "HrSe" = "HrSe",
                         "HrSp" = "HrSp",
                         "SongA1" = "SongA1",
                         "SongA2" = "SongA2",
                         "MaxLLF" = "MaxLLF",
                         "MaxNLF" = "MaxNLF",
                         "MaxNLFAllCases" = "MaxNLFAllCases",
                         "ExpTrnsfmSp" = "ExpTrnsfmSp",
                         "JAFROC" = "JAFROC",
                         "Weighted JAFROC" = "wJAFROC",
                         "JAFROC1" = "JAFROC1",
                         "Weighted JAFROC1" = "wJAFROC1")  
        fomBtn <- list(
          radioButtons("fom", "Figure of Merit",
                       fomBtn[[1]],
                       inline = TRUE,
                       selected = "wJAFROC"
          )
        )
      }else if (dataset$dataType == "ROI"){
        fomBtn[[1]] <- c("ROI" = "ROI")
        fomBtn <- radioButtons("fom", "Figure of Merit",
                               fomBtn[[1]],
                               inline = TRUE
        )
      }     
      trtGroup <- as.list(1:length(dataset$modalityID))
      rdrGroup <- as.list(1:length(dataset$readerID))
      names(trtGroup) <- dataset$modalityID
      names(rdrGroup) <- dataset$readerID 
      values$trtGroup <- trtGroup
      values$rdrGroup <- rdrGroup
      fluidPage(
        tabsetPanel(type = "tabs",
                    tabPanel("Data Viewer", 
                             tabsetPanel(type = "tabs", 
                                         tabPanel("Truth", tableOutput("truthTable")), 
                                         tabPanel("NL", tableOutput("nlTable")), 
                                         tabPanel("LL", tableOutput("llTable")),
                                         tabPanel("Data Conversion", uiOutput("dataCnvt"))
                             )
                    ),
                    
                    tabPanel("Analysis",
                             div(sidebarLayout(
                               sidebarPanel(
                                 wellPanel(
                                   p(sprintf("Data type: %s", dataset$dataType)),
                                   p(sprintf("Number of modalities: %d", length(dataset$modalityID))),
                                   p(sprintf("Number of readers: %d", length(dataset$readerID))),
                                   p(sprintf("Number of normal cases: %d", dim(dataset$NL)[3])),
                                   p(sprintf("Number of abnormal cases: %d", dim(dataset$LL)[3]))
                                 ),
                                 radioButtons("mthd", "Analysis Methods:",
                                              c("DBMH" = "DBMH",
                                                "ORH" = "ORH"),
                                              inline = TRUE
                                 ),
                                 fomBtn,
                                 uiOutput("covEstBtn"),
                                 uiOutput("nBootsBtn"),
                                 textInput("alpha", HTML("Significance Level (&alpha;)"), value = 0.05),
                                 actionButton("analyzeBtn", "Analyze"),
                                 downloadButton("downloadReportBtn", "Save Report", class = NULL),
                                 width = 4
                               ),
                               mainPanel(
                                 tags$style(type='text/css', '#report {font-size: 8pt}'), 
                                 verbatimTextOutput("report"),
                                 width = 8                          
                               )
                             ),
                             style = 'width:1200px;'
                             )
                    ),
                    
                    tabPanel("Plotting",
                             div(sidebarLayout(
                               sidebarPanel(  
                                 uiOutput("plotboxes"),
                                 actionButton("nGroup", "Add a plotting group"),
                                 width = 4
                               ),
                               mainPanel(                                 
                                 uiOutput("plots")                                 
                               )
                             ),
                             style = 'width:1200px;'
                             )
                    ),
                    
                    tabPanel("Sample Size",
                             tabsetPanel(type = "tabs", 
                                         tabPanel("For input data file",
                                                  div(sidebarLayout(
                                                    sidebarPanel(
                                                      textInput("alphaDataPower", HTML("Significance Level (&alpha;)"), value = 0.05),
                                                      textInput("effectSizeDataPower", "Effect Size", value = 0.05),
                                                      textInput("desiredDataPower", "Desired Power", value = 0.8),
                                                      radioButtons("randomDataPower", "Random Option:",
                                                                   c("ALL" = "ALL",
                                                                     "READERS" = "READERS",
                                                                     "CASES" = "CASES"),
                                                                   inline = TRUE
                                                      ),
                                                      actionButton("calculateDataPowerBtn", "Calculate"),
                                                      width = 4
                                                    ),
                                                    mainPanel(
                                                      wellPanel(
                                                        tableOutput("powerTable")
                                                      ),
                                                      width = 8                          
                                                    )
                                                  ),
                                                  style = 'width:1100px;'
                                                  )
                                         ),
                                         
                                         tabPanel("Use variance components",
                                                  div(sidebarLayout(
                                                    sidebarPanel(
                                                      textInput("JSampleSize", "Number of Readers"),                                            
                                                      radioButtons("compType", "Analysis Methods:",
                                                                   c("DBMH" = "DBMH",
                                                                     "ORH" = "ORH"),
                                                                   inline = TRUE
                                                      ),
                                                      uiOutput("varComp"),
                                                      textInput("alphaSampleSize", HTML("Significance Level (&alpha;)"), value = 0.05),
                                                      textInput("effectSizeSampleSize", "Effect Size", value = 0.05),
                                                      textInput("desiredSampleSize", "Desired Power", value = 0.8),
                                                      radioButtons("randomSampleSize", "Random Option:",
                                                                   c("ALL" = "ALL",
                                                                     "READERS" = "READERS",
                                                                     "CASES" = "CASES"),
                                                                   inline = TRUE
                                                      ),
                                                      actionButton("calculateSampleSizeBtn", "Calculate"),
                                                      width = 4
                                                    ),
                                                    mainPanel(
                                                      verbatimTextOutput("sampleSize"),
                                                      width = 8                          
                                                    )
                                                  ),
                                                  style = 'width:1200px;'
                                                  )
                                         )
                             )
                    )
        )
      )
    }
  })
  
  output$dataCnvt <- renderUI({
    dataset <- values$dataset
    if (dataset$dataType == "FROC"){
      list(
        p(sprintf("Data type: %s", dataset$dataType)),
        checkboxInput("ifHr", "Highest Rating Inferred ROC"),
        uiOutput("FROCFormat")
      )
    }else if (dataset$dataType == "ROC"){
      values$cnvtedDataset <- values$dataset
      list(
        p(sprintf("Data type: %s", dataset$dataType)),        
        radioButtons("saveFormat", "Save As:",
                     c("JAFROC" = "JAFROC",
                       "MRMC(.csv)" = "MRMCCSV",
                       "MRMC(.lrc)" = "MRMCLRC",
                       "iMRMC" = "iMRMC"),
                     inline = TRUE
        ),
        downloadButton("saveCnvted", "Save")
      )
    }    
  })
  
  output$FROCFormat <- renderUI({
    if (input$ifHr){
      values$cnvtedDataset <- FROC2HrROC(values$dataset)
      list(
        radioButtons("saveFormat", "Save As:",
                     c("JAFROC" = "JAFROC",
                       "MRMC(.csv)" = "MRMCCSV",
                       "MRMC(.lrc)" = "MRMCLRC",
                       "iMRMC" = "iMRMC"),
                     inline = TRUE
        ),
        downloadButton("saveCnvted", "Save")
      )
    }else{
      values$cnvtedDataset <- values$dataset
      list(
        radioButtons("saveFormat", "Save As:",
                     c("JAFROC" = "JAFROC"),
                     inline = TRUE
        ),
        downloadButton("saveCnvted", "Save")
      )
    }
  })
  
  output$saveCnvted <- downloadHandler(
    filename = function() {
      if (input$saveFormat == "JAFROC"){
        paste("data", ".xlsx", sep="")
      }else if (input$saveFormat == "MRMCCSV"){
        paste("data", ".csv", sep="")
      }else if (input$saveFormat == "MRMCLRC"){
        paste("data", ".lrc", sep="")
      }else if (input$saveFormat == "iMRMC"){
        paste("data", ".imrmc", sep="")
      }      
    },
    
    content = function(file) {
      if (input$saveFormat == "JAFROC"){
        SaveDataFile(values$cnvtedDataset, fileName = file, format = "JAFROC")
      }else if (input$saveFormat == "MRMCCSV"){
        SaveDataFile(values$cnvtedDataset, fileName = file, format = "MRMC")
      }else if (input$saveFormat == "MRMCLRC"){
        SaveDataFile(values$cnvtedDataset, fileName = file, format = "MRMC")
      }else if (input$saveFormat == "iMRMC"){
        SaveDataFile(values$cnvtedDataset, fileName = file, format = "iMRMC")
      }       
    }
  )
  
  output$powerTable <- renderTable({
    if (input$calculateDataPowerBtn == 0){
      NULL
    }else{
      input$calculateDataPowerBtn
      isolate(PowerTable(values$dataset, alpha = as.numeric(input$alphaDataPower), 
                         effectSize = as.numeric(input$effectSizeDataPower),
                         desiredPower = as.numeric(input$desiredDataPower), 
                         randomOption = input$randomDataPower))
    }
  })
  
  output$sampleSize <- renderText({
    if (input$calculateSampleSizeBtn == 0){
      NULL
    }else{
      input$calculateSampleSizeBtn
      isolate({
        sampSize <- SampleSizeGivenJ(J = as.numeric(input$JSampleSize), 
                                     varYTR = as.numeric(input$varYTR), varYTC = as.numeric(input$varYTC), 
                                     varYEps = as.numeric(input$varYEps), cov1 = as.numeric(input$cov1),
                                     cov2 = as.numeric(input$cov2), cov3 = as.numeric(input$cov3), 
                                     varEps = as.numeric(input$varEps), msTR = as.numeric(input$msTR),
                                     KStar = as.numeric(input$KStar),
                                     alpha = as.numeric(input$alphaSampleSize), 
                                     effectSize = as.numeric(input$effectSizeSampleSize),
                                     desiredPower = as.numeric(input$desiredSampleSize), 
                                     randomOption = input$randomSampleSize)
        sprintf("The required number of cases for desired power %f is %d.", sampSize$power, sampSize$K)
      })
    }
  })
  
  output$varComp <- renderUI({
    if (input$compType == "DBMH"){
      list(
        textInput("varYTR", "VAR(TR)"),
        textInput("varYTC", "VAR(TC)"),
        textInput("varYEps", "VAR(ERROR)")
      )
    }else{
      list(
        textInput("cov1", "COV1"),
        textInput("cov2", "COV2"),
        textInput("cov3", "COV3"),
        textInput("varEps", "VAR(ERROR)"),
        textInput("msTR", "MS(T*R)"),
        textInput("KStar", "Number of Cases in Pilot Study")
      )
    }
  })
  
  output$plots <- renderUI({
    if (values$dataset$dataType == "ROC"){
      wellPanel(
        p("ROC Curve"),
        p(
          downloadButton("downloadROCPlot", "Save ROC Plot"),
          downloadButton("downloadROCPoints", "Save ROC Points")
        ),
        plotOutput("ROCPlot")
      )
    }else{
      list(
        splitLayout(
          wellPanel(
            p("ROC Curve"),
            p(
              downloadButton("downloadROCPlot", "Save ROC Plot"),
              downloadButton("downloadROCPoints", "Save ROC Points")
            ),
            plotOutput("ROCPlot")
          ),
          wellPanel(
            p("AFROC Curve"),
            p(
              downloadButton("downloadAFROCPlot", "Save AFROC Plot"),
              downloadButton("downloadAFROCPoints", "Save AFROC Points")
              ),
            plotOutput("AFROCPlot")
          )        
        ),
        splitLayout(
          wellPanel(
            p("FROC Curve"),
            p(
              downloadButton("downloadFROCPlot", "Save FROC Plot"),
              downloadButton("downloadFROCPoints", "Save FROC Points")
              ),
            plotOutput("FROCPlot")
          ),
          cellWidths = "50%"
        )
      )
    }
  })
  
  output$downloadROCPlot <- downloadHandler(
    filename = function() {
      paste("ROCPlot", ".png", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, values$ROCPlot, device = device)
    }
  )
  
  output$downloadROCPoints <- downloadHandler(
    filename = function() {
      paste("ROCPoints", ".csv", sep="")
    },
    content = function(file) {
      ROCPoints <- values$ROCPoints
      modalities <-  unlist(lapply(strsplit(as.character(ROCPoints$class), split = "\n"), "[[", 1))
      readers <- unlist(lapply(strsplit(as.character(ROCPoints$class), split = "\n"), "[[", 2))
      ROCPoints <- data.frame(FPF = ROCPoints$FPF, TPF = ROCPoints$TPF, Modality = modalities, Reader = readers)
      write.csv(ROCPoints, file, row.names = FALSE)
    }
  )
  
  output$downloadAFROCPlot <- downloadHandler(
    filename = function() {
      paste("AFROCPlot", ".png", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, values$AFROCPlot, device = device)
    }
  )
  
  output$downloadAFROCPoints <- downloadHandler(
    filename = function() {
      paste("AFROCPoints", ".csv", sep="")
    },
    content = function(file) {
      AFROCPoints <- values$AFROCPoints
      modalities <-  unlist(lapply(strsplit(as.character(AFROCPoints$class), split = "\n"), "[[", 1))
      readers <- unlist(lapply(strsplit(as.character(AFROCPoints$class), split = "\n"), "[[", 2))
      AFROCPoints <- data.frame(FPF = AFROCPoints$FPF, TPF = AFROCPoints$TPF, Modality = modalities, Reader = readers)
      write.csv(AFROCPoints, file, row.names = FALSE)
    }
  )
  
  output$downloadFROCPlot <- downloadHandler(
    filename = function() {
      paste("FROCPlot", ".png", sep="")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, values$FROCPlot, device = device)
    }
  )
  
  output$downloadFROCPoints <- downloadHandler(
    filename = function() {
      paste("FROCPoints", ".csv", sep="")
    },
    content = function(file) {
      FROCPoints <- values$FROCPoints
      modalities <-  unlist(lapply(strsplit(as.character(FROCPoints$class), split = "\n"), "[[", 1))
      readers <- unlist(lapply(strsplit(as.character(FROCPoints$class), split = "\n"), "[[", 2))
      FROCPoints <- data.frame(FPF = FROCPoints$FPF, TPF = FROCPoints$TPF, Modality = modalities, Reader = readers)
      write.csv(FROCPoints, file, row.names = FALSE)
    }
  )
  
  output$ROCPlot <- renderPlot({ 
    trts <- values$trts
    rdrs <- values$rdrs
    if (length(trts) == 0 || length(rdrs) == 0){
      NULL
    }else{
      if (length(trts) > 5 || length(rdrs) > 5){
        ROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, lgdPos = "right", opChType = "ROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$ROCPlot <- ROCPlot$ROCPlot
        values$ROCPoints <- ROCPlot$ROCPoints
        values$ROCPlot
      }else{
        ROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, opChType = "ROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$ROCPlot <- ROCPlot$ROCPlot
        values$ROCPoints <- ROCPlot$ROCPoints
        values$ROCPlot
      }
    }
  })
  
  output$AFROCPlot <- renderPlot({    
    trts <- values$trts
    rdrs <- values$rdrs
    if (length(trts) == 0 || length(rdrs) == 0){
      NULL
    }else{
      if (length(trts) > 5 || length(rdrs) > 5){
        AFROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, lgdPos = "right", opChType = "AFROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$AFROCPlot <- AFROCPlot$AFROCPlot
        values$AFROCPoints <- AFROCPlot$AFROCPoints
        values$AFROCPlot
      }else{
        AFROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, opChType = "AFROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$AFROCPlot <- AFROCPlot$AFROCPlot
        values$AFROCPoints <- AFROCPlot$AFROCPoints
        values$AFROCPlot
      }
    }
  })
  
  output$FROCPlot <- renderPlot({   
    trts <- values$trts
    rdrs <- values$rdrs
    if (length(trts) == 0 || length(rdrs) == 0){
      NULL
    }else{
      if (length(trts) > 5 || length(rdrs) > 5){
        FROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, lgdPos = "right", opChType = "FROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$FROCPlot <- FROCPlot$FROCPlot
        values$FROCPoints <- FROCPlot$FROCPoints
        values$FROCPlot
      }else{
        FROCPlot <- EmpiricalOpCharac(values$dataset, trts, rdrs, opChType = "FROC")
        values$trts <- trts
        values$rdrs <- rdrs
        values$FROCPlot <- FROCPlot$FROCPlot
        values$FROCPoints <- FROCPlot$FROCPoints
        values$FROCPlot
      }
    }
  })
  
  output$downloadReportBtn <- downloadHandler(
    filename = function() {
      paste0(paste(file_path_sans_ext(input$dataFile$name), input$mthd, input$fom, sep="_"), ".txt")
    },
    content = function(file) {
      write(values$reportStrings, file)
    }
  )
  
  output$report <- renderText({
    if (input$analyzeBtn == 0){
      "Click \"Analyze\" to generate analysis report."
    }else{      
      input$analyzeBtn
      values$reportStrings <- isolate(ReportTextForGUI(dataset = values$dataset, method = input$mthd, fom = input$fom, 
                                                       alpha = as.numeric(input$alpha), covEstMethod = input$covEstMthd, 
                                                       nBoots = as.numeric(input$nBoots)))
    }
  })
  
  output$plotboxes <- renderUI({
    trtGroup <- values$trtGroup
    rdrGroup <- values$rdrGroup
    ret <- reactiveValuesToList(input)
    plotList <- list()
    for (i in 1:(input$nGroup + 1)){
      if (paste0("plotTrts", i) %in% names(ret)){
        plotList <- c(plotList,                     
                      list(
                        wellPanel(
                          checkboxGroupInput(paste0("plotTrts", i), 
                                             "Modalities: ",
                                             trtGroup,
                                             get(paste0("plotTrts", i), ret),
                                             inline = TRUE),
                          checkboxGroupInput(paste0("plotRdrs", i),
                                             "Readers: ",
                                             rdrGroup,
                                             get(paste0("plotRdrs", i), ret),
                                             inline = TRUE),  
                          checkboxInput(paste0("ifAvg", i),
                                        "Averaged Curve?",
                                        get(paste0("ifAvg", i), ret))
                        )
                      )
        )
      }else{
        plotList <- c(plotList,                     
                      list(
                        wellPanel(
                          checkboxGroupInput(paste0("plotTrts", i), 
                                             "Modalities: ",
                                             trtGroup,
                                             inline = TRUE),
                          checkboxGroupInput(paste0("plotRdrs", i),
                                             "Readers: ",
                                             rdrGroup,
                                             inline = TRUE),  
                          checkboxInput(paste0("ifAvg", i),
                                        "Averaged Curve?")
                        )
                      )
        )
      }
    }    
    selectGroup <- NULL
    for (i in 1:(1 + input$nGroup)){
      groupTrt <- as.numeric(ret[[paste0("plotTrts", i)]])
      groupRdr <- as.numeric(ret[[paste0("plotRdrs", i)]])
      if (length(groupTrt) == 0 || length(groupRdr) == 0) next
      selectGroup <- c(selectGroup, i)
    }
    
    if (!is.null(selectGroup) && length(selectGroup) == 1 && !ret[[paste0("ifAvg", selectGroup)]]){
      trts <- as.numeric(ret[[paste0("plotTrts", selectGroup)]])
      rdrs <- as.numeric(ret[[paste0("plotRdrs", selectGroup)]])
    }else{
      trts <- list()
      rdrs <- list()
      for (i in 1:(input$nGroup + 1)){
        ifAvg <- ret[[paste0("ifAvg", i)]]
        if (length(ifAvg) != 0 && ifAvg){
          groupTrt <- as.numeric(ret[[paste0("plotTrts", i)]])
          groupRdr <- as.numeric(ret[[paste0("plotRdrs", i)]])
          if (length(groupTrt) == 0 || length(groupRdr) == 0) next
          trts <- c(
            trts,
            list(groupTrt)
          )
          rdrs <- c(
            rdrs,
            list(groupRdr)
          )
        }else{
          allComb <- expand.grid(as.list(as.numeric(ret[[paste0("plotTrts", i)]])), 
                                 as.list(as.numeric(ret[[paste0("plotRdrs", i)]])))
          if (nrow(allComb) == 0) next
          trts <- c(
            trts,
            as.list(as.numeric(allComb[ , 1]))
          )
          rdrs <- c(
            rdrs,
            as.list(as.numeric(allComb[ , 2]))
          )
        }
      }
    }
    values$trts <- trts
    values$rdrs <- rdrs
    plotList
  })
  
  output$covEstBtn <- renderUI({
    if (is.null(input$mthd)){
      NULL
    }else if (input$mthd == "DBMH"){
      NULL
    }else if (input$mthd == "ORH") {
      if (input$fom %in% c("Wilcoxon", "HrAuc", "ROI")){
        if (is.null(input$covEstMthd)){
          radioButtons("covEstMthd", "Covariances Estimate Method:",
                       c("Jackknife" = "Jackknife",
                         "Bootstrap" = "Bootstrap",
                         "DeLong" = "DeLong"
                       ),
                       inline = TRUE
          )
        }else{
          radioButtons("covEstMthd", "Covariances Estimate Method:",
                       c("Jackknife" = "Jackknife",
                         "Bootstrap" = "Bootstrap",
                         "DeLong" = "DeLong"
                       ),
                       inline = TRUE,
                       selected = input$covEstMthd
          )
        }
      }else{
        if (is.null(input$covEstMthd) || input$covEstMthd == "DeLong"){
          radioButtons("covEstMthd", "Covariances Estimate Method:",
                       c("Jackknife" = "Jackknife",
                         "Bootstrap" = "Bootstrap"
                       ),
                       inline = TRUE
          )
        }else{
          radioButtons("covEstMthd", "Covariances Estimate Method:",
                       c("Jackknife" = "Jackknife",
                         "Bootstrap" = "Bootstrap"
                       ),
                       inline = TRUE,
                       selected = input$covEstMthd
          )          
        }
      }
    }
  })
  
  output$nBootsBtn <- renderUI({
    if (input$mthd == "DBMH"){
      NULL
    }else{
      if (is.null(input$covEstMthd)){
        NULL
      }else if (input$covEstMthd == "Bootstrap"){
        textInput("nBoots", "Number of Bootstrapping", value = 200)
      }else{
        NULL
      }
    }
  })
  
  output$truthTable <- renderTable({ 
    dataset <- values$dataset
    NL <- dataset$NL
    LL <- dataset$LL
    lesionNum <- dataset$lesionNum
    lesionID <- dataset$lesionID
    lesionWeight <- dataset$lesionWeight
    maxNL <- dim(NL)[4]
    dataType <- dataset$dataType
    modalityID <- dataset$modalityID
    readerID <- dataset$readerID
    I <- length(modalityID)
    J <- length(readerID)
    K <- dim(NL)[3]
    K2 <- dim(LL)[3]
    K1 <- K - K2
    
    caseIDs <- c(1:K1, rep(K1 + 1:K2, lesionNum))
    lesionIDs <- as.vector(t(lesionID))
    lesionIDs <- lesionIDs[lesionIDs != -Inf]
    lesionIDs <- c(rep(0, K1), lesionIDs)
    lesionWeights <- as.vector(t(lesionWeight))
    lesionWeights <- lesionWeights[lesionWeights != -Inf]
    lesionWeights <- c(rep(0, K1), lesionWeights)
    data.frame(CaseID = caseIDs, LesionID = as.integer(lesionIDs), Weight = lesionWeights)
  })
  
  output$nlTable <- renderTable({  
    dataset <- values$dataset
    NL <- dataset$NL
    LL <- dataset$LL
    lesionNum <- dataset$lesionNum
    lesionID <- dataset$lesionID
    lesionWeight <- dataset$lesionWeight
    maxNL <- dim(NL)[4]
    dataType <- dataset$dataType
    modalityID <- dataset$modalityID
    readerID <- dataset$readerID
    I <- length(modalityID)
    J <- length(readerID)
    K <- dim(NL)[3]
    K2 <- dim(LL)[3]
    K1 <- K - K2
    dataSheet <- NULL
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          for (l in 1:maxNL) {
            if (NL[i, j, k, l] != -Inf) {
              dataSheet <- rbind(dataSheet, c(j, i, k, NL[i, j, k, l]))
            }
          }
        }
      }
    }
    data.frame(ReaderID = readerID[dataSheet[, 1]], ModalityID = modalityID[dataSheet[, 2]], CaseID = as.integer(dataSheet[, 3]), NL_Rating = signif(dataSheet[, 4], 6))
  })
  
  output$llTable <- renderTable({   
    dataset <- values$dataset
    NL <- dataset$NL
    LL <- dataset$LL
    lesionNum <- dataset$lesionNum
    lesionID <- dataset$lesionID
    lesionWeight <- dataset$lesionWeight
    maxNL <- dim(NL)[4]
    dataType <- dataset$dataType
    modalityID <- dataset$modalityID
    readerID <- dataset$readerID
    I <- length(modalityID)
    J <- length(readerID)
    K <- dim(NL)[3]
    K2 <- dim(LL)[3]
    K1 <- K - K2
    dataSheet <- NULL
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K2) {
          for (l in 1:lesionNum[k]) {
            if (LL[i, j, k, l] != -Inf) {
              dataSheet <- rbind(dataSheet, c(j, i, k + K1, lesionID[k, l], LL[i, j, k, l]))
            }
          }
        }
      }
    }
    data.frame(ReaderID = readerID[dataSheet[, 1]], ModalityID = modalityID[dataSheet[, 2]], CaseID = as.integer(dataSheet[, 3]), LesionID = as.integer(dataSheet[, 4]), LL_Rating = signif(dataSheet[, 5], 6))
  })
})