fences <-
function(xx, units = "ppm", display = TRUE)
{
     # Function to display the mean, SD, median and MAD, and then the mean +- 2SD,
     # median +- 2MAD, Tukey whisker end fences and the 98th %ile; calculations
     # are made without (top lines) and with (middle lines) a log10 transformation,
     # and with a logit transformation (lowest lines).
     #
     # Note: the logit transformation requires an argument in the range 0 to 1,
     # therefore the input values have to be divided by 100, 10^6, 10^9 or 10^12
     # if the reporting units are percent ("pct"), "ppm" (micro_g/g, mg/kg),
     # "ppb" (pico_g/g, micro_g/kg), or "ppt" (pico_g/kg), respectively.  This is
     # done internally, with the user providing the units, the default is "ppm".
     # 
     # Function may be used with fences.summary to generate fences based on
     # some spatial or contextual framework, i.e. output grouped by a classification
     # factor, that are saved as a file.  The parameter display, when set to FALSE,
     # is used to suppress output to the current device when the function is used
     # with fences.summary.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required, to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     qntls <- numeric(7)
     table <- numeric(40)
     temp.x <- remove.na(xx, iftell = FALSE)
     x <- temp.x$x[1:temp.x$n]
     nx <- temp.x$n
     nna <- temp.x$nna
     if(nx > 2) {
         qntls <- quantile(x, prob = c(0, 0.02, 0.25, 0.5, 0.75, 0.98, 1))
         table[1] <- nx
         table[2] <- nna
         table[3] <- mean(x)
         table[4] <- sqrt(var(x))
         table[5] <- qntls[4]
         table[6] <- mad(x)
         table[7] <- table[3] + 2 * table[4]
         table[8] <- table[5] + 2 * table[6]
         iqr <- qntls[5] - qntls[3]
         table[9] <- qntls[5] + 1.5 * iqr
         table[25] <- max(x[x < table[9]])
         table[10] <- qntls[6]
         table[11] <- table[3] - 2 * table[4]
         table[12] <- table[5] - 2 * table[6]
         table[13] <- qntls[3] - 1.5 * iqr
         table[26] <- min(x[x > table[13]])
         table[14] <- qntls[2]
         #
         lx <- log10(x)
         table[15] <- mean(lx)
         table[16] <- sqrt(var(lx))
         table[17] <- log10(qntls[4])
         table[18] <- mad(lx)
         table[19] <- 10^(table[15] + 2 * table[16])
         table[20] <- 10^(table[17] + 2 * table[18])
         liqr <- log10(qntls[5]) - log10(qntls[3])
         table[21] <- 10^(log10(qntls[5]) + 1.5 * liqr)
         table[27] <- max(x[x < table[21]])
         table[22] <- 10^(table[15] - 2 * table[16])
         table[23] <- 10^(table[17] - 2 * table[18])
         table[24] <- 10^(log10(qntls[3]) - 1.5 * liqr)
         table[28] <- min(x[x > table[24]])
         #
         k <- 10^6
         if(units == "pct") k <- 100
         if(units == "ppb") k <- 10^9
         if(units == "ppt") k <- 10^12
         lx <- logit(x/k)
         table[29] <- mean(lx)
         table[30] <- sqrt(var(lx))
         table[31] <- logit(qntls[4] / k)
         table[32] <- mad(lx)
         table[33] <- expit(table[29] + 2 * table[30]) * k
         table[34] <- expit(table[31] + 2 * table[32]) * k
         liqr <- logit(qntls[5] / k) - logit(qntls[3] / k)
         table[35] <- expit(logit(qntls[5] / k) + 1.5 * liqr) * k
         table[39] <- max(x[x < table[35]])
         table[36] <- expit(table[29] - 2 * table[30]) * k
         table[37] <- expit(table[31] - 2 * table[30]) * k
         table[38] <- expit(logit(qntls[3] / k) - 1.5 * liqr) * k
         table[40] <- min(x[x > table[38]])
     }
     else {
         table[1] <- nx
         table[2] <- nna
         table[3:28] <- NA
         if(nx == 1) {
             table[5] <- table[3] <- temp.x$x[1:temp.x$n]
             table[17] <- table[15] <- log10(table[3])
             table[31] <- table[29] <- logit(table[3])
         }
     }
     table[3:40] <- signif(table[3:40], 3)
     if(display) {
         cat(" ", deparse(substitute(xx)), "(Units =", units, ") :  N =", table[1], "    NAs =",
          table[2], "\t2%ile =", table[14], "\t98%ile =", table[10], 
          "\n\t Mean\t SD\t  Median   MAD\t\t Mean\t\t Med\t       Tukey Fences\n",
          "\t\t\t\t\t\t \2612SD\t\t \2612MAD\t\t     (actual)\n\t", 
          table[3], "\t", table[4], "\t  ", table[5], "  ", table[6], "\t", table[7], "\t\t", table[8], 
          "\t      ", table[9], " (", table[25], ")",
          "\n\t\t\t\t\t\t", table[11], "\t\t", table[12], "\t      ", table[13], " (", table[26], ")",
          "\n  Log10\t", table[15], "\t", table[16], "\t  ", table[17], "  ", table[18], "\t", table[19],
          "\t\t", table[20],   "\t      ", table[21], " (", table[27], ")", 
          "\n\t\t\t\t\t\t", table[22], "\t\t", table[23], "\t      ", table[24], " (", table[28], ")",
          "\n  Logit\t", table[29], "\t", table[30], "\t  ", table[31], " ", table[32], "\t", table[33],
          "\t\t", table[34],   "\t      ", table[35], " (", table[39], ")", 
          "\n\t\t\t\t\t\t", table[36], "\t\t", table[37], "\t      ", table[38], " (", table[40], ")\n")
     }
     invisible(table)
}
