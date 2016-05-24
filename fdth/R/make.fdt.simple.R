make.fdt.simple <- function (x,
                             start,
                             end,
                             h,
                             right)

{
  ## Ivan B. Allaman: thank you to the bug fix
  f <- table(cut(x,
                 br=seq(start,
                        end,
                        h),
                 right=right,
                 dig.lab=nchar(as.character(round(max(x),
                                                  2)))))    # Absolute freq.
  rf  <- as.numeric(f/length(x))                            # Relative freq
  rfp <- as.numeric(100*(f/length(x)))                      # Relative freq, %
  cf  <- as.numeric(cumsum(f))                              # Cumulative freq
  cfp <- as.numeric(100*(cumsum(f/length(x))))              # Cumulative freq, %

  res <- data.frame(f,                                      # Make final table
                    rf,
                    rfp,
                    cf,
                    cfp)                   

  names(res) <- c('Class limits',
                  'f',
                  'rf',
                  'rf(%)',
                  'cf', 
                  'cf(%)')

  return(res)
}
