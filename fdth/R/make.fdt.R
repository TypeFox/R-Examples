## Make an full fdt from a minimal information: a vector of frequency and class limits
## Useful to reproduce any previous continuous fdt
make.fdt <- function(f,
                     start,
                     end,
                     right=FALSE)

{
  stopifnot(start < end)

  f   <- as.numeric(f)             # Frequency
  rf  <- f/sum(f)                  # Relative freq
  rfp <- 100*(f/sum(f))            # Relative freq, %
  cf  <- cumsum(f)                 # Cumulative freq
  cfp <- 100*(cumsum(f/sum(f)))    # Cumulative freq, %

  brk <- seq(start, 
             end, 
             len=length(f) + 1)

  if (right)
    cl <- paste('(',
                brk[1:length(brk) - 1],
                ',',
                brk[2:length(brk)],
                ']',
                sep='')
  else
    cl <- paste('[',
                brk[1:length(brk) - 1],
                ',',
                brk[2:length(brk)],
                ')',
                sep='')

  fdt <- data.frame(cl,            # Make final table
                    f,                                   
                    rf,
                    rfp,
                    cf,
                    cfp)                   

  names(fdt) <- c('Class limits',
                  'f',
                  'rf',
                  'rf(%)',
                  'cf', 
                  'cf(%)')

  breaks <- c(start,
              end,
              h=diff(brk)[1],
              ifelse (right,
                      1,
                      0))

  names(breaks) <- c('start',
                     'end',
                     'h',
                     'right')

  res <- list(table=fdt,
              breaks=breaks)

  class(res) <- c('fdt.default',
                  'fdt',
                  'list')

  invisible(res)
}
