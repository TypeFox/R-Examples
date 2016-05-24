fp.scale.coefs <- function
(
  x
  )
### pseudo function
### taken from package <mfp>
### 2008-06
{
        x.ok <- na.omit(x)
        
        scale <- 1
        shift <- 0
        if (min(x.ok) <= 0) {
                z <- diff(sort(x.ok))
                shift <- min(z[z > 0]) - min(x.ok)
                shift <- ceiling(shift * 10)/10
        }
        range <- mean(x.ok + shift)
        scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
        
        c(shift=shift, scale=scale, rep(-999, length(x) - 2))
}
