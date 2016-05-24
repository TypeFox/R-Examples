fp.scale <- function
(
 x
 )
### taken from package <mfp>
### 2008-06
{
        scale <- 1
        shift <- 0
        if (min(x) <= 0) {
                z <- diff(sort(x))
                shift <- min(z[z > 0]) - min(x)
                shift <- ceiling(shift * 10)/10
        }
        range <- mean(x + shift)
        scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
        
        list(
             shift=shift,
             scale=scale
             )
}
