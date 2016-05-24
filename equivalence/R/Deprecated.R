

## <entry>
## Deprecated in 2.5.0
tost.data <- function(x, null = 0, alpha = 0.05, Epsilon = 0.36, absolute = FALSE) {
    .Deprecated("tost")
    if (!absolute) stop("Relative intervals should use ptte, *not* tost.")
    tost(x, y=NULL, mu=null, alpha=alpha, epsilon=Epsilon)
}
## </entry>


