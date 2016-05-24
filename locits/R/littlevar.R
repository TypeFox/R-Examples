littlevar <-
function (WP, ll) 
{
    WPll <- accessD(WP, lev = ll)
    2 * mean(WPll^2)
}
