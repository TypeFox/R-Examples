as.Date.cal.yr <-
function( x, ... )
{
structure( round( ( x - 1970 ) * 365.25 ), class="Date" )
}
