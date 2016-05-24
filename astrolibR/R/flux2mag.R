flux2mag = function( flux, zero_pt=21.10, ABwave) {

if(missing(ABwave) )
     return( -2.5*log10(flux) - zero_pt)
else 
     return(-2.5*log10(flux) - 5*log10(ABwave) - 2.406)
}
