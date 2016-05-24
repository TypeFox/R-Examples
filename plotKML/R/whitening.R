# Purpose        : Derive whitenned color based on the uncertainty;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : this technique requires a special 2D legend;

whitening <- function(
   z,      # target variable
   zvar,   # associated uncertainty
   zlim = c(min(z, na.rm=TRUE), max(z, na.rm=TRUE)), 
   elim = c(.4,1), 
   global.var = var(z, na.rm=TRUE), 
   col.type = "RGB")  # output col.type can be "RGB" or "hex" 
   {
   
   # Derive the normalized error:
   er <- sqrt(zvar)/sqrt(global.var)
   # Strech the values (z) to the inspection range:
   tz <- (z-zlim[1])/(zlim[2]-zlim[1])
   tz <- ifelse(tz<=0, 0, ifelse(tz>1, 1, tz))
   # Derive the Hues:
   f1 <- -90-tz*300
   f2 <- ifelse(f1<=-360, f1+360, f1)
   H <- ifelse(f2>=0, f2, (f2+360))
   # Strech the error values (e) to the inspection range:
   er <- (er-elim[1])/(elim[2]-elim[1])
   er <- ifelse(er<=0, 0, ifelse(er>1, 1, er))
   # Derive the saturation and intensity images:
   S <- 1-er
   V <- 0.5*(1+er)
   
   # Convert the HSV values to RGB and put them as R, G, B bands:
   if(col.type=="hex"){
      out.cols <- hex(HSV(H, S, V))
      return(out.cols)
   }
   else { 
      out.cols <- as(HSV(H, S, V), "RGB")
      return(out.cols)
} 
} 

# end of script;
