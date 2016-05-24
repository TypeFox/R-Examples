# This substitutes for `poly' in IDL scripts
# See http://www.exelisvis.com/docs/POLY.html

polyidl = function(x,cc) { 
   polysum=0 
   for (i in 1:length(cc)) polysum = polysum + cc[i]*x^{i-1}
   return(polysum)
   }   
