
#' @export

p0G_Func <- function(p, status, relative, model= "dominant")
{
   # status = 0, non-carrier; = 1, carrier; = 2, Homozygous; = 3, Heterozygous.
   # relative = 1, parents; = 2, sibling; = 3, offspring.
   # model = "dominant" or "recessive".

   prob = 0;

   if(model == "dominant")
   {
     if(status == 0)
     {
       prob = (relative == 1)*p + 
              (relative == 2)*(p-0.25*p^2) +                # 0.25*p^2 + 0.5*(-p^2+2*p)
              (relative == 3)*p;
     }
     
     if(status == 1)
     {
       prob = (relative == 1)*( -0.5*p^2 + p + 0.5 ) +      # p*(p + 1 - p) + (1-p)*(0.5*p + 0.5)
              (relative == 2)*( -0.5*p^2 + 0.5 ) +          # p*(0.25*(p^2-2*p+1) + 0.5*(1-p^2)) + (1-p)*(0.25*(p^2+p) + 0.5*(-p^2 - p + 1))
              (relative == 3)*( -0.5*p^2 + p + 0.5 );       # p*(p + 1 - p) + (1-p)*(0.5*p + 0.5)
     }
     
     if(status == 2)
     {
       prob = (relative == 1)*1 +                           # ( p + 1 - p )
              (relative == 2)*(-0.25*p^2 -0.5*p +0.75) +    # 0.25*(p^2-2*p+1) + 0.5*(1-p^2)
              (relative == 3)*1;
     }

     if(status == 3)
     {
       prob = (relative == 1)*( 0.5*p + 0.5 ) + 
              (relative == 2)*( -0.25*p^2 -0.25*p +0.5 ) +   # 0.25*(p^2+p) + 0.5*(-p^2 - p + 1)
              (relative == 3)*( 0.5*p + 0.5 );
     }
   }
 
   if(model == "recessive")
   {
     if(status == 0)
     {
       prob = (relative == 1)*( 1-p ) + 
              (relative == 2)*( 0.25*p^2 - p + 1 ) + 
              (relative == 3)*( 1-p );
     }

     if(status == 1)
     {
       prob = (relative == 1)*( 0.5*(1-p)^2 ) + 
              (relative == 2)*( 0.5*(1-p)^2 ) +             # p*(0.25*(1-p)^2) + (1-p)*0.25*(p^2-3*p+2) 
              (relative == 3)*( 0.5*(1-p)^2 );
     }
     
     if(status == 2)
     {
       prob = (relative == 1)*0 + 
              (relative == 2)*( 0.25*(1-p)^2 ) + 
              (relative == 3)*0;     
     }

     if(status == 3)
     {
       prob = (relative == 1)*( 0.5*(1-p) ) + 
              (relative == 2)*( 0.25*(p^2-3*p+2) ) + 
              (relative == 3)*( 0.5*(1-p) );      
     }
   }

  return(prob)
}
