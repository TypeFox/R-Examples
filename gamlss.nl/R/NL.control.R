"NL.control" <- function(fscale = 1,   
                      typsize = NULL, #abs(p0), 
                      stepmax = NULL, # sqrt(p0 %*% p0),  
                      iterlim = 100,  
                       ndigit = 10,   
                      steptol = 1e-05, 
                      gradtol = 1e-05,
                  print.level = 0, 
            check.analyticals = TRUE,
                      hessian = TRUE )
{
##  Control iteration for GLIM
##  MS  Sunday, February 17, 2002 at 19:18
##
        if(fscale <= 0) 
        {
warning("the scale value supplied is zero or negative the default value of 1 was used instead")
          fscale <- 1
        }
 #       if(any(typsize <= 0)) 
 #       {
#warning("the value of typsize supplied is zero or negative the default value of abs(p0) was used instead")
#         typsize <- abs(p0)
#        }
#        if(stepmax < 0) 
#        {
#warning("the value of stepmax supplied is zero or negative the default value of  sqrt(p0 %*% p0) was used instead")
#         stepmax <-  sqrt(p0 %*% p0)
#        }
        if( iterlim <= 0) 
        {
warning("the value of  iterlim supplied is zero or negative the default value of 100 was used instead")
          iterlim <- 100
        } 
         if( ndigit < 0) 
        {
warning("the value of  ndigit supplied is zero or negative the default value of 10 was used instead")
           ndigit <- 10
        }   
         if(  steptol < 0) 
        {
warning("the value of  steptol supplied is zero or negative the default value of  1e-05 was used instead")
          steptol <- 1e-05
       }
         if(  gradtol < 0) 
       { 
warning("the value of  gradtol supplied is zero or negative the default value of  1e-05 was used instead")
          gradtol <- 1e-05
       } 
         if( print.level < 0) 
       {
 warning("the value of  print.level supplied is negative the default value of  0 was used instead")
      print.level <- 0
       }           
       list(fscale = fscale, typsize = typsize, stepmax = stepmax,  
            iterlim = iterlim, 
            ndigit = ndigit, steptol = steptol, gradtol = gradtol,  
            print.level =  print.level,
            check.analyticals = as.logical(check.analyticals)[1],
            hessian = as.logical(hessian)[1])
}
