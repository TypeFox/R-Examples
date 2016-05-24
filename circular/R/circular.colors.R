#############################################################
#                                                           
#   circular.colors function                                       
#   Author: Claudio Agostinelli                             
#   E-mail: claudio@unive.it                                
#   Date: October, 7, 2007                                  
#   Version: 0.7                                            
#                                                           
#   Copyright (C) 2007 Claudio Agostinelli                  
#                                                           
#############################################################

circular.colors <- function(n, m=0, M=2*pi, offset=0, ...) {
  hh <- seq(from=(m-offset)%%(2*pi), to=if((M-offset)==2*pi) (M-offset) else (M-offset)%%(2*pi), length.out=n)/(M-m)
  hsv(h=hh, ...)
}
