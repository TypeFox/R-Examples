#Function called within Gst.est and Dest and H.out

harmonic<-function(values,n){
                                     out <- n/sum(1/values)
                                     invisible(out)
                                      }
