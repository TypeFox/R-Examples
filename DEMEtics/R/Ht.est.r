# Function called within Dest.r and Gst.est.r

          Ht.est<-function(Ht,Hs.est,N,n){
                                            out <- Ht+ (Hs.est/(2*N*n))
                                            invisible(out)
                                            }
