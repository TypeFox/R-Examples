## Copyright (c) 2014, Pioneer Hi-Bred International, Inc.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.

##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.

##     Neither the name of Pioneer Hi-Bred International, Inc. nor the
##     names of its contributors may be used to endorse or promote
##     products derived from this software without specific prior written
##     permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


globalVariables(c('i','j','k','l','m'))
nodeProbabilities <- function(F,generation){
    
    prob.list <- list()
    
    prob.list[["F2"]] <- extractProbs(F)
    
    ## F3
    if(generation > 2){
        F2 <- extractProbs(F)
        prob.list[["F3"]] <- foreach(i = 1:length(F2),.combine = "c") %do% {
        
            F3 <- extractProbs(F[[i]][[3]])

            prob.temp <- paste(F2[[i]],"*",F3)
            names(prob.temp) <- names(F3)
            return(prob.temp)
        
        }
                
    }
    
    ## F4
    if(generation > 3){
        F2 <- extractProbs(F)
        prob.list[["F4"]] <- foreach(i = 1:length(F2),.combine = "c") %do% {
        
            F3 <- extractProbs(F[[i]][[3]])
        
            F4.temp <- foreach(j = 1:length(F3),.combine = "c") %do% {
            
                F4 <- extractProbs(F[[i]][[3]][[j]][[3]])
            
                prob.temp <- paste(F2[[i]],"*",
                                   F3[[j]],"*",
                                   F4)
            
                names(prob.temp) <- names(F4)
                return(prob.temp)
            }

            return(F4.temp)

        }
    }

    
    ## F5
    if(generation > 4){
        F2 <- extractProbs(F)
        prob.list[["F5"]] <- foreach(i = 1:length(F2),.combine = "c") %do% {
    
            F3 <- extractProbs(F[[i]][[3]])
    
            F4.temp <- foreach(j = 1:length(F3),.combine = "c") %do% {    

                F4 <- extractProbs(F[[i]][[3]][[j]][[3]])
            
                F5.temp <- foreach(k = 1:length(F4),.combine = "c") %do% {

                    F5 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]])
            
                    prob.temp <- paste(F2[[i]],"*",
                                       F3[[j]],"*",
                                       F4[[k]],"*",
                                       F5)
        
                    names(prob.temp) <- names(F5)
                    return(prob.temp)
    
                }

                return(F5.temp)
            }
    
            return(F4.temp)
        }
    }

    ## F6
    if(generation > 5){
        F2 <- extractProbs(F)
        prob.list[["F6"]] <- foreach(i = 1:length(F2),.combine = "c") %do% {
    
            F3 <- extractProbs(F[[i]][[3]])
    
            F4.temp <- foreach(j = 1:length(F3),.combine = "c") %do% {    
        
                F4 <- extractProbs(F[[i]][[3]][[j]][[3]])
        
                F5.temp <- foreach(k = 1:length(F4),.combine = "c") %do% {

                    F5 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]])
            
                    F6.temp <- foreach(l = 1:length(F5),.combine = "c") %do% {
                
                        F6 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]])
                
                        prob.temp <- paste(F2[[i]],"*",
                                           F3[[j]],"*",
                                           F4[[k]],"*",
                                           F5[[l]],"*",
                                           F6)
                
                        names(prob.temp) <- names(F6)
                
                        return(prob.temp)
                    }
                    return(F6.temp)
                }
                return(F5.temp)
            }
    
            return(F4.temp)
        }
    }


    ## F7
    if(generation > 6){
        F2 <- extractProbs(F)
        prob.list[["F7"]] <- foreach(i = 1:length(F2),.combine = "c") %do% {
    
            F3 <- extractProbs(F[[i]][[3]])
    
            F4.temp <- foreach(j = 1:length(F3),.combine = "c") %do% {    
        
                F4 <- extractProbs(F[[i]][[3]][[j]][[3]])
        
                F5.temp <- foreach(k = 1:length(F4),.combine = "c") %do% {

                    F5 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]])
            
                    F6.temp <- foreach(l = 1:length(F5),.combine = "c") %do% {
                
                        F6 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]])

                        F7.temp <- foreach(m = 1:length(F6),.combine = "c") %do% {
                
                            F7 <- extractProbs(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]][[m]][[3]])

                
                            prob.temp <- paste(F2[[i]],"*",
                                               F3[[j]],"*",
                                               F4[[k]],"*",
                                               F5[[l]],"*",
                                               F6[[m]],"*",
                                               F7)
                
                            names(prob.temp) <- names(F7)
                
                            return(prob.temp)
                        }
                        return(F7.temp)
                    }
                    return(F6.temp)
                }
                return(F5.temp)
            }    
            return(F4.temp)
        }
    }
        
    return(prob.list)
    
}
