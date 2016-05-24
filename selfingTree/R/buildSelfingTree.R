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


buildSelfingTree <- function(genF,generation,gam1,gam2){
    
    ## F2
    F <- genF(gam1,gam2)

    ## F3
    if(generation > 2){
        for(i in 1:length(F)){
            F[[i]][[3]] <- genF(F[[i]][[1]][1,],F[[i]][[1]][2,])
        }
    }

    ## F4
    if(generation > 3){
        for(i in 1:length(F))
            for(j in 1:length(F[[i]][[3]]))
                F[[i]][[3]][[j]][[3]] <- genF(F[[i]][[3]][[j]][[1]][1,],F[[i]][[3]][[j]][[1]][2,])
    }
    
    ## F5
    if(generation > 4){
        for(i in 1:length(F))
            for(j in 1:length(F[[i]][[3]]))
                for(k in 1:length(F[[i]][[3]][[j]][[3]]))
                    F[[i]][[3]][[j]][[3]][[k]][[3]] <- genF(F[[i]][[3]][[j]][[3]][[k]][[1]][1,],
                                                            F[[i]][[3]][[j]][[3]][[k]][[1]][2,])
    }

    ## F6
    if(generation > 5){
        for(i in 1:length(F))
            for(j in 1:length(F[[i]][[3]]))
                for(k in 1:length(F[[i]][[3]][[j]][[3]]))
                    for(l in 1:length(F[[i]][[3]][[j]][[3]][[k]][[3]]))
                        F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]] <- genF(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[1]][1,],
                                                                          F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[1]][2,])
    }

    ## F7
    if(generation > 6){
        for(i in 1:length(F))
            for(j in 1:length(F[[i]][[3]]))
                for(k in 1:length(F[[i]][[3]][[j]][[3]]))
                    for(l in 1:length(F[[i]][[3]][[j]][[3]][[k]][[3]]))
                        for(m in 1:length(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]]))
                            F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]][[m]][[3]] <-
                                genF(F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]][[m]][[1]][1,],
                                     F[[i]][[3]][[j]][[3]][[k]][[3]][[l]][[3]][[m]][[1]][2,])
    }

    return(F)
}
