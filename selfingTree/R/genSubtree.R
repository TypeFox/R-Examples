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

genSubtree.2M <- function(gam1,gam2){

    
    F <- list()

    ## generate all possible gametes from parent genotype
    L1 <- paste(substr(gam1,1,1),substr(gam2,1,1),sep="")
    L2 <- paste(substr(gam1,2,2),substr(gam2,2,2),sep="")

    AB.series <- unique(paste(rep(c(substr(L1,1,1),substr(L1,2,2)),each = 2),
                              rep(c(substr(L2,1,1),substr(L2,2,2)),2),
                              sep = ""))

    for(AB1 in AB.series){
        for(AB2 in AB.series){

            geno.prob <- paste(haploProb.2M(gam1,gam2,AB1),
                               haploProb.2M(gam1,gam2,AB2),
                               sep = "*")

            F[[length(F) + 1]] <- list(rbind(AB1,AB2),geno.prob)
       
        }
    }

    return(F)    
}



genSubtree.3M <- function(gam1,gam2){
    
    F <- list()

    L1 <- paste(substr(gam1,1,1),substr(gam2,1,1),sep="")
    L2 <- paste(substr(gam1,2,2),substr(gam2,2,2),sep="")
    L3 <- paste(substr(gam1,3,3),substr(gam2,3,3),sep="")


    AB.series <- unique(paste(rep(c(substr(L1,1,1),substr(L1,2,2)),each = 4),
                              rep(rep(c(substr(L2,1,1),substr(L2,2,2)),each = 2),2),
                              rep(c(substr(L3,1,1),substr(L3,2,2)), 4),
                              sep = ""))

    for(AB1 in AB.series){
        for(AB2 in AB.series){

            geno.prob <- paste(haploProb.3M(gam1,gam2,AB1),
                               haploProb.3M(gam1,gam2,AB2),
                               sep = "*")

            F[[length(F) + 1]] <- list(rbind(AB1,AB2),geno.prob)
       
        }
    }

    return(F)
    
}
