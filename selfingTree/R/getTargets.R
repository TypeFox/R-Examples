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

globalVariables(c('g1','g2'))
getTargets <- function(target.geno){

    M <- nchar(target.geno)
    
    if(M == 2){
        AB.series <- unique(paste(rep(c("A","B"),each = 2),
                                  rep(c("A","B"),2),
                                  sep = ""))
    }

    if(M == 3){
        AB.series <- unique(paste(rep(c("A","B"),each = 4),
                                  rep(rep(c("A","B"),each = 2),2),
                                  rep(c("A","B"), 4),
                                  sep = ""))
    }
        
    targets <- foreach(g1 = AB.series, .combine = "c") %do% {
        targets <- foreach(g2 = AB.series, .combine = "c") %do% {

            ## position 1
            p1 <- (substr(target.geno,1,1) == "H" &  (substr(g1,1,1) != substr(g2,1,1))) |
                (substr(target.geno,1,1) == "A" &  (substr(g1,1,1) == "A" & substr(g2,1,1) == "A")) |
                    (substr(target.geno,1,1) == "B" &  (substr(g1,1,1) == "B" & substr(g2,1,1) == "B"))
    
            ## position 2
            p2 <- (substr(target.geno,2,2) == "H" &  (substr(g1,2,2) != substr(g2,2,2))) |
                (substr(target.geno,2,2) == "A" &  (substr(g1,2,2) == "A" & substr(g2,2,2) == "A")) |
                    (substr(target.geno,2,2) == "B" &  (substr(g1,2,2) == "B" & substr(g2,2,2) == "B"))

            ## position 3
            if(M == 3)
                p3 <- (substr(target.geno,3,3) == "H" &  (substr(g1,3,3) != substr(g2,3,3))) |
                    (substr(target.geno,3,3) == "A" &  (substr(g1,3,3) == "A" & substr(g2,3,3) == "A")) |
                        (substr(target.geno,3,3) == "B" &  (substr(g1,3,3) == "B" & substr(g2,3,3) == "B"))
            
            if(p1 & p2 & ifelse(M == 3, p3, TRUE))
                return(paste(g1,g2,sep = "-"))
        }
        return(targets)
    }

    return(targets)
}
