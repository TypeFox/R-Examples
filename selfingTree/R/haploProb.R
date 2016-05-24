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

haploProb.2M <- function(gam1,gam2,target){


    ## fixed positions
    fixed <- c((substr(target,1,1) == substr(gam1,1,1) & (substr(target,1,1) == substr(gam2,1,1))),
               (substr(target,2,2) == substr(gam1,2,2) & (substr(target,2,2) == substr(gam2,2,2))))
    
    
    {
        ## target not present in one gamete
        if((substr(target,1,1) != substr(gam1,1,1) & (substr(target,1,1) != substr(gam2,1,1))) |
           (substr(target,2,2) != substr(gam1,2,2) & (substr(target,2,2) != substr(gam2,2,2))))
            haplo.prob <- "0"
        ## target present at all loci in both gametes
        else if(target == gam1 & target == gam2)
            haplo.prob <- "1"
        ## one loci fixed for target
        else if(sum(fixed) == 1 & (target == gam1 | target == gam2))
            haplo.prob <- "1/2 * ((1 - z) + z)"
        ## none fixed:
        ## one gamete equal to target, the other completely different
        else if((substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2))  &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2)) |
                ##
                (substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2))  &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2)))
            haplo.prob <- "1/2 * (1 - z)"
        else if((substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2))  &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2))  |
                ##
                (substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2))  &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2)))
            haplo.prob <- "1/2 * z"
    }
    
    return(haplo.prob)
}




haploProb.3M <- function(gam1,gam2,target){
    ## fixed positions
    fixed <- c((substr(target,1,1) == substr(gam1,1,1) & (substr(target,1,1) == substr(gam2,1,1))),
               (substr(target,2,2) == substr(gam1,2,2) & (substr(target,2,2) == substr(gam2,2,2))),
               (substr(target,3,3) == substr(gam1,3,3) & (substr(target,3,3) == substr(gam2,3,3))))

            
    {
        ## target not present in one gamete
        if((substr(target,1,1) != substr(gam1,1,1) & (substr(target,1,1) != substr(gam2,1,1))) |
           (substr(target,2,2) != substr(gam1,2,2) & (substr(target,2,2) != substr(gam2,2,2))) |
           (substr(target,3,3) != substr(gam1,3,3) & (substr(target,3,3) != substr(gam2,3,3))))
            haplo.prob <- "0"
        ## target present at all loci in both gametes
        else if(target == gam1 & target == gam2)
            haplo.prob <- "1"
        ## two loci fixed for target
        else if(sum(fixed) == 2)
            haplo.prob <- "1/2"
        ## one loci fixed for target
        else if(sum(fixed) == 1 & (target == gam1 | target == gam2)){
            if(fixed[1])
                haplo.prob <- "1/2 * ((1 - x) * (1 - y) + x * (1 - y))"
            else if(fixed[2])
                haplo.prob <- "1/2 * ((1 - x) * (1 - y) + x * y)"
            else if(fixed[3])
                haplo.prob <- "1/2 * ((1 - x) * (1 - y) + (1 - x) * y)"
        }
        else if(sum(fixed) == 1){
            if(fixed[1])
                haplo.prob <- "1/2 * (x * y + (1- x) * y)"
            else if(fixed[2])
                haplo.prob <- "1/2 * (x * (1 - y) + (1 - x) * y)"
            else if(fixed[3])
                haplo.prob <- "1/2 * (x * y + x * (1 - y))"
        }
        ## none fixed:
        ## one gamete equal to target, the other completely different
        else if((substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2) &
                 substr(target,3,3) == substr(gam1,3,3)) &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2) &
                 substr(target,3,3) != substr(gam2,3,3)) |
                ##
                (substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2) &
                 substr(target,3,3) != substr(gam1,3,3)) &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2) &
                 substr(target,3,3) == substr(gam2,3,3)))
            haplo.prob <- "1/2 * (1 - x) * (1 - y)"
        else if((substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2) &
                 substr(target,3,3) == substr(gam1,3,3)) &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2) &
                 substr(target,3,3) != substr(gam2,3,3)) |
                ##
                (substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2) &
                 substr(target,3,3) != substr(gam1,3,3)) &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2) &
                 substr(target,3,3) == substr(gam2,3,3)))
            haplo.prob <- "1/2 * x * y"
        else if((substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2) &
                 substr(target,3,3) != substr(gam1,3,3)) &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2) &
                 substr(target,3,3) == substr(gam2,3,3)) |
                ##
                (substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2) &
                 substr(target,3,3) == substr(gam1,3,3)) &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2) &
                 substr(target,3,3) != substr(gam2,3,3)))
            haplo.prob <- "1/2 * (1 - x) * y"
        else if((substr(target,1,1) != substr(gam1,1,1) &
                 substr(target,2,2) == substr(gam1,2,2) &
                 substr(target,3,3) == substr(gam1,3,3)) &
                (substr(target,1,1) == substr(gam2,1,1) &
                 substr(target,2,2) != substr(gam2,2,2) &
                 substr(target,3,3) != substr(gam2,3,3)) |
                (substr(target,1,1) == substr(gam1,1,1) &
                 substr(target,2,2) != substr(gam1,2,2) &
                 substr(target,3,3) != substr(gam1,3,3)) &
                (substr(target,1,1) != substr(gam2,1,1) &
                 substr(target,2,2) == substr(gam2,2,2) &
                 substr(target,3,3) == substr(gam2,3,3)))
            haplo.prob <- "1/2 * x * (1 - y)"
    }


    return(haplo.prob)

}
