## written in R by Moreno I. Coco, 2013, (moreno.cocoi@gmail.com)
## crqa, adapted from a Matlab code developed at
## summer school of: Nonlinear Methods for Psychological Science
## organized by the University of Cincinnati, 2012

## next time i check change sd(matrix) with sapply(matrix, sd)

## arguments to pass to crqa:
## ts1, ts2: times series of integers indicating the states
## delay = nr. of lags
## embed = the embedding dimension, i.e., the lag intervals
## rescale = the normalization for the distance; 
##           if 1 (Mean Distance); if 2 (Max Distance)
## radius =  distance to accept two points as recurrent (set it very 
##           small, if the series are categorical in nature)
## normalize = rescale factor for source data; 
##           if 1 (Unit interval); if 2 (z-score) 
## mindiagline = set a minimum diagonal line length
## mindiagline = set a minimum vertical line length

##  whiteline = FALSE # - flag to compute or not white vertical lines
##                    in the recurrence plot. Note, white lines are not
##                    yet used to derive any particular measure
##  recpt = FALSE # - flag to indicate whether the input ts1 is already
##                a recurrence plot

## tw = the size of the Theiler Window, the default is 0
## side = a string indicating whether the recurrence measures
## should be calculated in the "upper" triangle of the matrix
## "lower" triangle of the matrix, on the "whole" matrix

## checkl = a list with four arguments:
##       $do = TRUE|FALSE; normalize (or not) the length of ts
## if $do = TRUE, then the arguments of checkts() needs to be passed.
##    $datatype = (numerical, categorical) - nature of ts
##    $thrshd: number of timepoints we accept two ts to be different
##    $pad: whether the two series have to be padded (TRUE) or chopped

## try below

## ts1 = c(0,0,1,1,0,0)
## ts2 = c(2,2,1,1,2,2)
#  ts1 = c(0, 0, 1, 1, 0, 0, 2, 2, 1, 1)
#  ts2 = c(1,1, 2, 2, 0, 0, 1, 2)
# delay = 1; embed = 1; rescale = 1; radius = 0.001;
# normalize = 0; mindiagline = 2; minvertline = 2;
# tw = 0; whiteline = FALSE; recpt = FALSE; side = "both"
# checkl = list(do = TRUE, thrshd = 3, datatype = "categorical",
#     pad = TRUE)

#crqa(ts2, ts1, delay, embed, rescale, radius, normalize, mindiagline, minvertline, tw, whiteline, recpt, side, checkl)

.packageName <- 'crqa'

crqa <- function(ts1, ts2, delay, embed, rescale,
                 radius, normalize, mindiagline, minvertline,
                 tw = 0, whiteline = F, recpt = F, side = "both",
                 checkl = list(do = F)){

    ## passing a few default above
    # if( missing(tw) ){ tw = 0} -> not working

    v11 = v21 = NULL ## stupid initializations to please CRAN
    
    # require("fields")  ## to compute the Euclidean distance matrix
    # require("Matrix")  ## to manipulate sparse matrices 


    ## check if the input is a recurrence plot 
    if (recpt == FALSE){
    
        ts1 = as.vector(as.matrix(ts1)) ## make sure data is a vector
        ts2 = as.vector(as.matrix(ts2))
    
        if (is.matrix(ts1)){ stop("Your data must consist of a single column of data.")}  
        if (is.matrix(ts2)){ stop("Your data must consist of a single column of data.")}      

        
        ## check the length of the sequences and decide if they have
        ## to be normalized to the same length.
        
        if (checkl$do == TRUE){
            
            tsnorm = checkts(ts2, ts1, checkl$datatype,
                checkl$thrshd, checkl$pad)

            if (tsnorm[[2]] == FALSE){
                stop("Time-series difference longer than threshold. Increase threshold, or set checkl$do = FALSE avoiding normalization of ts")
                 } else {
                     ts1 = tsnorm[[1]][,1]
                     ts2 = tsnorm[[1]][,2]
                 }
        
        }

        ## check that the length of the series is not shorter than the phase embed*delay

        if (length(ts1) < embed*delay){ stop("Phase-space (embed*delay) longer than ts1")}  
        if (length(ts2) < embed*delay){ stop("Phase-space (embed*delay) longer than ts2")}      
    
        ##rescale the data if really necessary
        
        if (normalize > 0){
            switch (normalize,
                    {1
                     ## unit-interval
                     ts1 = (ts1 - min(ts1));
                     ts1 = ts1 / max(ts1);
                     ts2 = (ts2 - min(ts2));
                     ts2 = ts2 / max(ts2);},
                    
                    {2                     
                     ## z-score                    
                     ts1 = (ts1 - mean(ts1))/sd(ts1)                    
                     ts2 = (ts2 - mean(ts2))/sd(ts2)
                 }
                    )
        }
        
        ## start to compute recurrence
        ## do it twice for the two series
        
        for (loop in 1:embed){
            vectorstart = (loop-1) * delay + 1;
            vectorend = length(ts1) - ( (embed-loop)* delay);
            assign(paste("v1", loop, sep =""), ts1[vectorstart:vectorend]);
        }
        
        for (loop in 1:embed){
            vectorstart = (loop-1) * delay + 1;
            vectorend = length(ts2) - ( (embed-loop)* delay);
            assign(paste("v2", loop, sep =""), ts2[vectorstart:vectorend]);
        }
        
        ## Create matrix from vectors to use for distance matrix calcs

        dimts1 = dimts2 = vector() ## initialize dims of embedding 
        
        for (loop in 1:embed){
            if (loop == 1){ dimts1 = v11 }
            else{
                eval(
                    parse(
                        text = paste("dimts1 = cbind(dimts1,",
                            paste( "v1", loop, sep = ""),
                            ", deparse.level = 0)", sep = "")
                        )
                    )
            }
        }
    
        for (loop in 1:embed){
            if (loop == 1){ dimts2 = v21 }
            else{
                eval(
                    parse(
                        text = paste("dimts2 = cbind(dimts2,",
                            paste( "v2", loop, sep = ""),
                            ", deparse.level = 0)", sep = "")
                        )
                    )
            }
        }
        
        
        ## Compute euclidean distance matrix
        v1l = length(v11)
        v2l = length(v21)

        
        ## just to have the length of matrix saved
        dm = rdist(dimts1,dimts2);
        
        ## Find indeces of the distance matrix that fall
        ## within prescribed radius.
        if (rescale > 0){
            switch(rescale,
                   {1  ## Create a distance matrix that is re-scaled
                    ## to the mean distance
                    
                    rescaledist = mean( dm )    
                    dmrescale=(dm/rescaledist)*100},
                   
                   {2  ## Create a distance matrix that is re-scaled
                    ## to the max distance
                    
                    rescaledist = max(dm);
                    dmrescale = (dm/rescaledist)*100}
                   )
        } else { dmrescale = dm }
        ## Compute recurrence matrix
        
        ind = which(dmrescale <= radius, arr.ind = TRUE);
        r = ind[,1]; c = ind[,2]

    } else { ## take as input an RP directly

        ## as usual R needs fiddly code to make sure about identify of data
        ts1 = matrix(as.logical(ts1), ncol = ncol(ts1))
        v1l = nrow(ts1); v2l = ncol(ts1)        

        ## matrix needs to be logical
        ind = which(ts1 > 0, arr.ind = TRUE)

        ## just a trick to reduce the number of lines
        ## of the code
        r = ind[,1]; c = ind[,2]
        
    }
    
    if (length(r) != 0 & length(c) != 0){ ##avoid cases with no recurrence

        S = sparseMatrix(r, c, dims = c(v1l, v2l))
        ## this is the recurrent plot
        ## transpose it to make identical to Marwan
        S = t(S)

        ## apply the theiler argument here to recurrence matrix
        ## Marwan blanks out the recurrence along the diag
        S = theiler(S, tw)

        if (side == "upper"){
            ## if only the upper side is of interest
            ## it blanks out the lowest part
            S = as.matrix(S)
            S[lower.tri(S, diag = TRUE)] = 0
            S = Matrix(S, sparse = TRUE)
        }

        if (side == "lower"){
            ## viceversa
            S = as.matrix(S)
            S[upper.tri(S, diag = TRUE)] = 0
            S = Matrix(S, sparse = TRUE)
        }

        if (side == "both"){
            ## just keep it as is.
            S = S}
        
        spdiagonalize = spdiags(S) ##  spdiags should have decent speed 
        B = spdiagonalize$B
        
        ##calculate percentage recurrence by taking all non-zeros
        
        numrecurs = length(which(B == TRUE));
        percentrecurs = (numrecurs/((v1l*v2l)))*100;
        
####################################################################
####################################################################
        
        ## Computing the line counts

        ## This section finds the index of the zeros in the matrix B,
        ## which contains the diagonals of one triangle of the
        ## recurrence matrix (the identity line excluded).

        ## The find command indexes the matrix sequentially
        ## from 1 to the total number of elements.
        ## The element numbers for a 2X2 matrix would be [1 3; 2 4].
        ## You get a hit for every zero. If you take the difference
        ## of the resulting vector, minus 1, it yields the length of an
        ## interceding vector of ones, a line. Here is an e.g.
        ## using a row vector rather than a col. vector, since it types
        ## easier: B=[0 1 1 1 0], a line of length 3.
        ## find( B == 0 ) yields [1 5], diff( [1 5] ) -1 = 3,
        ## the line length.
        ## So this solution finds line lengths in the interior of
        ## the B matrix, BUT fails if a line butts up against either
        ## edge of the B matrix, e.g. say  B = [0 1 1 1 1],
        ## which( B == 0) returns a 1, and you miss the line of length 4.
        ## A solution is to "bracket" B with a row of zeros at each
        ## top and bottom.

        ## Bracket B with zeros
        if (is.vector(B)) {
            false = rep(FALSE, length(B)) ##cases where B is a vector
            B = rbind(false, B, false, deparse.level = 0)
        } else  {
            false = rep(FALSE, ncol(B))
            B = as.matrix(B)
            ## need to transform the sparseMat into normal to bracket it
            B = rbind(false, B, false, deparse.level = 0)
        }
        
        ## Get list of line lengths, sorted from largest to smallest
        diaglines = sort( diff(which(B == FALSE) ) -1, decreasing = TRUE)
        
        ## Delete line counts less than the minimum diagonal.
        diaglines = diaglines[-which(diaglines < mindiagline)]
        ## diaglines(diaglines>200)=[]; # Can define a maximum line length too.

        ## exlude the rare cases where there are no diaglines
        
        if(length(diaglines) != 0){
            
            numdiaglines = length(diaglines) ## extract the length of diag
            maxline = max(diaglines)
            meanline = mean(diaglines)
            
            tabled = as.data.frame(table(diaglines))
            
            total = sum(tabled$Freq)       
            p = tabled$Freq/total
            
            ##remove zero probability..it should not be necessary
            del = which(p == 0 )
            if (length(del) > 0) {
                p = p[-del]
            }
            
            ## entropy log2, and relative entropy divided by max
            entropy = - sum(p*log(p))    
            relEntropy = entropy/(-1*log(1/nrow(tabled)))

            ## entropy/max entropy: comparable across contexts and conditions.
        
            pdeter = sum(diaglines)/numrecurs*100
            ## percent determinism: the predictability of the dynamical system 
        
            ## calculate laminarity and trapping time
            restt = tt(S, minvertline, whiteline)            
            lam = restt$lam; TT = restt$TT
                    
        } else {
            
            numdiaglines = 0; maxline = 0; pdeter = NA;
            entropy = NA; relEntropy = NA; meanline = 0
            lam = 0; TT = 0; RP = NA
        }
        
        results = list(RR = percentrecurs, DET = pdeter, 
            NRLINE = numdiaglines, maxL = maxline, 
            L = meanline, ENTR = entropy, 
            rENTR = relEntropy,
            LAM = lam, TT = TT, RP = S)
        
    } else { # print (paste ("No recurrence found") )
             results = list(RR = 0, DET = NA, NRLINE = 0,
                 maxL = 0, L = 0, ENTR = NA, rENTR = NA,
                 LAM = NA, TT = NA, RP = NA)}  
    
    return (results)
    
}
