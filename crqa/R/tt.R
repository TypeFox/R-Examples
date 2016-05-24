## extract vertical lines from a recurrence plot,
## calculate laminarity and trapping time.
## part of this code was borrowed from the original
## tt.m function of the crptool written by Norbert Marwan

## build a random x matrix
# r = 100; c = 100
# x = round(matrix(runif(r*c), r, c))
# whiteline = FALSE
# minvertline = 2
# ans = tt(x, minvertline, whiteline)
# x = rbind(c(1,0,1,1,1,1), c(1,1,0,0,1,1), c(0,1,0,1,0,1))

tt <- function(x, minvertline, whiteline){
    # require(Matrix)
    

##########################################
    ## for black vertical lines
    xb = x ## just copy over x in .m use of double
    ## pad the list row with zeros
    
    xb = rBind(xb, rep(0, ncol(x)), deparse.level = 0) #xb(end+1,:) = 0;
    xb = as.vector(xb) ## make it a vector
    z = diff(xb)
    z0 = which(z == 1)  ## begin of black sequence
    z1 = which(z == -1) ## end of black sequence
    
    ## measure the length of black lines
    ## do some padding on the series
    
    if (z0[1] > z1[1]){ 
        z0 = c(0, z0)  # add one zero on the top = z0(1:end); z0(1,1) = 0
    }
    
    if (length(z0) > length(z1)){ 
        z0 = z0[-length(z0)] #(end)=[];
    }
    
    
    t = sort(z1-z0);
    t1 = t[which(t-1 > 0)]

    TT = mean(t1) ## trapping time
    
    ## calculate laminarity: the amount of laminar phases in the system.
    nrbln = sum(x) ## nr. of recurrent points
    
    ## total number of vertical lines
    ## including those below the vertline threshold.
    mnvert = which(t1 < minvertline)

    if (length(mnvert) != 0){ ## there are vertical lines smaller than minimum    
        t1tr = t1[-mnvert]
    } else {
        t1tr = t1
    }

    ## else just return t1tr

    if(length(t1tr) > 0){ ## there are vertical lines
        lam = (sum(t1tr)/nrbln)*100
    } else {
        lam = 0 }
    
    
####################################
    ## for white vertical lines
    
    if (whiteline == TRUE){
        
        xw = as.matrix(x) ## copy over x another time
        ind = which(xw > 0, arr.ind = TRUE) ## extract the indeces with rec. points
        
        r = ind[,1]; c = ind[,2]
        
        vertline = tapply(r, c,
            function(x){
                i1 = min(x); i2 = max(x)
                ind = cbind(i1, i2)
                return( ind ) } ) ## take min and max indeces of the matrix
        ## where recurrent points start(min) and end(max)
        
        ## as rec. indeces are repeated across columns address them using
        ## vector operations
    
        colind = as.numeric( names(vertline) )
        matind = matrix(as.numeric( unlist(vertline) ),
            byrow = TRUE, ncol = 2)
        ## put the indeces in a matrix
        matind = rBind(matind, rep(0,ncol(matind))) ## add one column for padding
        reprw = which(diff(matind)[,1] != 0) ## get the length of column
        unind = matind[reprw, ] 
        
        
        for (rp in 1:length(reprw) ){
            
            if (rp == 1){
                
                xw[1:unind[rp,1], colind[1:reprw[rp]] ] = 1
                xw[unind[rp,2]:nrow(xw), colind[1:reprw[rp]] ] = 1
                
            } else {
                
                xw[1:unind[rp,1], colind[reprw[rp -1]:reprw[rp] ] ] =  1            
                xw[unind[rp,2]:nrow(xw), colind[reprw[rp -1]:reprw[rp] ] ] = 1
                
            }
        }
        
        xw = rBind(xw, rep(1, ncol(xw)), deparse.level = 0)    
        
        zw = diff(as.vector(xw)) ;
        
        z0w = which(zw == -1); # begin of white sequence
        z1w = which(zw == 1);  # end of white sequence
        
        
        ## measure the length of white lines
        if (z0w[1] > z1w[1]){
            z0w = z0w[-1]
            if (length(z1w) > length(z0w) ){ 
                z1w = z1w[-length(z1w)]
            }
        }
        
        
        if ( length(z1w) > length(z0w) ){ 
            z0w = c(1,z0w)
        }
        
        tw = sort(z1w-z0w)
        t1w = tw[which(tw-1 > 0)]
        
        tw = tw     
        
    } else {  tw = NA }
    
    return(list (TT = TT, lam = lam, tw = tw, tb = t) )
    
}
