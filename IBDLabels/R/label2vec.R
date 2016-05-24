label2vec <-
function(label, ngam){
    
    ## checking
    if( length( label )>1) stop( "Label must be scalar" )
    if( label > maxlabel( ngam )) stop( "Label is too large for given number gametes \n" )
    
    ##initialise
    quo = label
    flag = 0
    vec = numeric(ngam)

    for( i in (ngam-1):0 ){
        if (quo==0) break
        
        ##find remainder
        rem = quo %% i
        
        ##continue if remainder is non-zero, or quotient < maxlabel(i+1)
        if (rem !=0 | quo < maxlabel(i+1) ){ 
            vec[i+1] = rem+1                    ## add (rem+1) to (i+1) gamete
            quo = quo %/% i                     ## update quotient
        }else{
            flag =1; break
        }

       ## cat( "i=", i, "quo=", quo, "rem=", rem, "flag=", flag,
       ## "max=", maxlabel(i+1), "vec=", vec, "\n" )
  }

    ## cat( "break out on ", i, "quo = ", quo, "rem=", rem, "flag=", flag, "vec=", vec, "\n") 

    ## Gives warning on e.g. label2vec( 2,5)  but not label2vec( 43, 5)
    ##if( rem >1 ) warning( "Label gives invalid state vector", "\n" ) 
    
    ##fill rest of s
    if(flag == 0) vec[1:(i+1)] = 1             ## never broke out: put in ones
    if(flag == 1) vec[1:(i+1)] = 1:(i+1)       ## broke out: put in max values


    if( !all( vec==fgl2vec(vec) ) ){ ## INVALID STATE ###
        vec = rep( NA, ngam )
    }
    return(vec)
}

######### OLD VERSION ##########################################
## Updated to give warnings and errors for invalid output/inputs

## label2vec <-
## function(label, ngam){

##     ##initialise
##     cc = label ##(lefotver quotient)
##     flag = 0
##     vec = numeric(ngam)
    
##     for( i in (ngam-1):0 ){
##         ##f cc = 0 break (flag=0)
##         if (cc==0) break 
##         ##find gg (remainder)
##         gg = cc - (cc %/% i)*i
##         ##continue if non-zero or cc < CDgg
##         if (gg !=0 | cc< maxlabel(i+1) ){ 
##             vec[i+1] = gg+1
##             cc = cc %/% i
##         }else{
##             flag =1; break
##         }
##         i=i-1 
##   }
##     ##fill rest of s
##     if(flag == 0) vec[1:(i+1)] = 1
##     if(flag == 1) vec[1:(i+1)] = 1:(i+1)
##     return(vec)
## }
