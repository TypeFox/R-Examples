label2lex <-
function(label, ngam){

    ## label can be a vector
    
    all.lex <- allLex( ngam ) 
    lex <- all.lex[ label+1 ] ## need to +1 as label starts at zero 
        
    return( lex )  
}
