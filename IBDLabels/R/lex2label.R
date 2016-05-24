lex2label <-
function( lex, ngam ){

    ## lex can be a vector
    
    all.lex <- allLex( ngam )

    label <- numeric( length( lex ))
    for(i in 1:length(lex))  label[i] <- which( all.lex == lex[i] )
   
    label <- label -1 ## indexing needs to start at 0
                     
    return( label ) 
}
