vec2label <-
function(vec ){

    ## always renumber
    vec = fgl2vec( vec )
    ## calculate valid label
    n = length(vec)
    label = 0
    for( i in 0:(n-1)){
        label = label*i + (vec[i+1] - 1)
    }
    return(label)
}
