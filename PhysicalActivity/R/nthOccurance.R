`nthOccurance` <-
function(dataVct, value, nth = NA, reverse = FALSE)
{
    loc = c()
    if(reverse){
        dataVct = rev(dataVct)
    }

    if(is.na(value)){
        value = "NA"
        dataVct[is.na(dataVct)] = "NA"
    }

    temp = 1:length(dataVct)
    if(length(nth)==1){
        if( is.na(nth)){
            loc = temp[match(dataVct, value, nomatch = 0)==1]
        }else{
            loc = temp[match(dataVct, value, nomatch = 0)==1][nth]
        }
    }else{
        loc = temp[match(dataVct, value, nomatch = 0)==1][nth]
    }

    if(reverse){ 
        loc = length(dataVct) - loc +1
    }

    if(sum(is.na(loc)) == length(loc)){
        loc = 0
    }

    return(loc)
}

