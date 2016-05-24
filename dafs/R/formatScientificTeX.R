formatScientificTeX  = function(x, width, digits){
    exponent = floor(log10(x))

    fmtString = paste('%',width,'.',digits,'f',sep='')

    str1 = sprintf(fmtString, 10^(-exponent)*x)
    str2 = paste('\\\\times 10^{',exponent,'}',sep='')

    return(paste(str1,str2,sep=''))
}

fmtST = function(x, width, digits)
    return(formatScientificTeX(x, width, digits))
