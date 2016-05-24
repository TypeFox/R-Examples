`getDistTable` <-
function(S.matrix)
{   
    freqTable = as.data.frame(table(apply(S.matrix, 1, function(x){paste(x, collapse=",")}))); 
    uniq = as.vector(freqTable[,1]);
    freq = as.vector(freqTable[,2],mode='numeric');
    freq = freq/sum(freq);
    
    ## Sort the joint sufficient statistics by the observed frequencies
    sorted = sort(freq, index.return=TRUE);
    uniq = uniq[sorted$ix];
    freq = freq[sorted$ix];
    
    if(length(uniq) == 1)
    {
        distTable = cbind(S.matrix[1,],1);
    }
    else
    {
        distTable = data.frame(matrix(t(unlist(strsplit(uniq,","))),nrow=length(uniq), byrow=T));
        distTable = cbind(distTable,freq);
        distTable = apply(distTable,2,as.numeric);
    }
    
    return(distTable);
}

