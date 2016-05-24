`adjustProb` <-
function(B.new, S.uniq, S.uniq.freq)
{   
    S.uniq.freq.new = NULL;

    k = range(S.uniq);
    S.uniq = (S.uniq - k[1])/(k[2]-k[1]);

    B.new = B.new;
    
    if(length(B.new) == 1)
    {
        norm.const = as.double(1.0 / t(as.matrix(S.uniq.freq))%*%exp(B.new*as.matrix(S.uniq)));
        
        S.uniq.freq.new = norm.const*as.matrix(S.uniq.freq)*exp(B.new*as.matrix(S.uniq));
    }
    else
    {   
        norm.const = exp(t(B.new)%*%t(as.matrix(S.uniq)));
        norm.const = norm.const%*%as.matrix(S.uniq.freq);
        norm.const = 1.0/norm.const;
        
        S.uniq.freq.new = as.vector(S.uniq.freq*exp(B.new%*%t(as.matrix(S.uniq))));
        
        S.uniq.freq.new = S.uniq.freq.new*norm.const;
    }
    
    return(S.uniq.freq.new);
}

