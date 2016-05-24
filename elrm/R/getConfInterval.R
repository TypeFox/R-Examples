`getConfInterval` <-
function(distTable,S.observed, alpha)
{      
    lower.bound = 'NA';
    upper.bound = 'NA';
    
    lower.func = function(beta.coeff)
    {
        freq = adjustProb(beta.coeff,distTable[,1],distTable[,2]);
        
        prob = sum(freq[distTable[,1]>=S.observed]);
        
        return(prob-(alpha/2));
    }
    
    upper.func = function(beta.coeff)
    {
        freq = adjustProb(beta.coeff,distTable[,1],distTable[,2]);
        
        prob = sum(freq[distTable[,1]<=S.observed]);
        
        return(prob-(alpha/2));
    }
    
    if(S.observed == min(distTable[,1]))
    {
        lower.bound = -Inf;
    }
    else
    {
        lower = 0;

        while(lower.func(lower) > 0)
        {
            lower = lower - 1;
        }

        upper = 0;

        while(lower.func(upper) < 0)
        {
            lower = upper;
            upper = upper + 1;
        }

	  if(lower == upper)
	  {
		lower.bound = lower;
	  }
	  else
	  {
        	lower.bound = uniroot(f=lower.func,interval=c(lower,upper))$root;
        }
    }
    
    if(S.observed == max(distTable[,1]))
    {
        upper.bound = Inf;
    }
    else
    {
        lower = 0;

        while(upper.func(lower) < 0)
        {
            lower = lower - 1;        
        }

        upper = 0;

        while(upper.func(upper) > 0)
        {
            lower = upper;
            upper = upper + 1;        
        }

        if(lower == upper)
	  {
		upper.bound = upper;
	  }
	  else
	  {
     	      upper.bound = uniroot(f=upper.func,interval=c(lower,upper))$root;
        }
    }

    confidence.interval = c(lower.bound,upper.bound);
    
    return(confidence.interval);
}

