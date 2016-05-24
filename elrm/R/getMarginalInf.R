`getMarginalInf` <-
function(S.matrix,S.observed,alpha)
{
    distribution = list();
    estimate = c();
    estimate.ci = c();
    pvalue = c();
    pvalue.se = c();
    mc.size = c();

    for(i in 1:ncol(S.matrix))
    {  
        ## GET FULL CONDITIONAL MARGINAL DISTRIBUTION TABLE
        
        if(ncol(S.matrix) == 1)
        {
            S.marginal = S.matrix;
        }
        else
        {
            matches = row.matches(as.vector(S.observed[-i]),as.matrix(S.matrix[,-i]));
            
            if(length(matches) != 0)
            {
                S.marginal = as.matrix(S.matrix[row.matches(as.vector(S.observed[-i]),as.matrix(S.matrix[,-i])),i]);
            }
            else
            {
                S.marginal = matrix(ncol = 0, nrow = 0);
            }
        }
        
        if(nrow(S.marginal) == 0)
        {
            distTable = NULL;
        }
        else
        {
            distTable = getDistTable(S.marginal);
            colnames(distTable) = c(names(S.observed)[i],"freq"); 
        }
        
        distribution[[names(S.observed)[i]]] = distTable;      
        mc.size = c(mc.size,nrow(S.marginal));
        
        index = match(S.observed[i],distTable[,1]);
        
        if(!is.null(distTable))
        {
            if(nrow(distTable) > 1  && length(S.marginal) >= 1000)
            {
                if(!is.na(index))
                {
                    ## GET P-VALUE ESTIMATE

                    observed.freq = distTable[index,2];
                    p.value = sum(distTable[distTable[,2] <= observed.freq,2]);

                    ## GET PARAMETER ESTIMATE CI

                    ci = getConfInterval(distTable,S.observed[i],alpha);

                    ## GET PARAMETER ESTIMATE

                    if(S.observed[i] == min(S.marginal[,1]))
                    {
                        lhs.func = function(beta.coeff)
                        {
                            freq = adjustProb(beta.coeff,distTable[,1],distTable[,2]);
                            prob = sum(freq[distTable[,1]<=S.observed[i]]);
                            return(prob-0.5);
                        }

                        lower = 0;

                        while(lhs.func(lower) < 0)
                        {
                            lower = lower - 1;
                        }
                        
                        upper = 0;

                        while(lhs.func(upper) > 0)
                        {
                            lower = upper;
                            upper = upper + 1;
                        }

                        param.estimate = uniroot(f=lhs.func,interval=c(lower,upper))$root;
                    }
                    else if(S.observed[i] == max(S.marginal[,1]))
                    {
                        rhs.func = function(beta.coeff)
                        {
                            freq = adjustProb(beta.coeff,distTable[,1],distTable[,2]);
                            prob = sum(freq[distTable[,1]>=S.observed[i]]);
                            return(prob-0.5);
                        }

                        lower = 0;

                        while(rhs.func(lower) > 0)
                        {
                            lower = lower - 1;
                        }

                        upper = 0;

                        while(rhs.func(upper) < 0)
                        {
                            lower = upper;
                            upper = upper + 1;
                        }

                        param.estimate = uniroot(f=rhs.func,interval=c(lower,upper))$root;
                    }
                    else
                    {
                        likelihood = function(beta.coeffs)
                        {
                            prob = adjustProb(beta.coeffs,distTable[,1],distTable[,2])[index];
                            return(prob);
                        }

                        lower = ci[1]
                        upper = ci[2]

                        param.estimate = optimize(f=likelihood,lower=lower,upper=upper,maximum=TRUE)$maximum;
                    }

                    ## GET STANDARD ERROR OF P-VALUE ESTIMATE

                    values = S.marginal;

                    for(j in 1:nrow(distTable))
                    {
                        values[values==distTable[j,1]]=distTable[j,2];
                    }

                    values=as.vector(values,mode='numeric');
                    batch.size = batchSize(values);

                    if(batch.size > length(values)/10)
                    {
                        batch.size=sqrt(length(values));
                    }

                    se = sqrt(bm(vals=values,bs=batch.size,g=function(x) { return(as.vector(x <= observed.freq, mode='numeric')); })$var);

                    ## get scaling constant used in adjustProb function
                    k = range(distTable[,1]);
                    k = k[2]-k[1];

                    ## fill the results with the correctly scaled estimates
                    estimate = c(estimate,param.estimate/k);
                    estimate.ci = rbind(estimate.ci,ci/k);
                    pvalue = c(pvalue,p.value);
                    pvalue.se = c(pvalue.se,se);
                }
                else
                {
                    ## observed value of sufficient statistic not sampled
                    estimate = c(estimate,NA);
                    estimate.ci = rbind(estimate.ci,c(NA,NA));
                    pvalue = c(pvalue,0);
                    pvalue.se = c(pvalue.se,0);

                    warning(paste("'",names(S.observed)[i],"'"," observed value of the sufficient statistic was not sampled",sep=""), call.=FALSE);
                }
            }
            else
            {
                ## distribution is degenerate
                estimate = c(estimate,NA);
                estimate.ci = rbind(estimate.ci,c(NA,NA));
                pvalue = c(pvalue,NA);
                pvalue.se = c(pvalue.se,NA);

                if(nrow(distTable) <= 1)
                {
                    warning(paste("'",names(S.observed)[i],"'"," conditional distribution of the sufficient statistic was found to be degenerate",sep=""), call. = FALSE);
                }
                else
                {
                    warning(paste("'",names(S.observed)[i],"'"," extracted sample is too small for inference (less than 1000)",sep=""), call. = FALSE);
                }
            }
        }
        else
        {
            ## distTable is null
            estimate = c(estimate,NA);
            estimate.ci = rbind(estimate.ci,c(NA,NA));
            pvalue = c(pvalue,NA);
            pvalue.se = c(pvalue.se,NA);

            warning(paste("'",names(S.observed)[i],"'"," conditional distribution of the sufficient statistic was found to be degenerate",sep=""), call. = FALSE);
        }
    }
    
    names(estimate)=names(S.observed);
    estimate.ci = as.data.frame(estimate.ci,row.names=names(S.observed));
    names(estimate.ci) = c("lower","upper");
    names(pvalue)=names(S.observed);
    names(pvalue.se)=names(S.observed);
    names(mc.size)=names(S.observed);
    marginal = list(estimate,estimate.ci,pvalue,pvalue.se,mc.size,distribution);
    names(marginal) = c("estimate","estimate.ci","pvalue","pvalue.se","mc.size","distribution");
    return(marginal);
}

