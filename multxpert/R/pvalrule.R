# PValRule function  is a secondary function which generates decision rules for the Bonferroni, Holm, 
# Hommel, Hochberg, fixed-sequence and fallback procedures in hypothesis testing problems with 
# equally weighted null hypotheses 
pvalrule<-function(rawp,alpha,proc)
# RAWP, Vector of raw p-values
# ALPHA, Familywise error rate
# PROC, Procedure name
{
	# Number of p-values
	m<-length(rawp)
	
	cat("Hypothesis testing problem\n\n")
	cat("Familywise error rate: alpha=",alpha,"\n",sep="")
	
	# List of original null hypotheses and raw p-values
	cat("Original null hypotheses: ",sep="")
	for (i in 1:m) cat("H",i," ",sep="")
	cat("(equally weighted)\n")    
	cat("Original raw p-values: ")
	for (i in 1:m) cat("p",i,"=",round(rawp[i],4)," ",sep="")        
	cat("\n")
	
	# List of ordered null hypotheses and raw p-values    
	o<-order(rawp)
	orderp<-rawp[o]        
	cat("Ordered null hypotheses: ",sep="")
	for (i in 1:m) cat("H(",i,") ",sep="")
	cat("\n")
	cat("Ordered raw p-values: ")
	for (i in 1:m) cat("p(",i,")=",round(orderp[i],4)," ",sep="")        
	cat("\n\n")
	
	# Print decision rules for each specified procedure
	proc.to.print <- proc
	for(p in proc.to.print) {
	
		# Bonferroni procedure (parameters needed: m, rawp, alpha)
		if (p=="Bonferroni")
		{
			cat("\n\n")
			cat("Decision rules for the Bonferroni procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Original null hypothesis H",i,sep="")
				level<-alpha/m
				if (rawp[i]<=level) cat(" is rejected since p",i,"<=alpha/",m,"=",round(level,4),"\n\n",sep="")
				else cat(" is accepted since p",i,">alpha/",m,"=",round(level,4),"\n\n",sep="")
			}
		}    
        
        # Fixed-sequence procedure    
        if (p=="Fixed-sequence")
            {
            cat("\n\n")
            cat("Decision rules for the fixed-sequence procedure:\n\n",sep="")        
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")
                cat("Original null hypothesis H",i,sep="")
                if (i==1)
                    {
                    # Current null hypothesis is rejected                
                    if (rawp[1]<=alpha) 
                        {
                        current<-1
                        cat(" is rejected since p1<=alpha=",round(alpha,4),"\n\n",sep="")
                        }
                    # Current null hypothesis is accepted
                    if (rawp[1]>alpha) 
                        {
                        current<-0                 
                        cat(" is accepted since p1>alpha=",round(alpha,4),"\n\n",sep="")
                        }
                    }
                if (i>1)
                    {
                    # Preceding null hypothesis
                    prec<-current
                    # Current null hypothesis is rejected                
                    if (rawp[i]<=alpha & prec==1) current<-1
                    # Current null hypothesis is accepted
                    if (rawp[i]>alpha | prec==0) current<-0                
                    if (rawp[i]<=alpha & prec==1) cat(" is rejected since p",i,"<=alpha=",round(alpha,4)," and the preceding original null hypothesis is rejected\n\n",sep="")
                    if (rawp[i]>alpha & prec==1) cat(" is accepted since p",i,">alpha=",round(alpha,4),"\n\n",sep="")
                    if (prec==0) cat(" is accepted since the preceding original null hypothesis is accepted\n\n",sep="")                		       
                }            
                }
            }
            
        # Fallback procedure    
        if (p=="Fallback")
            {
            cat("\n\n")
            cat("Decision rules for the fallback procedure:\n\n",sep="")
            numer<-1
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")
                cat("Original null hypothesis H",i,sep="")            
                level<-numer*alpha/m
                if (i<m)
                    {                               
                    # Current null hypothesis is rejected                
                    if (rawp[i]<=level) 
                        {
                        cat(" is rejected since p",i,"<=",numer,"*alpha/",m,"=",round(level,4)," and the significance level used in this test is carried over to the next original null hypothesis\n\n",sep="")
                        numer<-numer+1
                        }
                    # Current null hypothesis is accepted
                    if (rawp[i]>level) 
                        {
                        current<-0 
                        cat(" is accepted since p",i,">",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")
                        numer<-1
                        }
                    }
                if (i==m)
                    {
                    if (rawp[i]<=level) cat(" is rejected since p",i,"<=",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")
                    else cat(" is accepted since p",i,">",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")                    
                    }
                }        
            }                    
        
		# Holm procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Holm")
		{
			cat("Decision rules for the Holm procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Step ",i,"\n",sep="")
				cat("Ordered null hypothesis H(",i,") [Original null hypothesis H",o[i],"]",sep="")
				level<-alpha/(m-i+1)
				if (i==1)
				{                
					# Current null hypothesis is rejected                
					if (orderp[1]<=level) current<-1
					# Current null hypothesis is accepted
					if (orderp[1]>level) current<-0                
					if (orderp[1]<=level) cat(" is rejected since p(1)<=alpha/",m,"=",round(level,4),"\n\n",sep="")
					else cat(" is accepted since p(1)>alpha/",m,"=",round(level,4),"\n\n",sep="")
				}
				if (i>1)
				{
					# Preceding null hypothesis
					prec<-current
					# Current null hypothesis is rejected                
					if (orderp[i]<=level & prec==1) current<-1
					# Current null hypothesis is accepted
					if (orderp[i]>level | prec==0) current<-0                
					if (orderp[i]<=level & prec==1) cat(" is rejected since p(",i,")<=alpha/",m-i+1,"=",round(level,4)," and the preceding ordered null hypothesis is rejected\n\n",sep="")
					if (orderp[i]>level & prec==1) cat(" is accepted since p(",i,")>alpha/",m-i+1,"=",round(level,4),"\n\n",sep="")
					if (prec==0) cat(" is accepted since the preceding ordered null hypothesis is accepted\n\n",sep="")                		       
				}
			}
		}  
		
		# Hochberg procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Hochberg")
		{
			cat("Decision rules for the Hochberg procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Step ",i,"\n",sep="")            
				level<-alpha/i
				if (i==1)
				{                
					if (orderp[m]<=level) 
					{
						cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
						cat(" is rejected since p(",m,")<=alpha/1=",round(level,4),"\n\n",sep="")  
						if (m>1)    
						{
							for (j in (m-1):1)
							{
								cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
								cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                        
							}
						}    
						break
					}
					if (orderp[m]>level)
					{
						cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
						cat(" is accepted since p(",m,")>alpha/1=",round(level,4),"\n\n",sep="")
					}
				}
				if (i>1)
				{
					if (orderp[m-i+1]<=level) 
					{
						cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"]",sep="")
						cat(" is rejected since p(",m-i+1,")<=alpha/",i,"=",round(level,4),"\n\n",sep="")
						if (m-i+1>1)
						{
							for (j in (m-i):1)
							{
								cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
								cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                        
							}
						}    
						break    
					}
					if (orderp[m-i+1]>level)
					{
						cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"] ",sep="")
						cat("is accepted since p(",m-i+1,")>alpha/",i,"=",round(level,4),"\n\n",sep="")
					}    
				}
			}
		}         
		
		# Hommel procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Hommel")
		{
            cat("Decision rules for the Hommel procedure:\n\n",sep="")
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")            
                if (i==1)
                    {                
                    if (orderp[m]<=alpha) 
                        {
                        cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
                        cat(" is rejected since p(",m,")<=alpha/1=",round(alpha,4),"\n\n",sep="")  
                        if (m>1)    
                            {
                            for (j in (m-1):1)
                                {
                                cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
                                cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                            }
                            }    
                        break   
                        }
                    if (orderp[m]>alpha)
                        {
                        cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
                        cat(" is accepted since p(",m,")>alpha/1=",round(alpha,4),"\n\n",sep="")
                        }
                    }
                if (i>1)
                    {
                    # All p-values are greater than corresponding significance levels
                    all<-1
                    hyprej<-0
                    jrej<-0
                    irej<-0
                    for (j in 1:i)
                        {
                        if (orderp[m-i+j]<=j*alpha/i) 
                            {
                            if (hyprej==0)
                                {
                                all<-0 
                                hyprej<-m-i+j 
                                jrej<-j 
                                irej<-i
                                }
                            }
                        }
                    if (all==1)
                        {
                        cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"] ",sep="")
                        cat("is accepted since",sep="")
                        for (j in 1:(i-1)) cat(" p(",m-i+j,")>",j,"*alpha/",i,"=",round(j*alpha/i,4)," and",sep="")
                        cat(" p(",m,")>alpha/1=",round(alpha,4),"\n\n",sep="")
                        }
                    if (all==0)
                        {
                        cat("p(",hyprej,")<=",jrej,"*alpha/",irej,"=",round(jrej*alpha/irej,4),sep="")
                        cat(" and thus all remaining hypotheses are tested with alpha/",i-1,"=",round(alpha/(i-1),4),": \n\n",sep="")                      
                        allsub<-1
                        for (j in (m-i+1):1)
                            {                      
                            cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
                            if (orderp[j]<=alpha/(i-1)) 
                                {
                                if (allsub==1)
                                    {
                                    allsub<-0
                                    cat(" is rejected since p(",j,")<=alpha/",i-1,"=",round(alpha/(i-1),4),"\n\n",sep="") 
                                    }
                                else cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")
                                }
                            if (orderp[j]>alpha/(i-1)) cat(" is accepted since p(",j,")>alpha/",i-1,"=",round(alpha/(i-1),4),"\n\n",sep="")                              
                            }
                        break    
                        }
                    }
                }
			
		}         
		
		cat("\n\n")
	}
}
# End of pvalrule