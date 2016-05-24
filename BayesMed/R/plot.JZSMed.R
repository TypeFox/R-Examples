plot.JZSMed <-
  function(x,...){
    par(mai=c(0,0,0.4,0))
    
    # plot mediation figure
    plot(1:3,bty="n",type="n",axes=F,xlab="",ylab="",main="",...)
    
    points(2,2.8,pch="M",cex=2)
    points(1.5,1.5,pch="X",cex=2)
    points(2.5,1.5,pch="Y",cex=2)
    
    arrows(1.55,1.6,1.95,2.65,lwd=2,length=.15)
    arrows(2.05,2.65,2.45,1.6,lwd=2,length=.15)
    arrows(1.55,1.5,2.45,1.5,lwd=2,length=.15)
    
    # add values for paths
    
    text(1.5,2.2,substitute("p("*{alpha != 0}*"|D) = "*p_a,list(p_a=round(x[rownames(x)%in%"alpha",colnames(x)%in%"PostProb"],2))))
    
    if(any(colnames(x)%in%"Estimate")){
      text(1.5,2.1,substitute({hat(alpha)==a},list(a=round(x[rownames(x)%in%"alpha",colnames(x)%in%"Estimate"],2))))
    }
    
    text(2.5,2.2,substitute("p("*beta != 0*"|D) = "*p_b,list(p_b=round(x[rownames(x)%in%"beta",colnames(x)%in%"PostProb"],2))))
    
    if(any(colnames(x)%in%"Estimate")){
      text(2.5,2.1,substitute({hat(beta)==b},list(b=round(x[rownames(x)%in%"beta",colnames(x)%in%"Estimate"],2))))
    }
    
    text(2,1.3,substitute("p("*tau*{symbol("\242")} != 0*"|D) = "*p_t,list(p_t=round(x[rownames(x)%in%"tau_prime",colnames(x)%in%"PostProb"],2))))
    
    if(any(colnames(x)%in%"Estimate")){
      text(2,1.2,substitute(hat(tau*{symbol("\242")}) == t,list(t=round(x[rownames(x)%in%"tau_prime",colnames(x)%in%"Estimate"],2))))
    }
    
    # add Bayes factors for mediation and full mediation
    legend("topleft",legend=paste("BF_Mediation",round(x[rownames(x)%in%"Mediation (alpha*beta)",colnames(x)%in%"BF"],2),sep=": "))
    
    par(mai=c(1.360000, 1.093333, 1.093333, 0.560000)) # restore default settings
    
  }