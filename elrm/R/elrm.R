`elrm` <-
function(formula, interest, r=4, iter=1000, dataset, burnIn=0, alpha=0.05)
{   
    if(iter %% 1 != 0)
    {
        stop("'iter' must be an integer");
    }

    if(iter < 1000)
    {
        stop("'iter' must be atleast 1,000");
    }

    if(burnIn %% 1 != 0)
    {
        stop("'burnIn' must be an integer");
    }

    if(iter < burnIn + 1000)
    {
        stop("'burnIn' must be atleast 1,000 iterations smaller than 'iter'");
    }

    info = getDesignMatrix(formula,interest,dataset=dataset);
    design.matrix = info$design.matrix;
    yCol = info$yCol;
    mCol = info$mCol;
    zCols = info$zCols;
    wCols = info$wCols;
    
    if((r %% 2) != 0)
    {
        stop("'r' must be an even number");
    }
    
    if(r > length(as.vector(design.matrix[,yCol])))
    {
        msg = paste("'r' must be smaller than or equal to",length(as.vector(design.matrix[,yCol])),"(length of the response vector)");
        stop(msg);
    }
    
    out.filename = paste("elrmout",round(runif(1,0,10000),0),".txt",sep="");
    tempdata.filename = paste("elrmtempdata",round(runif(1,0,10000),0),".txt",sep="");
    write.table(design.matrix,tempdata.filename,quote=FALSE,row.names=FALSE,col.names=FALSE);
    
    Z.matrix = matrix(nrow = nrow(design.matrix), ncol = length(zCols));
    Z.matrix = as.matrix(design.matrix[,zCols]);
    observed.response = as.matrix(design.matrix[,yCol]);
    
    S.matrix = NULL;
    
    message("Generating the Markov chain ...");
    
    t1 = format(Sys.time(),"");
    
    sample.size = max(floor(iter*0.05),1);
    
    last.sample = observed.response;
    
    k = 0;
    
    while(k < iter)
    {
        progress((k/iter)*100);
        
        Sys.sleep(0.10);

        design.matrix[,yCol]=last.sample;
        write.table(design.matrix,tempdata.filename,quote=FALSE,row.names=FALSE,col.names=FALSE);
        
        if(k == 0)
	  {
        	.C("MCMC", as.integer(yCol), as.integer(mCol), as.integer(wCols), as.integer(length(wCols)), as.integer(zCols), as.integer(length(zCols)), as.integer(r), as.character(tempdata.filename), as.character(out.filename), as.integer(min(sample.size,iter-k)), as.integer(1), PACKAGE="elrm");
        }
        else
        {
        	.C("MCMC", as.integer(yCol), as.integer(mCol), as.integer(wCols), as.integer(length(wCols)), as.integer(zCols), as.integer(length(zCols)), as.integer(r), as.character(tempdata.filename), as.character(out.filename), as.integer(min(sample.size,iter-k)), as.integer(0), PACKAGE="elrm");
        }

        mc = matrix(scan(out.filename,what=numeric(),skip=0,n=min(sample.size,iter-k)*nrow(dataset),sep="\t",quiet=T),nrow=min(sample.size,iter-k),ncol=nrow(dataset),byrow=T);
    
        S.temp =  mc %*% Z.matrix
        k = k + sample.size;
        S.matrix = rbind(S.matrix,S.temp);
        
        last.sample = mc[nrow(mc),];
    }
    
    S.matrix = round(S.matrix,6);
    
    progress((k/iter)*100);
    
    cat('\n');
        
    Sys.sleep(0.10);
    
    unlink(tempdata.filename);
    unlink(out.filename);
    
    t2 = format(Sys.time(),"");
    
    dif1 = difftime(t2,t1,units="auto");
    message(paste("Generation of the Markov Chain required",  round(dif1,4), attr(dif1,"units")));
    
    message("Conducting inference ...");
    
    t3 = format(Sys.time(),"");
    
    S.observed = round(as.vector(t(Z.matrix)%*%observed.response,mode='numeric'),6);
    names(S.observed) = names(design.matrix)[zCols];
    
    S.begin  = as.matrix(S.matrix[0:burnIn,]);

    if(burnIn == 1)
    {
        S.begin = as.matrix(t(S.begin));
    }

    S.matrix = as.matrix(S.matrix[(burnIn+1):iter,]);
    
    marginal = getMarginalInf(S.matrix, S.observed, alpha);

    if(length(zCols) > 1)
    {
        joint = getJointInf(S.matrix, S.observed);
        
        t4 = format(Sys.time(),"");
        
        dif2 = difftime(t4,t3,units="auto");
        message(paste("Inference required", round(dif2,4), attr(dif2,"units")));
        
        coeffs = as.vector(c(NA,marginal$estimate),mode='numeric');
        names(coeffs) = c("joint",names(marginal$estimate));
        pvals = as.vector(c(joint$pvalue,marginal$pvalue),mode='numeric');
        pvals.se = as.vector(c(joint$pvalue.se,marginal$pvalue.se),mode='numeric');
        names(pvals) = names(coeffs);
        names(pvals.se) = names(coeffs);
        
        marginal$distribution[["joint"]] = joint$distribution;
        
        sizes = c(joint$mc.size,marginal$mc.size);
        names(sizes) = c("joint",names(marginal$mc.size));

        S.matrix = data.frame(rbind(S.begin,S.matrix));
        colnames(S.matrix) = names(design.matrix)[zCols];
        
        ## S3 definition
        object <- list(coeffs=coeffs,coeffs.ci=marginal$estimate.ci,p.values=pvals,p.values.se=pvals.se,mc=mcmc(S.matrix),mc.size=sizes,obs.suff.stat=S.observed,distribution=marginal$distribution,call.history=list(match.call()),dataset=dataset,last=last.sample,r=r,ci.level=(1-alpha)*100);
        class(object) <- "elrm";
    }
    else
    {
        t4 = format(Sys.time(),"");
        
        dif2 = difftime(t4,t3,units="auto");
        message(paste("Inference required", round(dif2,4), attr(dif2,"units")));

        S.matrix = data.frame(rbind(S.begin,S.matrix));
        colnames(S.matrix) = names(design.matrix)[zCols];
        
        ## S3 definition
        object <- list(coeffs=marginal$estimate,coeffs.ci=marginal$estimate.ci,p.values=marginal$pvalue,p.values.se=marginal$pvalue.se,mc=mcmc(S.matrix),mc.size=marginal$mc.size,obs.suff.stat=S.observed,distribution=marginal$distribution,call.history=list(match.call()),dataset=dataset,last=last.sample,r=r,ci.level=(1-alpha)*100);
        class(object) <- "elrm";
    }

    return(object);
}

