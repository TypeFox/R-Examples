PearsonICA<-function(X, n.comp=0,   row.norm = FALSE, maxit = 200, tol = 1e-04,border.base=c(2.6,4),border.slope=c(0,1) 
    ,verbose = FALSE, w.init = NULL, na.rm=FALSE, whitening.only=FALSE, PCA.only=FALSE) 
{
    dd <- dim(X)
    d <- dd[dd != 1]
    if (length(d) != 2) 
        stop("data must be matrix-conformal")
    X <- if (length(d) != length(dd)) 
        matrix(X, d[1], d[2])
    else as.matrix(X)
    Xna<-is.na(rowSums(X));
    if (na.rm) {
        if (verbose) {cat(sum(Xna)," rows omitted due to missing values.\n")}
        X<-X[!Xna,];
    }else {
        if (sum(Xna)>0) 
        { 
            stop("Pearson-ICA cannot proceed because the input data has missing values. 
            If you wish to omit missing values, call PearsonICA with option na.rm=TRUE") 
        }
    }
    n <- nrow(X)
    p <- ncol(X)
    if (n.comp==0) {n.comp <- min(n, p)}
    if (n.comp > min(n, p)) {
        cat("n.comp is too large\nn.comp set to", min(n, p), 
            "\n")
        n.comp <- min(n, p)
    }
    if (is.null(w.init)) {
        w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
    }
    else {
        if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
            stop("w.init is not a matrix or is the wrong size")
    }
    if (verbose) {cat("Centering\n")}
    Xmu<-apply(X,2,mean,na.rm=TRUE);    
    X <- scale(X, scale = FALSE)
   
    if (row.norm) {
        X <- scale(t(X))
    }
    else {
        X <- t(X)
    }
    if (PCA.only)
    {
        V <- X %*% t(X)/n
        s <- La.svd(V)
        K <- t(s$u)
        K <- matrix(K[1:n.comp, ], n.comp, p)
        X1 <- K %*% X
        w <- K
        S <- w %*% X
        A <- t(w) %*% solve(w %*% t(w))
        if (whitening.only) {warning("Both whitening.only and PCA.only are requested. PCA is performed.")}
        return(list(X = t(X), whitemat = NA, W = t(w), A = t(A), S = t(S),Xmu=Xmu,w.init=NA,maxit=NA,tol=NA,it=NA))
    } else {
        if (verbose) 
            cat("Whitening\n")
        V <- X %*% t(X)/n
        s <- La.svd(V)
        D <- diag(c(1/sqrt(s$d)))
        K <- D %*% t(s$u)
        K <- matrix(K[1:n.comp, ], n.comp, p)
        X1 <- K %*% X
    }
    if (whitening.only)
    {
        w <- K
        S <- w %*% X
        A <- t(w) %*% solve(w %*% t(w))
        return(list(X = t(X), whitemat = t(K), W = t(w), A = t(A), S = t(S),Xmu=Xmu,w.init=NA,maxit=NA,tol=NA,it=NA))  
    } else {
        alist <- fasticapearson(X1, n.comp, tol = tol, border.base=border.base,border.slope=border.slope,maxit = maxit, 
            verbose = verbose,w.init = w.init)
        a<-alist$a
        w <- a %*% K
        S <- w %*% X
        A <- t(w) %*% solve(w %*% t(w))
    return(list(X = t(X), whitemat = t(K), W = t(w), A = t(A), S = t(S),Xmu=Xmu,w.init=w.init,maxit=maxit,tol=tol,it=alist$it))
    }
}

fasticapearson<-function (X, n.comp, tol,border.base,border.slope, maxit, verbose, w.init) 
{
    p <- nrow(X)
    n <- ncol(X)
    W <- w.init
    sW <- La.svd(W)
    W <- sW$u %*% diag(1/sW$d) %*% t(sW$u) %*% W
    W1 <- W
    lim <- rep(1000, maxit)
    it <- 1
    mu3<-rep(0,p);
    mu4<-rep(3,p);
    Ftypeold<-rep(0,p);
    Ftype<-rep(1,p)
    tanhconst<-1;
    scoresaturation<-20;
    
    while (it < maxit & ((max(abs(Ftype-Ftypeold))!=0) | lim[it] > tol)) {
            y <- W %*% X
            yt<-t(y); # column vectors
            mu3<-apply(yt^3,2,FUN=mean);            
            mu4old<-mu4;
            mu4<-apply(yt^4,2,FUN=mean);
            if (verbose) { 
                cat('Estimated 3rd and 4th moment:\n')
                print(cbind(mu3,mu4)) }
            score1<-matrix(0,p,n);
            score2<-matrix(0,p,n);
            for (i in 1:p) {
                if ((mu4[i]<border.base[1]+border.slope[1]*mu3[i]^2) | (mu4[i]>border.base[2]+border.slope[2]*mu3[i]^2)) {
                    # tanh is used. See FastICA (http://www.cis.hut.fi/projects/ica/fastica/)
                    # and related papers for more details.
                    Ftype[i]<-3;
                    score1[i,]<-tanh(tanhconst*y[i,]);
                    score2[i,]<-tanhconst*(1-score1[i,]^2);
                } else {
                    # Pearson moment estimation is used. 
                    # See pearson_momentfit.r for more details.
                    Ftype[i]<-1;
                    pearson.par<-pearson_momentfit(mu3[i],mu4[i]);
                    #if (verbose) { print(pearson.par) }
                    pearson.scores<-pearson_score(y[i,],a=pearson.par$a,b=pearson.par$b,C=scoresaturation);
                    score1[i,]<-pearson.scores$score1;
                    score2[i,]<-pearson.scores$score2;
                }
            }  
            v1<-score1 %*% t(X)/n;
            v2<-diag(apply(score2,1,FUN=mean)) %*% W;      
            W1 <- v1 - v2
            sW1 <- La.svd(W1)
            W1 <- sW1$u %*% diag(1/sW1$d) %*% t(sW1$u) %*% W1
            lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
            W <- W1
            if (verbose) 
                cat("Iteration", it, "tol =", format(lim[it + 
                  1]), "\n")
            it <- it + 1
            Ftypeold<-Ftype;
    }
    fasticapearson<-list(a=W,it=it)
}

pearson_score<-function(x,b,a,C)
{
    # PEARSON_SCORE - Calculates the score function and its derivative
    # for a distribution from the Pearson system.
    #
    # Input:
    #   x     = the realization matrix (sample values)
    #   b     = the distribution parameters vector [b0 b1 b2]
    #   a     = the distribution parameters vector [a0 a1]    
    #   C     = the saturation limit (scalar); 
    #           After saturation -C < phi(x) < C for all x
    # Output:
    #   phi   = the score function matrix corresponding to x  
    #   phi2  = the derivative of the score function matrix corresponding to x  
    if (b[3]!=0 & b[3]!=-1) {
        delta<-b[2]^2/(4*b[1]*b[3]);
        if (a[2]==0) {
            ddelta<-1;
        } else {
            ddelta<-1+(delta-1)/(1+(4*b[3]/a[2]*(1+b[3]/a[2]))*delta);
        }
        #ddelta is used to determine the type of distribution. 
        if (ddelta>=1) { 
            if ((-C*b[2]-a[2])^2+4*C*(-C*b[1]+b[2])*b[3]>=0) { 
                C_1<-(-a[2]-C*b[2]-sqrt((-C*b[2]-a[2])^2+4*C*(-C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                C_2<-(-a[2]-C*b[2]+sqrt((-C*b[2]-a[2])^2+4*C*(-C*b[1]+b[2])*b[3]))/(2*C*b[3]);
            } else {
                C_1<-NA;
                C_2<-NA;
            }
            if ((C*b[2]-a[2])^2-4*C*(C*b[1]+b[2])*b[3]>=0) {
                C_3<-(a[2]-C*b[2]-sqrt((C*b[2]-a[2])^2-4*C*(C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                C_4<-(a[2]-C*b[2]+sqrt((C*b[2]-a[2])^2-4*C*(C*b[1]+b[2])*b[3]))/(2*C*b[3]);
            } else {
                C_3<-NA;
                C_4<-NA;
            }
            if (b[2]<0) { #left bound domain
                minscorelimit<-max(c(C_1,C_2,C_3,C_4));
                if (!is.na(minscorelimit)) {
                    x<-pmax(x,minscorelimit);
                }
            } else { #right bound domain
                maxscorelimit<-min(c(C_1,C_2,C_3,C_4));
                if (!is.na(maxscorelimit)) {
                    x<-pmin(x,maxscorelimit); 
                }
            }
        } else { 
            if (ddelta<0 | (ddelta==0 & b[3]>0)) { #if  (delta<=0 & delta>-1/(4*a[2]*b[3]*(1+a[2]*b[3])))   #beta of first kind
                C_1<-(-a[2]-C*b[2]-sqrt((-C*b[2]-a[2])^2+4*C*(-C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                C_2<-(-a[2]-C*b[2]+sqrt((-C*b[2]-a[2])^2+4*C*(-C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                C_3<-(a[2]-C*b[2]-sqrt((C*b[2]-a[2])^2-4*C*(C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                C_4<-(a[2]-C*b[2]+sqrt((C*b[2]-a[2])^2-4*C*(C*b[1]+b[2])*b[3]))/(2*C*b[3]);
                CC<-sort(c(C_1,C_2,C_3,C_4));
                minscorelimit<-CC[2];
                maxscorelimit<-CC[3];
                x<-pmin(pmax(x,minscorelimit),maxscorelimit);
            }   
        }
    } else {
        if (b[3]==0 & b[2]!=0) {  #special case: Gamma distribution
            C_1<-(b[1]*C+a[1])/(a[2]-C*b[2])
            C_2<-(-b[1]*C+a[1])/(a[2]+C*b[2])
            if (b[2]<0) { #left bound Gamma distribution
                minscorelimit<-max(c(C_1,C_2));
                if (!is.na(minscorelimit)) {
                    x<-pmax(x,minscorelimit);
                } 
            } else { #right bound Gamma distribution
                maxscorelimit<-min(c(C_1,C_2));
                if (is.na(maxscorelimit)) {
                    x<-pmin(x,maxscorelimit);
                }
            }
        }
    }
    #If none of the conditions above is fullfilled, no bounding is needed.
    
    #the actual score function is calculated here:
    denominator<-b[1]+b[2]*x+b[3]*x^2;
    phi<--(a[2]*x-a[1])/denominator;
    phi2<--(b[1]*a[2]+b[2]*a[1]+2*a[1]*b[3]*x-a[2]*b[3]*x^2)/(denominator^2);
    return(list(score1=phi,score2=phi2))
}

pearson_momentfit<-function(mu3,mu4)
{
    # PEARSON_MOMENTFIT - Estimates the parameters of the zero mean and unit
    # variance distribution from the Pearson system using the third and
    # forth central moments (the method of moments).
    #
    # Input
    #   mu3    = sample 3rd moment
    #   mu4    = sample 4th moment
    #  
    # Output
    #   pearsona = the distribution parameters vector [a0 a1]
    #   pearsonb = the distribution parameters vector [b0 b1 b2]
    if (mu4<mu3^2+1) {
        warning('Impossible values for moments: ',mu3,mu4,'\n')
    }
    denominator<-abs(10*mu4-12*mu3^2-18);
    b0<--(4*mu4-3*mu3^2);
    b1<--mu3*(mu4+3);
    b2<--(2*mu4-3*mu3^2-6);
    pearsonb<-c(b0,b1,b2);
    pearsona<-c(b1,denominator);
   return(list(a=pearsona,b=pearsonb))
}


PearsonICAdemo<-function(numsig=4,signal_length=5000)
{
    cat('\nThis is a demonstration of the Pearson-ICA algorithm in work.\n');
    cat('Pearson-ICA performs independent component analysis using nonlinearities from the Pearson system.\n')
    #############################################################################
    cat('\nRayleigh distributed source signals are generated:\n')
    cat('Number of signals is ',numsig,'\n')
    cat('Number of observations (signal length) is ',signal_length,'\n')
    lastplot<-round(signal_length/4000*numsig)*20;
    cat(lastplot,' observations from the beginning of the signals will be plotted in the graphs.\n')
    
    sources<-matrix(sqrt(rnorm(numsig*signal_length)^2+rnorm(numsig*signal_length)^2),signal_length,numsig);  
    
    lag_length<-round(lastplot/numsig);
    
    zero_beg<-numeric(numsig);
    zero_end<-numeric(numsig);
    
    for (signal_no in 1:numsig) {
    zero_beg[signal_no]<-(signal_no-1)*lag_length+1;
    zero_end[signal_no]<-signal_no*lag_length;
    source_mean<-mean(sources[-(zero_beg[signal_no]:zero_end[signal_no]),signal_no]);
    sources[zero_beg[signal_no]:zero_end[signal_no],signal_no]<-source_mean;
    }
    cat('\n A sequence of zeros (signal means) of length ',lag_length,'is substituted to\n')
    cat('every source signal such that the indices do not overlap. This is\n') 
    cat('done in order to visualize the separation result. If the final\n') 
    cat('separation is succesful, the same constant sequences should be again\n')
    cat('visible.\n');
    
    cat('\nPlotting the beginning of the source signals.\n')
    dev.new()
    par(mfrow = c(numsig, 1),mar=c(0,0,1,0));
    
    for (signal_no in 1:numsig)
    {
        plot(sources[1:lastplot,signal_no],type='l',ylab="",axes=FALSE);
        axis(1)
        if (signal_no==1) {title('The beginning of the source signals','The substituted means are drawn with red'); }
        interval<-zero_beg[signal_no]:(zero_end[signal_no]-1);
        lines(interval,sources[interval,signal_no],type='l',col='red');  
    }
    userinput<-readline('\nPress enter to continue...\n');
    
    ###############################################################################
    # Mixing matrix
    ###############################################################################
    
    MIX_MATRIX<-matrix(rnorm(numsig^2),numsig,numsig);
    fixed_matrix<-matrix(c(0.7396,    0.9084,    0.2994,    0.3089,
                0.4898,    0.2980,    0.5771,    0.4108,
                0.1096,    0.7808,    0.8361,    0.4669,
                0.4199,    0.8799,    0.2706,    0.7467),4,4,byrow=TRUE);
    #############################################################################
    # To use completely random matrix, comment out the next two lines.
    #############################################################################
    MIX_MATRIX[1:min(numsig,4),1:min(numsig,4)]<-fixed_matrix[1:min(numsig,4),1:min(numsig,4)];
    
    cat('\nThe source signals are linearily mixed together using the matrix\n');
    print(MIX_MATRIX);
    
    ###############################################################################
    # Mixed signals (observation matrix)
    ###############################################################################
    
    mixedsig<-sources %*% t(MIX_MATRIX);
    #############################################################################
    # Plotting the beginning of the mixed signals.
    #############################################################################
    dev.new()
    par(mfrow = c(numsig, 1),mar=c(0,0,1,0));
    for (signal_no in 1:numsig)
    {
        plot(mixedsig[1:lastplot,signal_no],type='l',ylab="",axes=FALSE);
        axis(1);
        if (signal_no==1) { title('The beginning of the observed mixtures'); }
    }
     
    ###############################################################################
    # Separation
    ###############################################################################
    
    cat('\nNow we will try to recover original signals from the\n')
    cat('mixture only using the Pearson-ICA algorithm.\n');
    userinput<-readline("Press enter to begin...");
    
    ############################################################################ 
    # the standard parameter values are used. 
    #############################################################################
    pearson_ica<-PearsonICA(mixedsig,n.comp=numsig,verbose=TRUE) 
    icasig<-pearson_ica$S;
    A<-pearson_ica$A;
    W<-pearson_ica$W;
    print(A)
    
    #############################################################################
    # Plotting the beginning of the separated signals.
    #############################################################################
    dev.new()
    par(mfrow = c(numsig, 1),mar=c(0,0,1,0));
    for (signal_no in 1:numsig)
    {
        plot(icasig[1:lastplot,signal_no],type='l',ylab="",col='blue',axes=FALSE);
        axis(1);
        if (signal_no==1) { title('The beginning of the separated signals','Compare the result to the originals source signals') }
    }
    
    userinput<-readline('\nPress enter to continue...\n');
    
    cat('\nIn order to compare directly the results, the sign and the permutation \n')
    cat('indeterminency must be resolved. This is done by selecting the largest \n')
    cat('absolute row values (and the corresponding signs) from \n')
    cat('the "result matrix" W*A. \n')
    Y<-t(MIX_MATRIX)%*%W;
    permutation<-max.col(abs(Y));
    source_sign<-sign(diag(Y[,permutation]));
    
    cat('\n The signals must have the same scale. This is done
    by normalizing the sources and the separation results. \n')
    normalized_sources<-scale(sources) 
    normalized_icasig<-scale(icasig)
    
    normalized_icasig<-normalized_icasig[,permutation];
    normalized_icasig<-t(source_sign*t(normalized_icasig));
    
    cat('\nPlotting the beginning of the normalized signals.\n') 
    dev.new()
    par(mfrow = c(numsig, 1),mar=c(0,0,1,0));
    for (signal_no in 1:numsig)
    {
        plot(normalized_sources[1:lastplot,signal_no],type='l',ylab="",col='red',axes=FALSE);
        axis(1);
        lines(normalized_icasig[1:lastplot,signal_no],type='l',col='blue')
        if (signal_no==1) { title('The beginning of the normalized source (red) and the separated signals (blue)') }
    }
}
   
