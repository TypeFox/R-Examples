RobPer <- function(ts, weighting, periods, regression, model, steps=10, tol=1e-3, var1=weighting,
    genoudcontrol=list(pop.size=50, max.generations=50,  wait.generations=5), LTSopt=TRUE,
    taucontrol=list(N=100, kk=2, tt=5, rr=2, approximate=FALSE), Scontrol=list(N=ifelse(weighting,200,50),
    kk=2, tt=5, b=.5, cc=1.547, seed=NULL)) {
    # Function to calculate periodograms based on robust M regression
    # Check arguments
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(abs(x) - round(x)) < tol
    # ts
    if(dim(ts)[2]!=3& weighting)     stop("if weighting==TRUE is chosen, ts needs to be a matrix or data frame with three columns.") 
    if(dim(ts)[2]<2 & (!weighting))    stop("ts has not enough columns")
    if(dim(ts)[2]>3 & (!weighting))    stop("ts has too many columns")
    if(!is.numeric(as.matrix(ts)))   stop("ts needs to be numeric")
    # weighting
    if(length(weighting)!=1)     stop("weighting is too long or not given")        
    if(!is.logical(weighting))   stop("weighting is not logical") 
    if(weighting){ if(any(ts[,3]<=0, na.rm=TRUE)) stop("When choosing weighting=TRUE, measurement accuracies have to be strictly positive.")}
    # periods
    if(length(periods)<1)      stop("please spezify trial periods using the argument periods") 
    if(!is.numeric(periods))   stop("periods needs to be numeric") 
    if(any(periods<=0, na.rm=TRUE))        stop("All trial periods need to be greater than 0") 
    if(any(is.na(periods)))    stop("NA periods are not allowed") 
    # regression    
    if(length(regression)!=1)                      stop("Please spezify one single method of regression") 
    if(any(!regression%in%c("huber", "bisquare", "LTS", "L2", "L1", "S", "tau"))) stop("regression method not implemented") 
    # model    
    if(any(!model%in%c("step", "2step", "sine", "fourier(2)", "fourier(3)", "splines"))) stop("model not implemented") 
    # steps
    if(model %in% c("step", "2step")){
        if(is.na(steps))           stop("Please specify steps")
        if(length(steps)!=1)       stop("steps has to be a single number")
        if(!is.wholenumber(steps)) stop("steps needs to be a positive integer number")
    }
    # tol
    ltslog <- regression=="LTS"
    if(ltslog) ltslog<- LTSopt==TRUE
    if(regression%in%c("huber", "bisquare")|ltslog){
        if(is.na(tol))          stop("Please specify tol") 
        if(length(tol)!=1)      stop("tol needs to be one positive number")        
        if(!is.numeric(tol))    stop("tol is an accuracy for convergence criteria and needs to be a positive number")
        if(tol<=0)              stop("tol has to be strictly positve")
    }
    # var1
    if(regression%in%c("huber", "bisquare")){
        if(is.na(var1))        stop("please specify var1")
        if(length(var1)!=1)    stop("var1 is too long or not given")        
        if(!is.logical( var1)) stop("var1 needs to be logical") 
    }
    # genoudcontrol
    ltslog<- regression=="LTS"
    if(ltslog) ltslog<- LTSopt==TRUE
    if(regression=="bisquare"|ltslog){
        # fill genoudcontrol with defaults in case list is not (completely) specified
        if(is.null(genoudcontrol$max.generations)) genoudcontrol$max.generations <- 100
        if(is.null(genoudcontrol$pop.size)) genoudcontrol$pop.size <- 100
        if(is.null(genoudcontrol$wait.generations)) genoudcontrol$wait.generations <- 5
        # reduce genoudcontrol on sensible entries
        genoudcontrol <- genoudcontrol[which(names(genoudcontrol)%in% c("max.generations", "wait.generations", "pop.size"))]
        # check arguments
        if(any(is.na(genoudcontrol))) stop("NA in genoudcontrol is not allowed")
        if(any(sapply(genoudcontrol, length)>1)) stop("All list entries for genoudcontrol have to be single numbers")
        if(any(!sapply(genoudcontrol, is.wholenumber))) stop("All list entries for genoudcontrol have to be positive integers")
        if(any(genoudcontrol<=0)) stop("All list entries for genoudcontrol have to be strictly positive")
    }
    # LTSopt
    if(regression=="LTS") {
        if(length(LTSopt)!=1)   stop("LTSopt is too long or not given")
        if(!is.logical(LTSopt)) stop("LTSopt is not logical")
    }
    # taucontrol
    if(regression =="tau") {
        # check list and fill it up
        if(is.null(taucontrol$N))     taucontrol$N<-500
        if(is.null(taucontrol$kk))    taucontrol$kk<-2
        if(is.null(taucontrol$tt))    taucontrol$tt<-5
        if(is.null(taucontrol$rr))    taucontrol$rr<-2
        if(is.null(taucontrol$approximate)) taucontrol$approximate<- FALSE
        # reduce to sensible entries
        taucontrol <- taucontrol[which(names(taucontrol)%in% c("N", "kk", "tt", "rr", "approximation"))]
        # check arguments
        if(any(is.na(taucontrol))) stop("NA in taucontrol is not allowed")
        if(any(sapply(taucontrol, length)>1)) stop("All list entries for taucontrol have to be onedimensional")
        if(any(!sapply(taucontrol[which(names(taucontrol)%in% c("N", "kk", "tt", "rr"))], is.wholenumber))) {
            stop("List entries N, kk, tt and rr in taucontrol have to be positive integers")
        }
        if(any(taucontrol[which(names(taucontrol)%in% c("N", "kk", "tt", "rr"))]<=0)) {
            stop("List entries N, kk, tt and rr in taucontrol have to be strictly positive")
        }
        if(!is.null(taucontrol$approximate)) {
            if(!is.logical(taucontrol$approximate)) stop("taucontrol$approximate is not logical")
        }
    }
    
    # Scontrol
    if(regression =="S") {
        # check list and fill it up
        Scontrol$int <- FALSE
        if(is.null(Scontrol$N))  {if(weighting) Scontrol$N <- 200 else Scontrol$N <- 50}
        if(is.null(Scontrol$kk))     Scontrol$kk <- 2
        if(is.null(Scontrol$tt))     Scontrol$tt <- 5
        if(is.null(Scontrol$b))      Scontrol$b  <-  .5
        if(is.null(Scontrol$cc))     Scontrol$cc <- 1.547
        # reduce to sensible entries
        Scontrol<-Scontrol[which(names(Scontrol)%in% c("int", "N", "kk", "tt", "b","cc", "seed"))]
        # check arguments
        if(any(is.na(Scontrol))) stop("NA in Scontrol is not allowed")
        if(any(sapply(Scontrol, length)>1)) stop("All list entries for taucontrol have to be onedimensional")
        if(any(!sapply(Scontrol[which(names(Scontrol)%in% c("N", "kk", "tt", "b","cc"))], is.numeric))) {
            stop("List entries N, kk, tt, b and cc in Scontrol have to be numeric")
        }
        if(any(!sapply(Scontrol[which(names(Scontrol)%in% c("N", "kk", "tt"))], is.wholenumber))) {
            stop("List entries N, kk and tt in Scontrol have to be positive integers")
        }
        if(any(Scontrol[which(names(Scontrol)%in% c("N", "kk", "tt", "b","cc"))]<=0)) {
            stop("List entries N, kk,tt, b and cc in Scontrol have to be strictly positive")
        }
    }        

            
    # Clean time series: Delete incomplete observations
    if(weighting) lookatcol <- 3 else lookatcol <- 2
    if(any(is.na(ts[,1:lookatcol]))) message("There exist incomplete observations which are ignored.")
    if(weighting)  ts<- ts[complete.cases(ts),] else ts<- ts[complete.cases(ts[,1:2]),]
	
    # Define some variables
    n <- dim(ts)[1]
    tt<- ts[,1]
    y <- ts[,2]
    periodogram<- numeric(length(periods))
	
    if(missing(steps)) steps<-10
	
    # Weight observation (or not)
	
    if(weighting) {
        s <- ts[,3]
        if(any(s<0)) stop("Negative measurement accuracies are not allowed, please change ts.")
        # it is assumed that s is not estimated for observations with s=0.
        # Therefor entrys of s equal to zero are replaced by the mean of the non-zero s-entries
        if(any(s==0)) {
            s[s==0]<-mean(s[s!=0])
            message("Measurement accuracies equal to zero are replaced by the mean of the nonzero measurement accuracies.")
        }
    } else s <- rep(1,n)
            
    if(isTRUE(all.equal(var(y),0))) stop("observations are constant")
            
    # Choose rho and W or zeta according to regression method
    if(regression=="L2"){ zeta <- function(r, q=1) sum(r^2)}        # q is dummy
    if(regression=="L1"){ zeta <- function(r, q=1) sum(abs(r))}     # q is dummy
    if(regression=="LTS"){ zeta <- function(r,q=1) sum(sort(r^2)[1:(floor(n/2)+floor((q+1)/2))]) } # default q=1 is arbitrary
    if(regression=="huber") { 
        # rho according to Maronna, Martin and Yohai (2006) , pages 26-30 
        # k chosen in order to have 95% efficiency 
        rho <- function(x) {
            k<-1.345
            ifelse(abs(x)<=k, x^2, 2*k*abs(x)-k^2)
        }
        W <- function(x) {
            k<- 1.345
            ifelse(abs(x)<=k, 1, k/abs(x))
        }
    }
    if(regression=="bisquare") {
        # rho according to Maronna, Martin and Yohai (2006) , pages 26-30
        # k chosen in order to have 95% efficiency
        rho <- function(x) {
            k <- 4.68
            ifelse(abs(x)<=k, 1-(1-(x/k)^2)^3, 1)
        }
        W <- function(x){
            k<-4.68
            ifelse(abs(x)<=k, (1-(x/k)^2)^2, 0)
        }
    }
    if(regression=="S" | regression=="tau"){ zeta <- numeric()}
                        
    yy<- y/s 
    ii<-cbind(rep(1,n)/s )
        
    # In case of "2step" first calculate a step model
    if(model=="2step") design <- "step" else design <- model
        
    # if possible, compute intercept estimate in advance                  
    gamma <- NA
    K_gamma<-NA
    if(regression=="LTS" & !(design %in%c("step", "stepB"))) {
        if(design=="sine") q<-3
        if(design=="fourier(2)") q<-5
        if(design=="fourier(3)") q<-7
        if(design=="splines") q<-4
        alpha <- (n+q-3)/(2*(n-2)) # use trimming h of full model
        if(length(unique(ii))==1) {
            tempINT  <- ltsReg(yy~1  , use.correction=FALSE, alpha=alpha)
        } else {
            tempINT  <- ltsReg(yy~0+ii,use.correction=FALSE, alpha=alpha)
        }
        ee <- tempINT$residuals
        K_gamma<- zeta(ee, q=q)
        gamma<- tempINT$coeff
    }
	
    if(regression=="L2") {
        ee <- lm(yy~0+ii)$residuals
        K_gamma<- zeta(ee)
    }
	
    if(regression=="L1") {
        ee <- suppressWarnings(rq(yy~0+ii, tau=.5, method="br")$residuals)
        K_gamma <- zeta(ee)
    }
	
    if(regression=="S") {
        temp_gamma<- FastS(ii,yy,Scontrol=Scontrol)
        K_gamma<-temp_gamma$scale
        gamma<- temp_gamma$beta
    }
	
    if(regression=="tau") {
        temp_gamma <- FastTau(ii,yy, taucontrol=taucontrol)
        K_gamma <- temp_gamma$scale
        gamma<- temp_gamma$beta
    }
	
    # Which function use to calculate the periodogram?
    if(regression %in% c("huber", "bisquare")) {
        singleFUN <- function(p_ind, design) {
            singlePerM(tt=tt,yy=yy,s=s,n=n,ii=ii, pp=periods[p_ind],W=W,rho=rho,regression=regression,
                design=design, steps=steps, var1=var1, tol=tol, genoudcontrol=genoudcontrol)
        }
    } else {
        singleFUN <- function(p_ind, design) {
            singlePernotM(tt=tt,yy=yy,s=s,n=n,ii=ii,pp=periods[p_ind], design=design, regression=regression, steps=steps,
                zeta=zeta, LTSopt=LTSopt, tol=tol, genoudcontrol=genoudcontrol, gamma=gamma, K_gamma=K_gamma,
                taucontrol=taucontrol, Scontrol=Scontrol)
        }
    }
	
    # calculate periodogram
	
    for(p_ind in seq(along=periods)) {
        periodogram[p_ind]<- singleFUN(p_ind, design)
    }
	
    # calculate second step model if model=="2step"
	
    if(model=="2step") {
        periodogram2<- numeric(length(periods))
        for(p_ind in seq(along=periods)) {
            periodogram2[p_ind]<- singleFUN(p_ind, design="stepB")
        }
        periodogram<-apply(cbind(periodogram, periodogram2),1, mean)
    }
	
    return(periodogram)
}
