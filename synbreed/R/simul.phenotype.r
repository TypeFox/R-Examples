simul.phenotype <- function(pedigree=NULL,A=NULL,mu=100,vc=NULL,Nloc=1,Nrepl=1){

        if(is.null(A) & is.null(pedigree)) step("either 'pedigree' or 'A' must be given")
        # create A matrix if missing
        if (is.null(A)) A <- kin(pedigree,ret="add")
        if (is.null(vc)) stop("missing variance components")
        
        # read data out of arguments
        N <- nrow(A)         
        n <- N*Nloc*Nrepl
        if(!is.null(vc$sigma2e)) sigmae <- sqrt(vc$sigma2e)
        else sigmae <- 1
        if(!is.null(vc$sigma2a)) sigmaa <- sqrt(vc$sigma2a)
        else sigmaa <- 0
        if(!is.null(vc$sigma2l)) sigmal <- sqrt(vc$sigma2l)
        else sigmal <- 0
        if(!is.null(vc$sigma2b)) sigmab <- sqrt(vc$sigma2b)
        else sigmab <- 0
        
        namesA <- rownames(A) 
        if(is.null(rownames(A))) namesA <- 1:N  
                                     
        # initialize data

        ID <- rep(namesA,each=Nrepl*Nloc)
        Loc <- rep(1:Nloc,length.out=n,each=Nrepl)
        Block <- rep(1:Nrepl,lengt.out=n)
        
         # as matrix
        A <- matrix(A,nrow=N)  
                       

        # simulate data for contribution parts of phenotypic values
        
        # true breeding values
        #tbv <- rmvnorm(1,rep(0,N),A*sigmaa^2)
        tbv <- sigmaa*(chol(A)%*%rnorm(N,0,1))
        tbv <- rep(tbv,each=Nloc*Nrepl)
        # location effect
        locEff <- rnorm(Nloc,0,sigmal)
        locEff <- rep(locEff,length.out=n,each=Nrepl)
        # block effects 
        blockEff <- rnorm(Nrepl,0,sigmab)
        blockEff <- rep(blockEff,length.out=n)
        # residual
        residual <- rnorm(n,0,sigmae)
        
        # simulate phenotypic value
        Trait <- mu + tbv + locEff + blockEff + residual

        # combine to a data.frame
        ret <- data.frame(ID=factor(ID),Loc=factor(Loc),Block=factor(Block),Trait=Trait,TBV=tbv)
        return(ret)

}

