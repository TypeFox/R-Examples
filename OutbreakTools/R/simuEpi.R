
if(getRversion() >= "2.15.1")  utils::globalVariables(c("date","number","category"))

## Simulate an epidemic following a SIR model
## N Size of the population
## D Duration of simulation
## beta Rate of infection
## nu Rate of recovery
## L Length of genetic sequences
## mu Probability of mutation per base per transmission event
## Returns to simulated epidemic as an obkData object
## Xavier Didelot
simuEpi <- function (N=1000,D=10,beta=0.001,nu=0.1,L=1000,mu=0.001, plot=TRUE, makePhyloTree=FALSE) {
    S<-matrix(0,D,3)
    T<-matrix(0,N,3)
    dates<-matrix("",N,1)
    seqs<-matrix('c',N,L)
    lets<-c('a','c','g','t')
    S[1,1]=N-1;S[1,2]=1;S[1,3]=0
    T[1,1]=1;T[1,2]=NA;T[1,3]=0
    seqs[1,]=sample(lets,L,replace=TRUE)
    dates[1,1]=as.character(as.Date(0,origin="2000-01-01"));
    curinf=1;
    ninf=1;
    for (i in 2:D) {
        inf=rbinom(1,S[i-1,1],1-(1-beta)^S[i-1,2])
        rec=rbinom(1,S[i-1,2],nu)
        if (inf>0) for (j in 1:inf) {
            T[ninf+1,1]=ninf+1;
            T[ninf+1,2]=sample(curinf,1);
            T[ninf+1,3]=i;
            seqs[ninf+1,]=seqs[T[ninf+1,2],]
            muts=runif(L)<mu
            muts=seq(along=muts)[muts!=FALSE]
            for (j in muts) seqs[ninf+1,j]=sample(setdiff(lets,seqs[ninf+1,j]),1)
            dates[ninf+1,1]=as.character(as.Date(i-1,origin="2000-01-01"));
            ninf=ninf+1;
        }
        curinf=.resample(curinf, length(curinf)-rec)
        if (inf>0) curinf=c(curinf,(ninf-inf+1):ninf)
        S[i,1]=S[i-1,1]-inf
        S[i,2]=S[i-1,2]+inf-rec
        S[i,3]=S[i-1,3]+rec
    }
    T=T[1:ninf,]
    dates=.process.Date(dates[1:ninf,])
    seqs=seqs[1:ninf,,drop=FALSE]
    rownames(seqs)<-1:ninf

    ## make contact matrix from -> to ##
    contacts <- na.omit(data.frame(from=T[,2], to=1:ninf))


    if (makePhyloTree){
        ## make a simple phylogeny from the sequences
        simplephylo=nj(dist.dna(as.DNAbin(seqs))) # quick and dirty but doesn't assume homochronous samples
        simplephylo=root(simplephylo,1)

        ## ## make 3 plots in this case: the epidemic, the trans tree and the phylo tree
        ## if (showPlots) {par(mfrow=c(1,3)); plotEpi(S); plot(infectorTableToNetwork(T),main="Transmission tree"); plot(simplephylo,main="NJ phylogeny")}

        ## convert to multiphylo
        simplephylo=list(simplephylo)
        class(simplephylo)="multiPhylo"  # this class is required by the obkData structure

        ## include the NJ tree in the obkData object that will be returned
        ret <- new("obkData",individuals=data.frame("individualID"=1:ninf,"infector"=T[,2],"DateInfected"=dates),
                   dna=as.DNAbin(seqs), dna.date=dates, dna.individualID=1:ninf, trees=simplephylo,
                   contacts=contacts, contacts.start=dates, contacts.end=dates, contacts.directed=TRUE)
    } else {     ## otherwise, just do 2 plots and don't make any phylogeny
        ## if (showPlots) {
        ##     par(mfrow=c(1,2)); plotEpi(S); plot(infectorTableToNetwork(T),main="Transmission tree")
        ## }
	ret <- new("obkData",individuals=data.frame("individualID"=1:ninf,"infector"=T[,2],"DateInfected"=dates),
                   dna=as.DNAbin(seqs), dna.date=dates, dna.individualID=1:ninf,
                   contacts=contacts, contacts.start=dates, contacts.end=dates, contacts.directed=TRUE)
    }


    ## MAKE POPULATION DYNAMICS TABLE ##
    S <- as.data.frame(S)
    names(S) <- c("Susceptible","Infected","Recovered")

    ## add dates
    S$date <- seq(from=min(dates), by=1, length=nrow(S))

    out <- list(x=ret, dynamics=S)

    ## MAKE PLOTS IF REQUESTED ##
    if(plot){
        ## build long form data.frame for ggplot
        df <- melt(S, id.var="date")
        names(df) <- c("date","category","number")

        ## make plot
        p <- ggplot(df, aes(x=date, y=number)) + geom_point(aes(colour=category),shape=20,size=4)
        p <- p + geom_line(aes(colour=category), size=3, alpha=I(.5))
        print(p)
        out$plot <- p
    }
    return(out)
} # end simuEpi




.resample <- function(x, ...) x[sample.int(length(x), ...)]
                                        # because sample can take an integer argument and then samples from 1:n,
                                        # need this resampler to make simuEpi work in cases where length(curinf) is 1 and the value of curinf is an integer.



                                        #Plot the number of susceptible, infected and recovered as a function of time
                                        #S Matrix containing the numbers to be plotted
                                        #Xavier Didelot
plotEpi <- function(S) {
    plot(c(0,dim(S)[1]),c(0,sum(S[1,])),type='n',xlab='Days',ylab='Individuals',main="Epidemic summary")
    lines(S[,1],col='black')
    lines(S[,2],col='red')
    lines(S[,3],col='blue')
    legend('right',lty=c(1,1),col=c('black','red','blue'),c('Susceptible','Infected','Recovered'))
}

##                                     #Convert transmission tree to a network
##                                         #transmissiontreeData Matrix of who infected whom
##                                         #Network of who infected whom
##                                         #Caroline Colijn
## infectorTableToNetwork <- function (transmissiontreeData)
## {
##     uniqueIDs <- sort(c(unique(as.character(transmissiontreeData[,1]),as.character(transmissiontreeData[,2]))))
##     nUniqueIDs <- length(uniqueIDs)
##     edgeList <- na.omit(transmissiontreeData[,1:2])
##     numEdges <- dim(edgeList)[1]
##     y <- network.initialize(nUniqueIDs)
##     network.vertex.names(y) <- uniqueIDs
##     for(i in 1:numEdges){
##         v1 <- match(as.character(edgeList[i,1]),uniqueIDs)
##         v2 <- match(as.character(edgeList[i,2]),uniqueIDs)

##         if (!is.na(v2)) add.edges(y,v2,v1)
##     }
##     return(y)
## }







## Create phylogenetic tree from transmission tree only, just in case there is no sequence data
## This is not relevant; now we simulate sequence data and the phylo tree returned by simuEpi is
## derived from these.
## param transmissiontreeData Matrix of who infected whom
## return phylogenetic tree representing how samples of the infectious agents may be related
## author Caroline Colijn
phylofromtranstree <- function(transmissiontreeData){

    ## use the transmission tree data to create an *undirected* network
    uniqueIDs <- sort(c(unique(as.character(transmissiontreeData[,1]),as.character(transmissiontreeData[,2]))))
    nUniqueIDs <- length(uniqueIDs)
    edgeList <- na.omit(transmissiontreeData[,1:2])

    numEdges <- dim(edgeList)[1]
    y <- network.initialize(nUniqueIDs)
    network.vertex.names(y) <- uniqueIDs
    for(i in 1:numEdges){
        v1 <- match(as.character(edgeList[i,1]),uniqueIDs)
        v2 <- match(as.character(edgeList[i,2]),uniqueIDs)

        if (!is.na(v2) & !is.na(v1)) add.edges(y,v2,v1)
        if (!is.na(v1) & !is.na(v2)) add.edges(y,v1,v2)
    }
    mynet <- y

    ## get pairwise shortest path distances in this network and use them to make a phylogeny in one of 2 very simple ways
    phylotree <- nj(geodist(mynet)$gdist)

    return(phylotree)
}

