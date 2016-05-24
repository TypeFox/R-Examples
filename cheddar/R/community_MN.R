# Community functions that operate on M and/or N
Log10MNBiomass <- function(community)
{
    # A matrix containing columns 'Log10M', 'Log10N' and 'Log10Biomass'
    np <- NPS(community, c('M','N'))
    res <- cbind(Log10M=log10(np$M), Log10N=log10(np$N), 
                 Log10Biomass=log10(np$M * np$N))
    rownames(res) <- np$node
    return (res)
}

Log10M <- function(community)
{
    # Returns a vector of log10-transformed. All elements are NA if 'M' is 
    # not a node property. Community() does not allow M==0. 
    return (log10(NP(community, 'M')))
}

Log10N <- function(community)
{
    # Returns a vector of log10-transformed. All elements are NA if 'N' is 
    # not a node property. Community() does not allow N==0. 
    return (log10(NP(community, 'N')))
}

Biomass <- function(community)
{
    # Returns a vector of biomasses. All elements are NA if 'M' or 'N' are 
    # not node properties. Elements are NA for nodes that have N and/or M of 
    # NA or N of 0.
    return (NP(community, 'M') * NP(community, 'N'))
}

Log10Biomass <- function(community)
{
    # Returns a vector of log10-transformmed biomasses. All elements are NA 
    # if 'M' or 'N' are not node properties.
    return (log10(Biomass(community)))
}

RCMRatio <- function(community)
{
    # Returns a vector containing the M of the resource divided by the M of the 
    # consumer for every consumer. All elements are NA if 'M' is not a node 
    # property.
    tlps <- TLPS(community, node.properties='M')
    return (tlps$resource.M / tlps$consumer.M)
}

Log10RCMRatio <- function(community)
{
    # log10(RCMRatio)
    return (log10(RCMRatio(community)))
}

CRMRatio <- function(community)
{
    # Returns a vector containin the M of the resource divided by the M of the 
    # consumer for every consumer. All elements are NA if 'M' is not a node 
    # property.
    tlps <- TLPS(community, node.properties='M')
    return (tlps$consumer.M / tlps$resource.M)
}

Log10CRMRatio <- function(community)
{
    # log10(CRMRatio)
    return (log10(CRMRatio(community)))
}

.NodesWithoutMOrN <- function(community)
{
    # Returns the names of nodes that have a M or N of NA
    return (names(which(is.na(NP(community, 'M')) | is.na(NP(community, 'N')))))
}

BodyMassBins <- function(community, 
                         lower=min(NP(community,'M'), na.rm=TRUE), 
                         upper=max(NP(community,'M'), na.rm=TRUE), 
                         n.bins=10)
{
    if(!is.Community(community)) stop('Not a Community')

    stopifnot(n.bins>2)
    stopifnot(upper>lower)

    np <- NPS(community)

    stopifnot(all(is.na(np$M) | (np$M>=lower & np$M<=upper)))

    # Bins
    breaks <- seq(log10(lower), log10(upper), length.out=1+n.bins)
    bin.centres <- head(breaks, -1) + diff(breaks)[1]/2
    binned <- cut(log10(np$M), breaks=breaks, labels=FALSE, include.lowest=TRUE)
    names(binned) <- np$node
    attr(binned, 'bin.centres') <- bin.centres
    attr(binned, 'breaks') <- breaks
    return (binned)
}

ResourceLargerThanConsumer <- function(community)
{
    # A matrix of consumer-resource pairs in which the resource has a larger 
    # body mass than the consumer. Columns are the same as 
    # those returned by LinkProperties().
    # Returns NULL if community does not have either M or trophic links
    if(!is.Community(community)) stop('Not a Community')

    if('M' %in% colnames(NPS(community)) && !is.null(TLPS(community)))
    {
        tlp <- TLPS(community, 'M')
        return (tlp[!is.na(tlp[,'resource.M']) & !is.na(tlp[,'consumer.M']) & 
                tlp[,'resource.M']>tlp[,'consumer.M'],,drop=FALSE])
    }
    else
    {
        return (NULL)
    }
}

SumMByClass <- function(community, class, na.rm=FALSE)
{
    return (ApplyByClass(community, 'M', class, sum, na.rm=na.rm))
}

SumNByClass <- function(community, class, na.rm=FALSE)
{
    return (ApplyByClass(community, 'N', class, sum, na.rm=na.rm))
}

SumBiomassByClass <- function(community, class, na.rm=FALSE)
{
    return (ApplyByClass(community, 'Biomass', class, sum, na.rm=na.rm))
}

NvMLinearRegressions <- function(community, class)
{
    # Returns a list of linear regressions through all data and per category. 
    # Nodes with M=NA, N=NA or N=0 are ignored

    # Similar check in LinearRegressionByClass() but the call to missing() in 
    # LinearRegressionByClass() won't work
    if(!is.Community(community)) stop('Not a Community')

    if(missing(class))
    {
        if('category' %in% NodePropertyNames(community))
        {
            class <- 'category'
        }
        else 
        {
            class <- NULL
        }
    }

    LinearRegressionByClass(community, 'Log10M', 'Log10N', class)
}

NvMSlopeAndIntercept <- function(community)
{
    # The slope and intercept of a line through log10(N) versus log10(M) data
    models <- NvMLinearRegressions(community, class = NULL)
    return(c(slope=unname(coef(models[["all"]])[2]), 
             intercept=unname(coef(models[["all"]])[1])))
}

NvMSlope <- function(community)
{
    # The slope of a line through log10(N) versus log10(M) data
    return (unname(NvMSlopeAndIntercept(community)['slope']))
}

NvMIntercept <- function(community)
{
    # The intercept of a line through log10(N) versus log10(M) data
    return (unname(NvMSlopeAndIntercept(community)['intercept']))
}

NvMSlopeAndInterceptByClass <- function(community, class)
{
    # The slopes and intercepts of lines through log10(N) versus log10(M) data
    models <- NvMLinearRegressions(community, class)
    slopes <- sapply(models, function(m) ifelse(is.null(m), NA, coef(m)[2]))
    names(slopes) <- paste('slope.', names(slopes), sep='')
    intercepts <- sapply(models, function(m) ifelse(is.null(m), NA, coef(m)[1]))
    names(intercepts) <- paste('intercept.', names(intercepts), sep='')
    return (c(slopes, intercepts))
}

NvMSlopeByClass <- function(community, class)
{
    # The slopes of lines through log10(N) versus log10(M) data
    res <- NvMSlopeAndInterceptByClass(community, class)
    return (res[grep('slope*', names(res))])
}

NvMInterceptByClass <- function(community, class)
{
    # The intercepts of lines through log10(N) versus log10(M) data
    res <- NvMSlopeAndInterceptByClass(community, class)
    return (res[grep('intercept*', names(res))])
}

NvMTriTrophicStatistics <- function(community)
{
    # Exclude cannibals and all nodes with missing M and/or N
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    .RequireTrophicLinks(community)

    community <- RemoveNodes(community, remove=with(NPS(community), 
                                                    node[is.na(M) | is.na(N)]))
    community <- RemoveCannibalisticLinks(community)

    lp <- TLPS(community, link.properties='.NvMTrophicLinkProperties')
    tnc <- ThreeNodeChains(community, 
                           chain.properties='.NvMThreeNodeChainProperties')
    tc <- TrophicChains(community, 
                        chain.properties='.NvMTrophicChainProperties')
    return (list(links=lp, three.node.chains=tnc, trophic.chains=tc))
}

.NvMTrophicLinkProperties <- function(community)
{
    # Returns a data.frame containing columns resource and consumer, one row 
    # for each link in the food web for which resource and consumer have both 
    # M and N.
    # The returned data.frame contains columns length, angle and 
    # slope, as defined by Cohen et al 2010 PNAS, p 22335-22336.

    # Cohen et al 2010 PNAS, p 22335 and 22336
    # The angle of a link (or link angle) was the counter-
    # clockwise angle to the link from a horizontal arrow starting from R
    # and pointing right parallel to the positive log(M)-axis, and took
    # values in the interval [-180, 180) (Fig. 2A)]. (The angle is not
    # defined when MR = MC and NR = NC, as in cannibalism, for
    # example.) If the link angle equaled -45, then the link had slope -1
    # because tan(-45) = tan(-pi/4 radians) = -1

    chains <- TLPS(community)
    log10M <- Log10M(community)
    log10N <- Log10N(community)
    x <- log10M[chains[,'consumer']] - log10M[chains[,'resource']]
    y <- log10N[chains[,'consumer']] - log10N[chains[,'resource']]

    # Link lengths
    length <- abs(x) + abs(y)

    # Link angles
    # Cohen et al 2010 PNAS, Fig 2 A
    UnsafeArg <- function(n)
    {
        # Like Arg() but returns NA for complex numbers with 0==real and 
        # 0==imag 
        r <- Arg(n)
        r[0==Re(n) & 0==Im(n)] <- NA
        return (r)
    }

    angle <- UnsafeArg(complex(real=x, imaginary=y)) * 360 / (2*pi)
    angle[180==angle] <- -180
    stopifnot(all( (angle>=-180 & angle<180) | is.na(angle)))

    res <- cbind(length=length, angle=angle, slope=y/x)
    rownames(res) <- NULL
    return (res)
}

.NvMThreeNodeChainProperties <- function(community, chains)
{
    # Returns a matrix containing columns bottom, intermediate, top, one row 
    # for each three-species chain in the food web. Also contains the 
    # columns Llower, Lupper, two.span, Alower, Aupper and Abetween as 
    # defined by Cohen et al 2010 PNAS, Fig 2 B/C and p 22336-22337
    # Cohen et al 2010 PNAS, Fig 2 B/C and p 22336 - 22337

    # chains should be a data.frame as computed by ThreeNodeChains
 
    # B -> I -> T

    log10M <- Log10M(community)
    log10N <- Log10N(community)

    # Lengths
    Llower <- abs(log10M[chains[,1]] - log10M[chains[,2]]) +
              abs(log10N[chains[,1]] - log10N[chains[,2]])
    Lupper <- abs(log10M[chains[,2]] - log10M[chains[,3]]) +
              abs(log10N[chains[,2]] - log10N[chains[,3]])

    # 2 span - distance from R to C
    two.span <- abs(log10N[chains[,3]] - log10N[chains[,1]]) +
                abs(log10M[chains[,3]] - log10M[chains[,1]])

    # Angles
    lower <- complex(     real=log10M[chains[,2]] - log10M[chains[,1]], 
                     imaginary=log10N[chains[,2]] - log10N[chains[,1]])
    upper <- complex(     real=log10M[chains[,3]] - log10M[chains[,2]], 
                     imaginary=log10N[chains[,3]] - log10N[chains[,2]])

    UnsafeArg <- function(n)
    {
        # Like Arg() but returns NA for complex numbers with 0==real and 
        # 0==imag 
        r <- Arg(n)
        r[0==Re(n) & 0==Im(n)] <- NA
        return (r)
    }

    Alower <- UnsafeArg(lower) * 360/(2*pi)
    Alower[180==Alower] <- -180
    Aupper <- UnsafeArg(upper) * 360/(2*pi)
    Aupper[180==Aupper] <- -180
    Abetween <- UnsafeArg(upper * Conj(lower)) * 360/(2*pi)
    Abetween[180==Abetween] <- -180

    stopifnot(all( (Abetween>=-180 & Abetween<180) | is.na(Abetween)))
    x <- which(chains[,1]==chains[,3])

    stopifnot(all.equal(Abetween[x], rep(-180, length(x))))

    res <- cbind(Llower, Lupper, two.span, Alower, Aupper, Abetween)
    rownames(res) <- NULL
    return (res)
}

.NvMTrophicChainProperties <- function(community, chains)
{
    # Properties by Cohen et al 2010 PNAS, p 22337, of all the trophic chains 
    # through the food web

    # chains should be a data.frame as computed by TrophicChains

    log10M <- Log10M(community)
    log10N <- Log10N(community)

    # Name of the last node in each chain
    chains <- as.matrix(chains)
    count.chain.length <- apply(chains, 1, function(row) max(which(""!=row)))
    last.in.chain <- sapply(1:length(count.chain.length), 
                            function(l) chains[l, count.chain.length[l]])
    count.chain.length <- count.chain.length-1

    chain.span <- abs(log10N[chains[,1]] - log10N[last.in.chain]) +
                  abs(log10M[chains[,1]] - log10M[last.in.chain])

    # Link lengths
    chains.log10M <- sapply(chains, function(col) log10M[col])
    dim(chains.log10M) <- dim(chains)
    chains.log10N <- sapply(chains, function(col) log10N[col])
    dim(chains.log10N) <- dim(chains)

    # Sum chain lengths
    upper.cols <- 2:ncol(chains)
    lower.cols <- 1:(ncol(chains)-1)
    l <- abs(chains.log10M[,upper.cols,drop=FALSE] - 
             chains.log10M[,lower.cols,drop=FALSE]) +
         abs(chains.log10N[,upper.cols,drop=FALSE] - 
             chains.log10N[,lower.cols,drop=FALSE])
    sum.chain.length <- apply(l, 1, sum, na.rm=TRUE)

    res <- cbind(chain.span, count.chain.length, sum.chain.length)
    rownames(res) <- NULL
    return (res)
}

.PolygonArea <- function(x, y)
{
    # htncp://en.wikipedia.org/wiki/Polygon#Area_and_centroid
    stopifnot(length(x)==length(y))
    upper <- 2:length(x)
    lower <- upper-1
    return (abs(0.5 * sum(x[lower]*y[upper] - x[upper]*y[lower])))
}

.ConvexHull <- function(community, x, y)
{
    p <- cbind(x, y)
    p <- p[!is.na(p[,1]) & !is.na(p[,2]),]
    hull <- chull(p)

    # Close polygon for calculation of area
    return (list(nodes=unname(NP(community, 'node')[hull]),
                 points=p[hull,], 
                 area=.PolygonArea(p[c(hull,hull[1]),1], p[c(hull,hull[1]),2])))
}

NvMConvexHull <- function(community)
{
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireN(community)
    return (.ConvexHull(community, Log10M(community), Log10N(community)))
}
