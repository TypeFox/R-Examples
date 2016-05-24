# Functions requiring M and/or N
TestBodyMassBins <- function()
{
    AssertEqual(c(1,3,10), BodyMassBins(c6), check.attributes=FALSE)
    AssertEqual(c(1,1,3), BodyMassBins(c6, n.bins=3), check.attributes=FALSE)
    AssertEqual(c(2,2,2,1,2,1,1,5,2,2,1,1,2,3,2,3,3,2,2,1,4,3,2,3,1,2,1,2,3,
                    4,1,4,5,5,3,6,5,6,4,4,4,7,3,3,3,6,6,4,4,4,4,6,7,10,10,10), 
                  BodyMassBins(TL84), check.attributes=FALSE)

    AssertRaises(BodyMassBins(c1))
    AssertRaises(BodyMassBins(c2))
    AssertRaises(BodyMassBins(TL84, n.bins=1))
    AssertRaises(BodyMassBins(TL84, upper=min(NP(TL84,'M')), 
                              lower=max(NP(TL84,'M'))))
    AssertRaises(BodyMassBins(TL84, lower=mean(NP(TL84,'M')), 
                              upper=max(NP(TL84,'M'))))
    AssertRaises(BodyMassBins(TL84, lower=min(NP(TL84,'M')), 
                              upper=mean(NP(TL84,'M'))))
}

TestResourceLargerThanConsumer <- function()
{
    AssertEqual(NULL, ResourceLargerThanConsumer(c1))
    AssertEqual(NULL, ResourceLargerThanConsumer(c2))
    AssertEqual(0, nrow(ResourceLargerThanConsumer(c6)))
    check <- ResourceLargerThanConsumer(TL84)
    AssertEqual(c('Cyclops varians rubellus', 'Leptodiaptomus siciloides'), 
                  check$resource)
    AssertEqual(c('Tropocyclops prasinus', 'Tropocyclops prasinus'), 
                  check$consumer)
}

TestSumNByClass <- function()
{
    AssertRaises(SumNByClass(c1))
    AssertRaises(SumNByClass(c2))
    AssertEqual(unname(NP(c6,'N')), 
                unname(SumNByClass(c6, class=c('a','b','c'))))
    AssertEqual(sum(NP(c6,'N')), unname(SumNByClass(c6, class=c('a','a','a'))))
    AssertEqual(sum(NP(c6,'N')[1:2]), 
                unname(SumNByClass(c6, class=c('a','a',NA))))
    AssertEqual(unname(c(NP(c6,'N')[3], sum(NP(c6,'N')[1:2]))), 
                unname(SumNByClass(c6, class=c('a','a',''))))

    AssertEqual(c(invertebrate=1.8678e+06, producer=3.2107e+09, 
                  vert.ecto=2.2350e+00), SumNByClass(TL84))

    AssertEqual(as.numeric(c(NA,NA,NA)), 
                unname(SumNByClass(BroadstoneStream)))
    AssertEqual(c('<unnamed>'=0, invertebrate=32081.3, producer=0.0), 
                SumNByClass(BroadstoneStream, na.rm=TRUE))

    AssertRaises(SumNByClass(c6, class=NULL))
}

TestSumBiomassByClass <- function()
{
    AssertRaises(SumNByClass(c1))
    AssertRaises(SumNByClass(c2))
    AssertEqual(unname(Biomass(c6)), 
                unname(SumBiomassByClass(c6, class=c('a','b','c'))))
    AssertEqual(sum(Biomass(c6)), 
                unname(SumBiomassByClass(c6, class=c('a','a','a'))))
    AssertEqual(sum(Biomass(c6)[1:2]), 
                unname(SumBiomassByClass(c6, class=c('a','a',NA))))
    AssertEqual(unname(c(Biomass(c6)[3], sum(Biomass(c6)[1:2]))), 
                unname(SumBiomassByClass(c6, class=c('a','a',''))))

    AssertEqual(as.numeric(c(NA,NA,NA)), 
                unname(SumBiomassByClass(BroadstoneStream)))
    AssertEqual(c('<unnamed>'=0, invertebrate=830.4337, producer=0), 
                round(SumBiomassByClass(BroadstoneStream, na.rm=TRUE),4))

    AssertRaises(SumBiomassByClass(c6, class=NULL))
}

TestNvMLinearRegressions <- function()
{
    res <- NvMLinearRegressions(c1)
    AssertEqual('all', names(res))
    AssertEqual(NULL, res[[1]])
    res <- NvMLinearRegressions(c2)
    AssertEqual('all', names(res))
    AssertEqual(NULL, res[[1]])

    # Using default class of category - all NA as c6 does not have category
    res <- NvMLinearRegressions(c6)
    AssertEqual('all', names(res))
    AssertEqual(c(2.58,-1.04), unname(round(coef(res[['all']]), 2)))

    # class, NULL should result in a single lm object being returned
    res <- NvMLinearRegressions(c6, class=NULL)
    AssertEqual('all', names(res))
    AssertEqual(c(2.58,-1.04), unname(round(coef(res[['all']]), 2)))

    # Each node is in it's own class. Can't fit lm through one data point 
    # so should have 'all' and three NULLs
    res <- NvMLinearRegressions(c6, class=c('a','b','c'))
    AssertEqual(c('all','a','b','c'), names(res))
    AssertEqual(c(2.58,-1.04), unname(round(coef(res[['all']]), 2)))
    AssertEqual(list(a=NULL, b=NULL, c=NULL), res[2:4])

    # Node 1 in a different class to nodes 2 and 3
    res <- NvMLinearRegressions(c6, class=c('a','b','b'))
    AssertEqual(c('all', 'a', 'b'), names(res))
    AssertEqual(c(2.58,-1.04), unname(round(coef(res[['all']]), 2)))
    AssertEqual(c(1.14,-0.20), unname(round(coef(res[['b']]), 2)))
    AssertEqual(NULL, res[['a']])

    # category is default class
    res <- NvMLinearRegressions(TL84)
    AssertEqual(c('all', 'producer', 'invertebrate', 'vert.ecto'), names(res))
    AssertEqual(c(-2.69, -0.83), unname(round(coef(res[['all']]), 2)))
    AssertEqual(c(2.56, -0.41), unname(round(coef(res[['producer']]), 2)))
    AssertEqual(c(1.47, -0.32), unname(round(coef(res[['invertebrate']]), 2)))
    AssertEqual(c(-34.66, -11.63), unname(round(coef(res[['vert.ecto']]), 2)))

    # Some nodes lack both a category and N and M
    res <- NvMLinearRegressions(BroadstoneStream)
    AssertEqual(c(1.08, -0.74), unname(round(coef(res[['all']]), 2)))
    AssertEqual(c(1.08, -0.74), unname(round(coef(res[['invertebrate']]), 2)))
    AssertEqual(NULL, res[['<unnamed>']])
    AssertEqual(NULL, res[['producer']])
}

TestNvMTriTrophic1 <- function()
{
    # Recreates Table 1 from Cohen et al 2010 PNAS
    data(TL84, TL86, YthanEstuary)
    collection <- CommunityCollection(list(TL84, TL86, YthanEstuary))
    table <- NvMTriTrophicTable(collection)

    # Exclude the rows that include all nodes and links
    table <- head(table, -4)

    # Check to 2 dp
    table <- round(table, 2)

    # The data from the paper
    check <- matrix(c(
               6.33,    5.90,    7.29, 
               5.41,    3.43,    5.06, 
               5.99,    5.69,    6.15, 
              12.67,   11.79,   14.57, 
              11.02,    8.65,   10.51, 
              11.40,    9.12,   11.20, 
               1.15,    1.36,    1.39, 
               1.03,    1.05,    1.07, 
               4.84,    4.84,    4.43, 
              30.62,   28.56,   32.31, 
              20.78,   22.66,   21.98, 
               1.47,    1.26,    1.47, 
              19.96,   23.33,   16.88, 
              18.71,   20.63,   13.18,  # TL86 mean chain span 20.63 vs 20.62
               0.90,    0.91,    0.60, 
               1.07,    1.13,    1.28, 
               0.96,    1.03,    0.77, 
             264.00,  236.00,  379.00, 
            2500.00, 2601.00, 8281.00, 
               0.11,    0.09,    0.05, 
               5.28,    4.63,    4.16), ncol=3, byrow=TRUE)
    dimnames(check) <- dimnames(table)
    AssertEqual(table, check)
}

TestNvMTriTrophic2 <- function()
{
    # Recreates Table S3 from Cohen et al 2010 PNAS
    communities <- list(TL84, TL86, RemoveNodes(YthanEstuary,'POM (detritus)'))
    communities <- lapply(communities, RemoveIsolatedNodes)
    communities <- lapply(communities, RemoveCannibalisticLinks)

    res <- lapply(communities, function(community)
    {
        chains <- ThreeNodeChains(community, node.properties='M')
        MR <- chains$bottom.M
        MI <- chains$intermediate.M
        MC <- chains$top.M

        lp <- TLPS(community, node.properties='M')

        return (c('MR<=MI<=MC'=sum(MR<=MI & MI<=MC), 
                  'MR<=MC<MI'=sum(MR<=MC & MC<MI), 
                  'MI<MR<=MC'=sum(MI<MR  & MR<=MC), 
                  'MI<=MC<MR'=sum(MI<=MC & MC<MR), 
                  'MC<MR<MI'=sum(MC<MR  & MR<MI), 
                  'MC<MI<MR'=sum(MC<MI  & MI<MR), 
                  'All 2-chains'=nrow(chains),
                  'MR<MC'=sum(lp$resource.M<lp$consumer.M), 
                  'MR=MC'=sum(lp$resource.M==lp$consumer.M), 
                  'MR>MC'=sum(lp$resource.M>lp$consumer.M), 
                  'All links'=nrow(lp)))
    })
    res <- do.call('cbind', res)

    colnames(res) <- c('TL84', 'TL86', 'Ythan Estuary')

    # The data from the paper
    check <- matrix(c(  1001,  577,  1232,
                          30,   59,    65,
                          12,   10,    68,
                           0,    1,     3,
                           1,    3,     0,
                           0,    1,     3,
                        1044,  651,  1371,
                         262,  232,   368,
                           0,    0,     2,
                           2,    4,     9,
                         264,  236,   379), ncol=3, byrow=TRUE)

    colnames(check) <- colnames(res)
    rownames(check) <- rownames(res)
    AssertEqual(check, res)
}

TestNvMTriTrophic3 <- function()
{
    # Some of Doris' Icelandic streams communities have only producers and 
    # herbivores, which caused an earlier version of 
    # .NvMTrophicChainProperties() to fail with a nasty error.
    # Let's check some cases for simple communities
    AssertRaises(NvMTriTrophicStatistics(c1))
    AssertRaises(NvMTriTrophicStatistics(c2))
    AssertRaises(NvMTriTrophicStatistics(c3))
    AssertRaises(NvMTriTrophicStatistics(c4))
    AssertRaises(NvMTriTrophicStatistics(c5))

    # Resource-consumer
    rc <- Community(nodes=data.frame(node=c('R', 'C'),
                                     M=c(1, 10), 
                                     N=c(100, 1)), 
                    trophic.links=data.frame(resource='R', consumer='C'), 
                    properties=list(title='Resource-consumer', M.units='g', 
                                    N.units='m^-2'))
    tts <- NvMTriTrophicStatistics(rc)
    check <- data.frame(resource='R',
                        consumer='C',
                        length=3,
                        angle=-63.43494882292200998108,
                        slope=-2, 
                        stringsAsFactors=FALSE, 
                        row.names="1")
    AssertEqual(check, tts$links)

    AssertEqual(NULL, tts$three.node.chains)

    check <- data.frame(Node.1='R',
                        Node.2='C', 
                        chain.span=3,
                        count.chain.length=1,
                        sum.chain.length=3,
                        stringsAsFactors=FALSE)
    AssertEqual(check, tts$trophic.chains)

    # Tri-trophic chain
    tts <- NvMTriTrophicStatistics(c6)
    check <- data.frame(resource=c('R','C'),
                        consumer=c('C','P'),
                        length=c(2.52287874528033739807,1.56066730616973736723),
                        angle=c(-75.34856290594565564334,-11.28584933963889191944),
                        slope=c(-3.82497857878639635487, -0.19956289353132869446),
                        stringsAsFactors=FALSE, 
                        row.names=c("1","2"))
    AssertEqual(check, tts$links)

    check <- data.frame(bottom='R', 
                        intermediate='C', 
                        top='P',
                        Llower=2.52287874528033739807,
                        Lupper=1.56066730616973736723, 
                        two.span=4.08354605145007454325,
                        Alower=-75.34856290594565564334,
                        Aupper=-11.28584933963889191944,
                        Abetween=64.06271356630675484212,
                        stringsAsFactors=FALSE)
    AssertEqual(check, tts$three.node.chains)

    check <- data.frame(Node.1='R',
                        Node.2='C', 
                        Node.3='P', 
                        chain.span=4.08354605145007454325,
                        count.chain.length=2.00000000000000000000,
                        sum.chain.length=4.08354605145007454325,
                        stringsAsFactors=FALSE)
    AssertEqual(check, tts$trophic.chains)
}

TestNodesWithoutMOrN <- function()
{
    NodesWithoutMOrN <- cheddar:::.NodesWithoutMOrN
    AssertEqual('S', NodesWithoutMOrN(c1))
    AssertEqual('S', NodesWithoutMOrN(c2))
    AssertEqual(c('R','C'), NodesWithoutMOrN(c3))
    AssertEqual(c('R','C','P'), NodesWithoutMOrN(c4))
    AssertEqual(c('R','C','O'), NodesWithoutMOrN(c5))
    AssertEqual(0, length(NodesWithoutMOrN(c6)))
    AssertEqual(0, length(NodesWithoutMOrN(TL84)))
    AssertEqual('POM (detritus)', NodesWithoutMOrN(YthanEstuary))
}

TestNvMSlopeAndIntercept <- function()
{
    AssertEqual(NULL, NvMSlope(c1))
    AssertEqual(NULL, NvMSlope(c2))
    AssertEqual(NULL, NvMSlope(c3))
    AssertEqual(NULL, NvMSlope(c4))
    AssertEqual(NULL, NvMSlope(c5))
    AssertEqual(-1.04009302605305009592, NvMSlope(c6))
    AssertEqual(-0.82711453844476423569, NvMSlope(TL84))

    AssertEqual(NULL, NvMIntercept(c1))
    AssertEqual(NULL, NvMIntercept(c2))
    AssertEqual(NULL, NvMIntercept(c3))
    AssertEqual(NULL, NvMIntercept(c4))
    AssertEqual(NULL, NvMIntercept(c5))
    AssertEqual(2.57689795300774049380, NvMIntercept(c6))
    AssertEqual(-2.68627529180856106095, NvMIntercept(TL84))

    AssertEqual(NULL, NvMSlopeAndIntercept(c1))
    AssertEqual(NULL, NvMSlopeAndIntercept(c2))
    AssertEqual(NULL, NvMSlopeAndIntercept(c3))
    AssertEqual(NULL, NvMSlopeAndIntercept(c4))
    AssertEqual(NULL, NvMSlopeAndIntercept(c5))
    AssertEqual(c(slope=-1.04009302605305009592, 
                  intercept=2.57689795300774049380), 
                NvMSlopeAndIntercept(c6))
    AssertEqual(c(slope=-0.82711453844476423569, 
                  intercept=-2.68627529180856106095), 
                NvMSlopeAndIntercept(TL84))
}

TestNvMSlopeAndInterceptByClass <- function()
{
    AssertEqual(c(slope.all=-1.04009302605305009592), NvMSlopeByClass(c6))
    AssertEqual(c(slope.all=-0.82711453844476423569, 
                  slope.producer=-0.40715089579609509141, 
                  slope.invertebrate=-0.32431675051076519489, 
                  slope.vert.ecto=-11.62787086769622924010), 
                NvMSlopeByClass(TL84))

    AssertEqual(c(intercept.all=2.57689795300774049380),NvMInterceptByClass(c6))
    AssertEqual(c(intercept.all=-2.68627529180856106095, 
                  intercept.producer=2.55833642234355362888, 
                  intercept.invertebrate=1.46560619016986248830, 
                  intercept.vert.ecto=-34.66097278952445748246), 
                NvMInterceptByClass(TL84))

    AssertEqual(c(slope.all=-1.04009302605305009592, 
                  intercept.all=2.57689795300774049380), 
                NvMSlopeAndInterceptByClass(c6))
    AssertEqual(c(slope.all=-0.82711453844476423569, 
                  slope.producer=-0.40715089579609509141, 
                  slope.invertebrate=-0.32431675051076519489, 
                  slope.vert.ecto=-11.62787086769622924010, 
                  intercept.all=-2.68627529180856106095, 
                  intercept.producer=2.55833642234355362888, 
                  intercept.invertebrate=1.46560619016986248830, 
                  intercept.vert.ecto=-34.66097278952445748246), 
                NvMSlopeAndInterceptByClass(TL84))

    AssertEqual(c(slope.all=-0.73916494129624910059, 
                  slope.invertebrate=-0.73916494129624910059, 
                  'slope.<unnamed>'=NA,
                  slope.producer=NA,
                  intercept.all=1.07951558005480174884, 
                  intercept.invertebrate=1.07951558005480174884,
                  'intercept.<unnamed>'=NA,
                  intercept.producer=NA), 
                NvMSlopeAndInterceptByClass(BroadstoneStream))

    AssertEqual(c(slope.all=-1.13010867280292459647,
                  slope.vert.endo=-2.46246611880339116851, 
                  'slope.<unnamed>'=-1.87594996135075309240, 
                  slope.vert.ecto=2.32192809488736218171,
                  intercept.all=1.43776079409099510897, 
                  intercept.vert.endo=2.95188840529695939452, 
                  'intercept.<unnamed>'=1.16171181260177314165, 
                  intercept.vert.ecto=0.00000000000000005773), 
                NvMSlopeAndInterceptByClass(c8))
}

TestNvMConvexHull <- function()
{
    AssertEqual(0.5, cheddar:::.PolygonArea(c(0, 1, 0.5, 0), c(0, 0, 1, 0)))
    AssertEqual(1, cheddar:::.PolygonArea(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))
    res <- NvMConvexHull(TL84)
    # Can't guarantee in ordering of points on the hull
    nodes <- c('Chaoborus punctipennis','Chromulina sp.',
               'Chrysosphaerella longispina','Keratella testudo','Phoxinus eos',
               'Phoxinus neogaeus','Umbra limi','Unclassified flagellates')
    AssertEqual(nodes, sort(res$nodes))
    AssertEqual(30.68266, round(res$area, 5))
    AssertEqual(matrix(c( -6.522879,  4.0791812,
                         -13.518557,  8.1731863,
                          -9.080399,  6.6020600,
                         -11.000000,  3.7781513,
                          -2.995679,  0.2944662,
                          -2.931814, -0.8761484,
                          -2.889410, -0.8794261,
                         -12.460924,  9.2741578), byrow=TRUE, ncol=2), 
                unname(round(res$points[nodes,], 7)), 
                tolerance=1e-7)
}
