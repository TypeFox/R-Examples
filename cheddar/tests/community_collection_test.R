# Community collections
TestCommunityCollectionSingle <- function()
{
    # Test collections containing a single community
    for(community in list(c1,c2,c3,c4,c5,c6,TL84,TL86,SkipwithPond))
    {
        print(community)
        collection <- CommunityCollection(list(community))
        AssertEqual(1, length(collection))
        AssertEqual(CP(community,'title'), names(collection))
        AssertEqual(data.frame(summary(community), 
                               row.names=CPS(community)$title), 
                    summary(collection))

        AssertEqual(collection, CommunityCollection(collection))

        cpa <- CollectionCPS(collection)
        cpb <- data.frame(CPS(community), stringsAsFactors=FALSE)
        rownames(cpb) <- cpb$title
        AssertEqual(cpa, cpb)

        npa <- CollectionNPS(collection)
        npb <- cbind(community=CP(community,'title'), NPS(community))
        rownames(npb) <- NULL
        AssertEqual(npa, npb)
    
        if('c1'==CP(community, 'title'))
        {
            AssertEqual(NULL, CollectionTLPS(collection))
        }
        else
        {
            cpa <- CollectionTLPS(collection)
            cpb <- cbind(community=CP(community,'title'), TLPS(community))
            AssertEqual(cpa, cpb)
        }
    }
}

TestCommunityCollectionProperties <- function()
{
    # TODO Test these cases
    col <- CommunityCollection(list(TL84, TL86, YthanEstuary))

    # Functions that return a single value
    CollectionCPS(col)
    CollectionCPS(col, 'title')
    CollectionCPS(col, c(S='NumberOfNodes'))
    CollectionCPS(col, c(S='NumberOfNodes', L='NumberOfTrophicLinks'))

    # A function that return more than one value
    CollectionCPS(col, 'SumBiomassByClass')
    CollectionCPS(col, c(B='SumBiomassByClass'))
    CollectionCPS(col, list(list('SumBiomassByClass', class='kingdom')))
    CollectionCPS(col, list(B=list('SumBiomassByClass', class='kingdom')))
    CollectionCPS(col, list(S='NumberOfNodes', B=list('SumBiomassByClass', class='kingdom'), L='NumberOfTrophicLinks'))

    # Functions that return more than one value
    CollectionCPS(col, c('SumBiomassByClass', 'SumNByClass'))
    CollectionCPS(col, c(B='SumBiomassByClass', 'SumNByClass'))
    CollectionCPS(col, list(list('SumBiomassByClass', class='kingdom'), 'SumNByClass'))
    CollectionCPS(col, list(B=list('SumBiomassByClass', class='kingdom'), N='SumNByClass'))
    CollectionCPS(col, list(B=list('SumBiomassByClass', class='kingdom'), S='NumberOfNodes', N='SumNByClass'))
    CollectionCPS(col, list(B=list('SumBiomassByClass', class='kingdom'), N=list('SumNByClass', class='kingdom')))
    CollectionCPS(col, list(B=list('SumBiomassByClass', class='kingdom'), S='NumberOfNodes', N=list('SumNByClass', class='kingdom')))

    head(CollectionNPS(col))
    head(CollectionNPS(col, 'Biomass'))
    head(CollectionNPS(col, c(B='Biomass')))
    head(CollectionNPS(col, c(B='Biomass', 'M')))
    head(CollectionNPS(col, 'Log10MNBiomass'))

    head(CollectionTLPS(col))
    head(CollectionTLPS(col, node.properties='Biomass'))
    head(CollectionTLPS(col, node.properties=c(B='Biomass')))
    head(CollectionTLPS(col, node.properties=c(B='Biomass', 'M')))

    S2 <- Community(properties=list(title='Skipwith (no links)'), 
                                    nodes=NPS(SkipwithPond))
    col <- CommunityCollection(list(SkipwithPond, S2))
    head(CollectionTLPS(col))
    head(CollectionTLPS(col, link.properties='link.evidence'))
    NULL
}

TestCommunityCollectionFailures <- function()
{
    # Not communities
    AssertRaises(CommunityCollection(list()))
    AssertRaises(CommunityCollection(1))
    AssertRaises(CommunityCollection(''))
    AssertRaises(CommunityCollection(list(c1,'')))
    AssertRaises(CommunityCollection(list(c1,1)))

    # Duplications
    AssertRaises(CommunityCollection(list(c1,c1)))
    AssertRaises(CommunityCollection(list(c1,c2,c1)))

    # Inconsistent units
    AssertRaises(CommunityCollection(list(c5, c6)))

    # Modifications are illegal
    collection <- CommunityCollection(list(TL84))
    AssertRaises(collection[1] <- 1)
    AssertRaises(collection[[1]] <- 1)
    AssertRaises(collection$silly <- 1)
    AssertRaises(names(collection) <- letters[1:3])
    AssertRaises(length(collection) <- 1)
}

TestAggregateCommunitiesSingle <- function()
{
    # Aggregating a collection containing one community should not change 
    # anything
    for(community in list(c1,c2,c3,c4,c5,c6,TL84,TL86,SkipwithPond))
    {
        print(community)
        a <- community
        b <- AggregateCommunities(CommunityCollection(list(community)), 
                                  title=CP(community,'title'))

        # Node, trophic link and whole-community properties should not be 
        # different.
        npa <- NPS(a)
        npa <- npa[order(npa$node),]
        npb <- NPS(b)
        npb <- npb[order(npb$node),]
        AssertEqual(npa, npb)

        cpsa <- CPS(a)
        cpsa <- cpsa[order(names(cpsa))]
        cpsb <- CPS(b)
        cpsb <- cpsb[order(names(cpsb))]
        AssertEqual(cpsa, cpsb)

        if('c1'==CP(a, 'title'))
        {
            AssertEqual(0, NumberOfTrophicLinks(a))
            AssertEqual(0, NumberOfTrophicLinks(b))
            AssertEqual(NULL, TLPS(a))
            AssertEqual(NULL, TLPS(b))
        }
        else
        {
            tlpa <- TLPS(a)
            tlpa <- tlpa[order(tlpa$resource, tlpa$consumer),]
            rownames(tlpa) <- NULL
            tlpb <- TLPS(b)
            tlpb <- tlpb[order(tlpb$resource, tlpb$consumer),]
            rownames(tlpb) <- NULL
            AssertEqual(tlpa, tlpb)
        }
    }
}

TestAggregateCommunitiesFailures <- function()
{
    # Some failure cases
    collection <- CommunityCollection(list(c1,c2,c3,c4,c5))
    AssertRaises(AggregateCommunities(collection, 0))
    AssertRaises(AggregateCommunities(collection, 6))
    AssertRaises(AggregateCommunities(collection, 0:1))
    AssertRaises(AggregateCommunities(collection, 1:6))
    AssertRaises(AggregateCommunities(collection, ''))
    AssertRaises(AggregateCommunities(collection, c('c1','x')))
}

TestAggregateCommunitiesProperties <- function()
{
    # Aggregation of three very simple communities with no node or link
    # properties and only the 'title' community property
    collection <- CommunityCollection(list(c1,c2,c3))
    a <- AggregateCommunities(collection)
    AssertEqual(data.frame(node=c('S','R','C'), row.names=c('S','R','C')),NPS(a))
    AssertEqual(data.frame(resource=c('S','R'), consumer=c('S','C')), TLPS(a))
    AssertEqual('Aggregation of c1,c2,c3', CPS(a)$title)

    # Aggregation of two simple communities with a single node, a single link 
    # and a single community property
    x <- Community(nodes=data.frame(node='n', np='np'),
                   trophic.links=data.frame(resource='n', consumer='n', tlp='tlp'),
                   properties=list(title='x', cp='cp'))
    y <- Community(nodes=data.frame(node='n', np='np'),
                   trophic.links=data.frame(resource='n', consumer='n', tlp='tlp'),
                   properties=list(title='y', cp='cp'))
    collection <- CommunityCollection(list(x,y))
    a <- AggregateCommunities(collection)
    AssertEqual(data.frame(node='n', np='np', row.names='n'), NPS(a))
    AssertEqual(data.frame(resource='n', consumer='n', tlp='tlp'), TLPS(a))
    AssertEqual('Aggregation of x,y',  CPS(a)$title)
    AssertEqual('cp', CPS(a)$cp)
    
    # Aggregation of two simple communities with two node, two link 
    # and two community properties
    x <- Community(nodes=data.frame(node='n', np1='npa', np2='npb'),
                   trophic.links=data.frame(resource='n', consumer='n',
                                            tlp1='tlpa', tlp2='tlpb'),
                   properties=list(title='x', cp1='cpa', cp2='cpb'))
    y <- Community(nodes=data.frame(node='n', np1='npc'),
                   trophic.links=data.frame(resource='n', consumer='n',
                                            tlp2='tlpc'),
                   properties=list(title='y', cp1='cpc'))
    collection <- CommunityCollection(list(x,y))
    a <- AggregateCommunities(collection)
    AssertEqual(data.frame(node='n', np1='npa,npc', np2='npb', row.names='n'), NPS(a))
    AssertEqual(data.frame(resource='n', consumer='n', tlp1='tlpa', tlp2='tlpb,tlpc'), TLPS(a))
    AssertEqual('Aggregation of x,y', CPS(a)$title)
    AssertEqual('cpa,cpc', CPS(a)$cp1)
    AssertEqual('cpb', CPS(a)$cp2)

    # Aggregation of properties for TL84 and TL86
    a <- CommunityCollection(list(TL84, TL86))

    # tapply coerces INDEX parameter to a factor, sorting alphabetically. This 
    # line retains the existing ordering.
    a.nodes <- factor(CollectionNPS(a)[,'node'], 
                      levels=unique(CollectionNPS(a)[,'node']))

    b <- AggregateCommunities(a)

    # Sanity check
    all.spp <- sort(unique(c(NP(TL84,'node'), NP(TL86,'node'))))
    AssertEqual(all.spp, sort(unique(CollectionNPS(a)[,'node'])))
    AssertEqual(unname(sort(NP(b, 'node'))), all.spp)

    # Check category
    AssertEqual(unname(NP(b,'category')),
                as.vector(tapply(CollectionNPS(a)[,'category'], a.nodes, unique)))

    # TODO Test that union of trophic links taken
}

TestAggregateCommunityMeans <- function()
{
    # Cases for aggregating numeric values

    CheckResults <- function(A, B, expected.M, expected.N, weight.by)
    {
        collection <- CommunityCollection(list(A, B))
        res <- NPS(AggregateCommunities(collection, weight.by=weight.by))
        check <- data.frame(node=c('S1','S2'), M=expected.M, N=expected.N, 
                            row.names=c('S1','S2'), stringsAsFactors=FALSE)
        AssertEqual(check, res)
    }

    # 1. Both species are in both communities. Valid M and N in both communities.
    #  community node  M N
    #         A   S1 20 1
    #         A   S2 20 1
    #         B   S1  5 2
    #         B   S2  5 2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c(1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5, 5), N=c(2, 2)))

    # Weight by N
    #   node  M   N
    #     S1 10 1.5
    #     S2 10 1.5
    CheckResults(A, B, c(10,10), c(1.5, 1.5), 'N')

    # No weighting
    #   node  M   N
    #     S1 12.5 1.5
    #     S2 12.5 1.5
    CheckResults(A, B, c(12.5,12.5), c(1.5, 1.5), NULL)


    # 2. Both species are in both communities. Species 2 lacks N in community B.
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 20  1
    #         B   S1  5  2
    #         B   S2  5 NA
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c( 1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5, 5), N=c( 2,NA)))

    # Weight by N
    #   node  M   N
    #     S1 10  1.5
    #     S2 NA  NA
    CheckResults(A, B, c(10, NA), c(1.5, NA), 'N')

    # No weighting
    #   node  M    N
    #     S1 12.5  1.5
    #     S2 12.5  NA
    CheckResults(A, B, c(12.5, 12.5), c(1.5, NA), NULL)


    # 3. Both species are in both communities. Species 2 lacks M in community B
    #  community node  M N
    #         A   S1 20 1
    #         A   S2 20 1
    #         B   S1  5 2
    #         B   S2 NA 2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c( 1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5,NA), N=c( 2, 2)))

    # Weight by N
    #   node  M   N
    #   S1   10  1.5
    #   S2   NA  1.5
    CheckResults(A, B, c(10, NA), c(1.5, 1.5), 'N')

    # No weighting
    #   node  M   N
    #   S1 12.5  1.5
    #   S2   NA  1.5
    CheckResults(A, B, c(12.5, NA), c(1.5, 1.5), NULL)


    # 4. Both species are in both communities. Species 2 lacks M and N in 
    #    community B
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 20  1
    #         B   S1  5  2
    #         B   S2 NA NA
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c( 1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5,NA), N=c( 2,NA)))

    # Weight by N
    #   node  M   N
    #   S1   10  1.5
    #   S2   NA  NA
    CheckResults(A, B, c(10, NA), c(1.5, NA), 'N')

    # No weighting
    #   node  M   N
    #   S1 12.5  1.5
    #   S2   NA  NA
    CheckResults(A, B, c(12.5, NA), c(1.5, NA), NULL)

    # 5. Species 2 is in community A only. Valid M and N.
    #  community node  M N
    #         A   S1 20 1
    #         A   S2 20 1
    #         B   S1  5 2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c(1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1'),      M=c( 5),    N=c(2)))

    # Weight by N
    #   node  M   N
    #   S1   10 1.5
    #   S2   20 0.5
    CheckResults(A, B, c(10, 20), c(1.5, 0.5), 'N')

    # No weighting
    #   node  M   N
    #   S1  12.5 1.5
    #   S2  10   0.5
    CheckResults(A, B, c(12.5, 10), c(1.5, 0.5), NULL)


    # 6. Species 2 is in community A only and lacks N.
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 20 NA
    #         B   S1  5  2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c(1, NA)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1'),      M=c( 5),    N=c(2)))
    CollectionNPS(CommunityCollection(list(A, B)))
    NPS(AggregateCommunities(CommunityCollection(list(A, B))))

    # Weight by N
    #   node  M   N
    #   S1   10 1.5
    #   S2   NA  NA
    CheckResults(A, B, c(10, NA), c(1.5, NA), 'N')

    # No weighting
    #   node  M   N
    #   S1  12.5 1.5
    #   S2  10   NA
    CheckResults(A, B, c(12.5, 10), c(1.5, NA), NULL)


    # 7. Species 2 is in community A only and lacks M.
    #  community node M N
    #         A   S1 20 1
    #         A   S2 NA 1
    #         B   S1  5 2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,NA), N=c(1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1'),      M=c( 5),    N=c(2)))
    CollectionNPS(CommunityCollection(list(A, B)))
    NPS(AggregateCommunities(CommunityCollection(list(A, B))))

    # Weight by N
    #   node  M   N
    #   S1   10  1.5
    #   S2   NA  0.5
    CheckResults(A, B, c(10, NA), c(1.5, 0.5), 'N')

    # No weighting
    #   node  M   N
    #   S1   12.5 1.5
    #   S2   NA   0.5
    CheckResults(A, B, c(12.5, NA), c(1.5, 0.5), NULL)


    # 8. Species 2 is in community A only and lacks M and N.
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 NA NA
    #         B   S1  5  2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,NA), N=c(1,NA)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1'),      M=c( 5),    N=c(2)))
    CollectionNPS(CommunityCollection(list(A, B)))
    NPS(AggregateCommunities(CommunityCollection(list(A, B))))

    # Weight by N
    #   node  M   N
    #     S1  10  1.5
    #     S2  NA  NA
    CheckResults(A, B, c(10, NA), c(1.5, NA), 'N')

    # No weighting
    #   node  M   N
    #   S1   12.5 1.5
    #   S2   NA   NA
    CheckResults(A, B, c(12.5, NA), c(1.5, NA), NULL)


    # 9. Both species are in both communities. Species 2 has missing M in both 
    #    communities
    #  community node  M N
    #         A   S1 20 1
    #         A   S2 NA 1
    #         B   S1  5 2
    #         B   S2 NA 2
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,NA), N=c(1, 1)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5,NA), N=c(2, 2)))

    # Weight by N
    #   node  M   N
    #   S1  10  1.5
    #   S2  NA  1.5
    CheckResults(A, B, c(10, NA), c(1.5, 1.5), 'N')

    # No weighting
    #   node  M   N
    #   S1   12.5 1.5
    #   S2   NA   1.5
    CheckResults(A, B, c(12.5, NA), c(1.5, 1.5), NULL)


    # 10. Both species are in both communities. Species 2 has missing N in both 
    #     communities
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 20 NA
    #         B   S1  5  2
    #         B   S2  5 NA
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,20), N=c(1,NA)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5, 5), N=c(2,NA)))

    # Weight by N
    #   node  M   N
    #    S1  10  1.5
    #    S2  NA  NA
    CheckResults(A, B, c(10, NA), c(1.5, NA), 'N')

    # No weighting
    #   node  M   N
    #   S1   12.5 1.5
    #   S2   12.5 NA
    CheckResults(A, B, c(12.5, 12.5), c(1.5, NA), NULL)


    # 11. Both species are in both communities. Species 2 has missing M and N in 
    #     both communities
    #  community node  M  N
    #         A   S1 20  1
    #         A   S2 NA NA
    #         B   S1  5  2
    #         B   S2 NA NA
    A <- Community(properties=list(title='A', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c(20,NA), N=c(1,NA)))
    B <- Community(properties=list(title='B', M.units='g', N.units='m^-2'), 
                   nodes=data.frame(node=c('S1','S2'), M=c( 5,NA), N=c(2,NA)))

    # Weight by N
    #   node  M   N
    #    S1  10  1.5
    #     S2 NA  NA
    CheckResults(A, B, c(12.5, NA), c(1.5, NA), NULL)

    # No weighting
    #   node  M   N
    #   S1   12.5 1.5
    #   S2   NA   NA
    CheckResults(A, B, c(12.5, NA), c(1.5, NA), NULL)
}

TestAggregateCommunitiesBy <- function()
{
    # TODO Test these cases
    AggregateCommunitiesBy(pHWebs, 'pH')
    AggregateCommunitiesBy(pHWebs, 'lat')
    AggregateCommunitiesBy(pHWebs, 'NumberOfTrophicLinks')
    NULL
}

TestPersistCollection <- function()
{
    path <- tempfile()
    on.exit(unlink(path, recursive=TRUE))
    for(collection in list(pHWebs, Millstream))
    {
        print(collection)
        SaveCollection(collection, path)
        loaded <- LoadCollection(path)[names(collection)]
        AssertEqual(collection, loaded)
        unlink(path, recursive=TRUE)

        SaveCollection(collection, path, fn='write.table', sep='\t')
        loaded <- LoadCollection(path, fn='read.table', sep='\t')
        loaded <- loaded[names(collection)]
        AssertEqual(collection, loaded)
        unlink(path, recursive=TRUE)
    }
}

TestSiteBySpeciesMatrix <- function()
{
    res <- SiteBySpeciesMatrix(Millstream)
    AssertEqual(c(67,2), dim(res))
    AssertEqual(117, sum(res))
    AssertTrue(all(res %in% 0:1))
    AssertEqual('Valvata piscinalis', rownames(res)[67])
    AssertEqual(c(c4=1, d4=1), res[1,])
    AssertEqual(c(c4=1, d4=0), res[67,])

    # Same data but missing species represented by NA
    res <- SiteBySpeciesMatrix(Millstream, na.missing=TRUE)
    AssertEqual(c(67,2), dim(res))
    AssertEqual(117, sum(res, na.rm=TRUE))
    AssertTrue(all(res %in% c(NA,1)))
    AssertEqual('Valvata piscinalis', rownames(res)[67])
    AssertEqual(c(c4=1, d4=1), res[1,])
    AssertEqual(c(c4=1, d4=NA), res[67,])

    # log10 biomass abundance
    res <- SiteBySpeciesMatrix(Millstream, abundance='Log10Biomass', 
                               na.missing=TRUE)
    AssertEqual(c(67,2), dim(res))
    AssertEqual('Valvata piscinalis', rownames(res)[67])
    AssertEqual(c(c4=-5.0662, d4=-4.9292),res[1,])
    AssertEqual(c(c4=2.629, d4=NA), res[67,])
}
