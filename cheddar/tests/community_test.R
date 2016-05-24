TestCommunityBasic <- function()
{
    # Minimum community representation
    community <- Community(nodes=data.frame(node='S'), 
                           properties=list(title='test'))
    AssertEqual('test', CPS(community)$title)
    AssertEqual('S', NPS(community)$node)
    AssertEqual('node', NodePropertyNames(community))
    AssertEqual(TRUE, is.Community(community))

    # Modifications are illegal
    AssertRaises(community[1] <- 1)
    AssertRaises(community[[1]] <- 1)
    AssertRaises(community$silly <- 1)
    AssertRaises(names(community) <- letters[1:3])
    AssertRaises(length(community) <- 1)

    # Non-existent node properties
    AssertEqual(c(S=NA), NP(community, 'z'))
    AssertEqual(names(NP(community, 'z')), 'S')

    # No nodes
    AssertRaises(Community(nodes=data.frame(), properties=list(title='test')))

    # Nodes not a data.frame    
    AssertRaises(Community(node='S', properties=list(title='test')))

    # No node columns
    AssertRaises(Community(nodes=data.frame(a=LETTERS), properties=list(title='test')))

    # The node column appears twice
    AssertRaises(Community(nodes=data.frame(node=LETTERS, node=LETTERS, 
                                 check.names=FALSE), 
                 properties=list(title='test')))

    # Column names must not start with 'resource.' or 'consumer.'
    AssertRaises(Community(nodes=data.frame(node='S', resource.a=TRUE), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node='S', consumer.a=TRUE), 
                properties=list(title='test')))
 
    # Node names that are numbers are illegal
    AssertRaises(Community(nodes=data.frame(node=1:10), properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=c(1, LETTERS)), 
                properties=list(title='test')))

    # Ensure duplicate nodes are deteced
    AssertRaises(Community(nodes=data.frame(node=c('A', 'A')), 
                                 properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=c('A', 'A ')), 
                                 properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=c('A', 'a')), 
                                 properties=list(title='test')))

    # Empty node name
    AssertRaises(Community(nodes=data.frame(node=''), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=c('', LETTERS)), 
                properties=list(title='test')))

    # Some illegal node properties
    AssertRaises(Community(nodes=data.frame(node=LETTERS, title=LETTERS), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, M.units=LETTERS), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, N.units=LETTERS), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, resource=LETTERS), 
                properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, consumer=LETTERS), 
                properties=list(title='test')))

    # Factors should be converted to strings
    community <- Community(nodes=data.frame(node=LETTERS, np=letters,
                                            stringsAsFactors=TRUE), 
                           properties=list(title='test'))
    AssertEqual(sapply(NPS(community), class), 
                c(node='character', np='character'))
}

TestCommunityNPS <- function()
{
    # Ensure node properties are picked up
    community <- Community(nodes=data.frame(node=LETTERS, x=letters, M=1:26, 
                                            N=26:1), 
                           properties=list(title='test', M.units='kg', 
                                           N.units='m^-2'))
    AssertEqual(c('node', 'x', 'M', 'N'), NodePropertyNames(community))
    AssertEqual(LETTERS, unname(NP(community, 'node')))
    AssertEqual(letters, unname(NP(community, 'x')))
    AssertEqual(1:26, unname(NP(community, 'M')))
    AssertEqual(26:1, unname(NP(community, 'N')))
    AssertEqual(LETTERS, names(NP(community, 'node')))
    AssertEqual(LETTERS, names(NP(community, 'x')))
    AssertEqual(LETTERS, names(NP(community, 'M')))
    AssertEqual(LETTERS, names(NP(community, 'N')))

    # Non-existent node properties
    AssertEqual(rep(NA, 26), unname(NP(community, 'z')))
    AssertEqual(LETTERS, names(NP(community, 'z')))

    AssertRaises(NP(community, ''))
}

TestCommunityProperties <- function()
{
    # No properties at all
    AssertRaises(Community(nodes=data.frame(node='S')))

    # Properties does not have names
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list('test')))

    # Properties does not include title
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(xyz='test')))

    # Title in not in lowercase
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(Title='test')))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(TITLE='test')))

    # title duplicated
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', title='test')))

    # Properties is not a vector
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=cbind(title='test')))

    # One or more properties of length!=1
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title=paste('title', 1:10))))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', weasel=1:10)))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title=paste('title', 1:10), 
                                           weasel=1:10)))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', x=NULL)))

    # Some illegal names
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', node='a')))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', resource='a')))
    AssertRaises(Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', consumer='a')))

    # Ensure properties are picked up
    community <- Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', M.units='kg', pH=5.6, 
                                           lat=12, x=NA))
    AssertEqual('test', CP(community, 'title'))
    AssertEqual('kg', CP(community, 'M.units'))
    AssertEqual(5.6, CP(community, 'pH'))
    AssertEqual(12, CP(community, 'lat'))
    AssertEqual(NA, CP(community, 'x'))
    AssertEqual(list(title='test', pH=5.6, lat=12), 
                CPS(community, c('title', 'pH', 'lat')))
    AssertEqual(list(title='test', M.units='kg', pH=5.6, lat=12, x=NA), 
                CPS(community))
    AssertEqual(list(title='test', M.units='kg', pH=5.6, lat=12, x=NA), 
                CPS(community))

    # Computed properties
    AssertEqual(list(NumberOfTrophicLinks=0, pH=5.6), 
                CPS(community, c('NumberOfTrophicLinks', 'pH')))
    AssertEqual(list(NumberOfTrophicLinks=0), 
                CPS(c1, 'NumberOfTrophicLinks'))
    AssertEqual(list(NumberOfTrophicLinks=269), 
                CPS(TL84, 'NumberOfTrophicLinks'))

    # Computed property names
    AssertEqual(list(L=269), CPS(TL84, c(L='NumberOfTrophicLinks')))

    # Computed property names length > 1
    AssertEqual(list(invertebrate=0.00537614600000000005, 
                     producer=0.00743826670000000047, 
                     vert.ecto=0.00231559000000000018), 
                CPS(TL84, 'SumBiomassByClass'))

    # Properties as using c() should work just as well as a list()
    community <- Community(nodes=data.frame(node='S'), 
                           properties=list(title='test', a=1, b=TRUE))
    AssertEqual('test', CP(community, 'title'))
    AssertEqual(1, CP(community, 'a'))
    AssertEqual(TRUE, CP(community, 'b'))
}

TestCommunityTrophicLinks <- function()
{
    # No trophic Links
    community <- Community(nodes=data.frame(node=LETTERS), 
                           properties=list(title='test'))
    AssertEqual(c(26,26), dim(PredationMatrix(community)))
    AssertEqual(LETTERS, colnames(PredationMatrix(community)))
    AssertEqual(LETTERS, rownames(PredationMatrix(community)))
    AssertEqual(0, sum(PredationMatrix(community)))
    AssertEqual(0, NumberOfTrophicLinks(community))

    # Empty trophic links
    community <- Community(nodes=data.frame(node=LETTERS), 
                 trophic.links=matrix(0, nrow=0, ncol=2, 
                                      dimnames=list(NULL, c('resource', 
                                                            'consumer'))), 
                 properties=list(title='test'))
    AssertEqual(c(26,26), dim(PredationMatrix(community)))
    AssertEqual(LETTERS, colnames(PredationMatrix(community)))
    AssertEqual(LETTERS, rownames(PredationMatrix(community)))
    AssertEqual(0, sum(PredationMatrix(community)))
    AssertEqual(0, NumberOfTrophicLinks(community))

    # Make sure a minimum set of trophic.links is OK
    community <- Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A',consumer='A'),
                           properties=list(title='test'))
    AssertEqual(LETTERS, colnames(PredationMatrix(community)))
    AssertEqual(LETTERS, rownames(PredationMatrix(community)))
    AssertEqual(1, sum(PredationMatrix(community)))
    AssertEqual(1, NumberOfTrophicLinks(community))
    AssertEqual('integer', class(PredationMatrix(community)[1]))
    AssertEqual(0:1, sort(unique(as.vector(PredationMatrix(community)))))
    AssertEqual(1, PredationMatrix(community)[1,1])
    AssertEqual(rep(0, 26*26-1), PredationMatrix(community)[-1])
    AssertEqual(1, sum(PredationMatrix(community)))

    # resource and/or consumer not in nodes
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='a', consumer='A'),
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='a'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='a', consumer='a'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='', consumer='A'), 
                           properties=list(title='test')))

    # Check white space is removed
    community <- Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A ',consumer='A'),
                           properties=list(title='test'))
    AssertEqual(data.frame(resource='A',consumer='A'), TLPS(community))

    community <- Community(nodes=data.frame(node='A '), 
                           trophic.links=data.frame(resource='A',consumer='A'),
                           properties=list(title='test'))
    AssertEqual(data.frame(resource='A',consumer='A'), TLPS(community))
    community <- Community(nodes=data.frame(node=' A '), 
                           trophic.links=data.frame(resource=' A',consumer='A '),
                           properties=list(title='test'))
    AssertEqual(data.frame(resource='A',consumer='A'), TLPS(community))

    # Incorrect column names
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=cbind(x='A', t='A'), 
                           properties=list(title='test')))

    # Producers appear in 'consumer' column
    AssertRaises(Community(nodes=data.frame(node=LETTERS, 
                                   category=c(rep('producer', 25), 'consumer')), 
                           trophic.links=cbind(resource='A', consumer='W'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, 
                                   category=c(rep('producer', 25), 'consumer')),
                           trophic.links=cbind(resource=c('A', 'B'), 
                                               consumer=c('W', 'W')),
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS, 
                                   category=c(rep('producer', 25), 'consumer')),
                           trophic.links=cbind(resource=c('A', 'A'), 
                                               consumer=c('W', 'Z')),
                           properties=list(title='test')))

    # Some illegal trophic.link properties
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    node='A'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    title='A'),
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    M.units='A'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    N.units='A'), 
                           properties=list(title='test')))

    # trophic.link properties start with either 'resource.' or 'consumer.'
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    resource.a='A'), 
                           properties=list(title='test')))
    AssertRaises(Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource='A', consumer='A', 
                                                    consumer.='A'), 
                           properties=list(title='test')))

    # All cannibals - no other trophic links
    community <- Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource=LETTERS,
                                                    consumer=LETTERS),
                           properties=list(title='test'))
    pm <- PredationMatrix(community)
    AssertEqual(LETTERS, colnames(pm))
    AssertEqual(LETTERS, rownames(pm))
    AssertEqual(rep(1,26), unname(diag(pm)))
    AssertEqual(rep(0, 325), pm[upper.tri(pm, diag=FALSE)])
    AssertEqual(rep(0, 325), pm[lower.tri(pm, diag=FALSE)])

    # No cannibals, Z consumes everything
    community <- Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=data.frame(resource=LETTERS,
                                                    consumer='Z'),
                           properties=list(title='test'))
    pm <- PredationMatrix(community)
    AssertEqual(LETTERS, colnames(pm))
    AssertEqual(LETTERS, rownames(pm))
    AssertEqual(rep(1,26), unname(pm[,'Z']))
    AssertEqual(rep(0,650), as.vector(pm[,LETTERS[1:25]]))

    # A few specific links
    tl <- data.frame(resource=c('A', 'B', 'Z'), 
                     consumer=c('B', 'C', 'Z'), stringsAsFactors=FALSE)
    community <- Community(nodes=data.frame(node=LETTERS), 
                           trophic.links=tl,
                           properties=list(title='test'))
    pm <- PredationMatrix(community)
    AssertEqual(LETTERS, colnames(pm))
    AssertEqual(LETTERS, rownames(pm))
    invisible(apply(tl, 1, function(row)
    {
        res <- row['resource']
        con <- row['consumer']
        AssertEqual(1, pm[res, con])
    }))
    AssertEqual(0:1, sort(unique(as.vector(pm))))

    # Factors should be converted to strings
    community <- Community(nodes=data.frame(node=LETTERS),
                           trophic.links=data.frame(resource=LETTERS,
                                                    consumer=LETTERS,
                                                    tlp=LETTERS,
                                                    stringsAsFactors=TRUE),
                           properties=list(title='test'))
    AssertEqual(c(resource='character', consumer="character", tlp="character"), 
                sapply(TLPS(community), class))
}

TestNodePropertyNames <- function()
{
    for(community in list(c1,c2,c3,c4,c5))
    {
        print(community)
        AssertEqual('node', NodePropertyNames(community))
    }

    AssertEqual(c('node','M','N','order','family'), NodePropertyNames(c6))
}

TestNPS <- function()
{
    AssertEqual(c(S='S'), NP(c1, 'node'))
    AssertEqual(c(S='S'), NP(c2, 'node'))
    AssertEqual(c(R='R',C='C'), NP(c3, 'node'))
    AssertEqual(c(R='R',C='C',P='P'), NP(c4, 'node'))
    AssertEqual(c(R='R',C='C',O='O'), NP(c5, 'node'))
    AssertEqual(c(R='R',C='C',P='P'), NP(c6, 'node'))

    AssertEqual(data.frame(node='S', row.names=c('S')), NPS(c1))
    AssertEqual(data.frame(node='S', row.names=c('S')), NPS(c2))
    AssertEqual(data.frame(node=c('R','C'), row.names=c('R','C')), NPS(c3))
    AssertEqual(data.frame(node=c('R','C','P'), row.names=c('R','C','P')), 
                NPS(c4))
    AssertEqual(data.frame(node=c('R','C','O'), row.names=c('R','C','O')), 
                NPS(c5))
    AssertEqual(data.frame(node=c('R','C','P'), 
                           M=c(1.5,5,100), 
                           N=c(1000,10,5.5), 
                           order=c('Order 1','Order 2','Order 2'), 
                           family=paste('Family',1:3), 
                           row.names=c('R','C','P')), 
                NPS(c6))

    # Subsets of properties
    AssertEqual(data.frame(M=c(1.5,5,100), row.names=c('R','C','P')), 
                NPS(c6, 'M'))
    AssertEqual(data.frame(M=c(1.5,5,100), 
                         family=paste('Family',1:3), 
                         row.names=c('R','C','P')), 
                NPS(c6, c('M','family')))

    # Computed properties
    AssertEqual(data.frame(Biomass=c(1500,50,550), row.names=c('R','C','P')), 
                NPS(c6, 'Biomass'))
    AssertEqual(data.frame(M=c(1.5,5,100), 
                           PreyAveragedTrophicLevel=c(1,2,3), 
                           IsBasalNode=c(TRUE,FALSE,FALSE),
                           IsTopLevelNode=c(FALSE,FALSE,TRUE), 
                           row.names=c('R','C','P')), 
                NPS(c6, c('M', 'PreyAveragedTrophicLevel', 'IsBasalNode', 
                          'IsTopLevelNode')))

    # Computed property names
    AssertEqual(data.frame(x=c(1500,50,550), row.names=c('R','C','P')), 
                NPS(c6, c(x='Biomass')))

    # Computed property names with extra params
    AssertEqual(data.frame(TS=TrophicSpecies(TL84, include.isolated=FALSE)), 
                NPS(TL84, list(TS=list('TrophicSpecies', include.isolated=FALSE))))

    AssertEqual(data.frame(TS=TrophicSpecies(TL84, include.isolated=TRUE)), 
                NPS(TL84, list(TS=list('TrophicSpecies', include.isolated=TRUE))))

    # Function returns incorrect length
    AssertRaises(NPS(TL84, 'NumberOfTrophicLinks'))

    # Missing node properties
    AssertRaises(NP(TL84, ''))
    AssertEqual(rep(NA,56), unname(NP(TL84, 'x')))
}

TestTrophicLinkPropertyNames <- function()
{
    AssertEqual(NULL, TrophicLinkPropertyNames(c1))
    for(community in list(c2,c3,c4,c5))
    {
        print(community)
        AssertEqual(c('resource', 'consumer'), 
                    TrophicLinkPropertyNames(community))
    }

    AssertEqual(c('resource','consumer','link.evidence','link.strength'), 
                TrophicLinkPropertyNames(c6))
}

TestTLPS <- function()
{
    # First-class link properties
    AssertEqual(NULL, TLPS(c1))
    AssertEqual(data.frame(resource='S', consumer='S'), TLPS(c2))
    AssertEqual(data.frame(resource='R', consumer='C'), TLPS(c3))
    AssertEqual(data.frame(resource=c('R','C'), consumer=c('C','P')), TLPS(c4))
    AssertEqual(data.frame(resource=c('R','R','C'), 
                           consumer=c('C','O','O')), 
                TLPS(c5))
    AssertEqual(data.frame(resource=c('R','C'), consumer=c('C','P'), 
                           link.evidence=c('Inferred','Known'), 
                           link.strength=c(0.5, 0.2)), 
                TLPS(c6))

    AssertEqual(269, nrow(TLPS(TL84)))
    AssertEqual(264, nrow(TLPS(RemoveCannibalisticLinks(TL84))))

    AssertEqual(c('Inferred','Known'), TLP(c6, 'link.evidence'))
    AssertEqual(c(0.5, 0.2), TLP(c6, 'link.strength'))
    AssertEqual(c(NA, NA), TLP(c6, 'x'))

    # Node properties
    tlps <- TLPS(TL84, c('M','N','Biomass'))
    AssertEqual(colnames(tlps), c('resource', 'consumer', 
                                  'resource.M', 'resource.N', 
                                  'resource.Biomass', 
                                  'consumer.M', 'consumer.N', 
                                  'consumer.Biomass'))
    AssertEqual(tlps[,'resource.M'], unname(NP(TL84,'M')[tlps[,'resource']]))
    AssertEqual(tlps[,'consumer.M'], unname(NP(TL84,'M')[tlps[,'consumer']]))
    AssertEqual(tlps[,'resource.N'], unname(NP(TL84,'N')[tlps[,'resource']]))
    AssertEqual(tlps[,'consumer.N'], unname(NP(TL84,'N')[tlps[,'consumer']]))
    AssertEqual(tlps[,'resource.Biomass'], 
                unname(Biomass(TL84)[tlps[,'resource']]))
    AssertEqual(tlps[,'consumer.Biomass'], 
                unname(Biomass(TL84)[tlps[,'consumer']]))

    # Missing node properties should be NA
    tlps <- TLPS(TL84, 'x')
    AssertEqual(rep(NA, 269), tlps[,'resource.x'])
    AssertEqual(rep(NA, 269), tlps[,'consumer.x'])
}

TestCommunityPropertyNames <- function()
{
    AssertEqual('title', CommunityPropertyNames(c1))
    AssertEqual('title', CommunityPropertyNames(c2))
    AssertEqual('title', CommunityPropertyNames(c3))
    AssertEqual('title', CommunityPropertyNames(c4))
    AssertEqual('title', CommunityPropertyNames(c5))
    AssertEqual(c('title','M.units','N.units'), CommunityPropertyNames(c6))
}

TestPersistCommunity <- function()
{
    # Tests LoadCommunity() and SaveCommunity()
    path <- tempfile()
    on.exit(unlink(path, recursive=TRUE))
    for(community in list(c1,c2,c3,c4,c5,c6,TL84,TL86,YthanEstuary,
                          SkipwithPond,BroadstoneStream,Benguela))
    {
        print(community)
        SaveCommunity(community, path)
        AssertEqual(community, LoadCommunity(path))
        unlink(path, recursive=TRUE)

        SaveCommunity(community, path, fn='write.table', sep='\t')
        AssertEqual(community, LoadCommunity(path, fn='read.table', sep='\t'))
        unlink(path, recursive=TRUE)
    }
}

TestResolveToNodeIndices <- function()
{
    ResolveToNodeIndices <- cheddar:::.ResolveToNodeIndices
    AssertEqual(c(S=1), ResolveToNodeIndices(c1, 'S'))
    AssertEqual(c(S=1), ResolveToNodeIndices(c1, 1))
    AssertEqual(c(S=1), ResolveToNodeIndices(c1, TRUE))
    AssertEqual(0, length(ResolveToNodeIndices(c1, FALSE)))

    AssertEqual(c(R=1), ResolveToNodeIndices(c4, 'R'))
    AssertEqual(c(R=1,P=3), ResolveToNodeIndices(c4, c('R','P')))
    AssertEqual(c(R=1), ResolveToNodeIndices(c4, 1))
    AssertEqual(c(R=1,C=2), ResolveToNodeIndices(c4, c(1,2)))
    AssertEqual(c(C=2,P=3), ResolveToNodeIndices(c4, c(FALSE,TRUE,TRUE)))
    AssertEqual(0, length(ResolveToNodeIndices(c4, c(FALSE,FALSE,FALSE))))

    AssertEqual(NULL, ResolveToNodeIndices(c4, NULL))

    AssertRaises(ResolveToNodeIndices(c4, ''))
    AssertRaises(ResolveToNodeIndices(c4, 'x'))
    AssertRaises(ResolveToNodeIndices(c4, 0))
    AssertRaises(ResolveToNodeIndices(c4, 4))
    AssertRaises(ResolveToNodeIndices(c4, TRUE))
    AssertRaises(ResolveToNodeIndices(c4, rep(TRUE,4)))
}

TestNodeNameIndices <- function()
{
    AssertEqual(c(S=1), NodeNameIndices(c1, 'S'))
    AssertEqual(c(S=1), NodeNameIndices(c2, 'S'))
    AssertEqual(c(C=2,R=1), NodeNameIndices(c3, c('C','R')))
    AssertEqual(c(R=1,C=2), NodeNameIndices(c3, c('R','C')))
    AssertEqual(c(C=2,C=2), NodeNameIndices(c3, c('C','C')))
    AssertEqual(c(R=1,C=2,P=3), NodeNameIndices(c4, c('R','C','P')))
    AssertEqual(c(P=3,C=2,R=1), NodeNameIndices(c4, c('P','C','R')))
    AssertEqual(c(C=2), NodeNameIndices(c4, 'C'))

    AssertRaises(NodeNameIndices(c1, ''))
    AssertRaises(NodeNameIndices(c1, 'x'))
}

TestNumberOfNodes <- function()
{
    AssertEqual(1, NumberOfNodes(c1))
    AssertEqual(1, NumberOfNodes(c2))
    AssertEqual(2, NumberOfNodes(c3))
    AssertEqual(3, NumberOfNodes(c4))
    AssertEqual(3, NumberOfNodes(c5))
    AssertEqual(3, NumberOfNodes(c6))
    AssertEqual(56, NumberOfNodes(TL84))
    AssertEqual(57, NumberOfNodes(TL86))
    AssertEqual(37, NumberOfNodes(SkipwithPond))
}

TestRemoveNodes <- function()
{
    AssertRaises(RemoveNodes(c1, 0))  # Illegal node index
    AssertRaises(RemoveNodes(c1, 1))  # Would result in empty community
    AssertRaises(RemoveNodes(c1, 2))  # Illegal node index
    AssertRaises(RemoveNodes(c2, 1))  # Would result in empty community

    AssertEqual(c1, RemoveNodes(c1, NULL))
    AssertEqual(c1, RemoveNodes(c1, vector(mode='character')))

    AssertEqual(RemoveNodes(c3, 1, title='c3'), 
                RemoveNodes(c3, 'R', title='c3'))
    AssertEqual(c(C='C'), NP(RemoveNodes(c3, 1),'node'))
    AssertEqual(NULL, TLPS(RemoveNodes(c3, 1)))
    AssertEqual(c(R='R'), NP(RemoveNodes(c3, 2),'node'))
    AssertEqual(NULL, TLPS(RemoveNodes(c3, 2)))

    AssertEqual(c(C='C',P='P'), NP(RemoveNodes(c4, 1),'node'))
    AssertEqual(c(R='R',P='P'), NP(RemoveNodes(c4, 2),'node'))
    AssertEqual(c(R='R',C='C'), NP(RemoveNodes(c4, 3),'node'))
    AssertEqual(c(R='R'), NP(RemoveNodes(c4, 2:3),'node'))
    AssertRaises(RemoveNodes(c4, 1:3))

    AssertEqual(TL84, RemoveNodes(TL84, NULL))
    AssertEqual(TL84, RemoveNodes(TL84, vector(mode='character')))

    TL84r <- RemoveNodes(TL84, 56)
    AssertEqual(NPS(TL84)[-56,], NPS(TL84r))
    TL84r <- RemoveNodes(TL84, 'Umbra limi')
    AssertEqual(NPS(TL84)[-56,], NPS(TL84r))

    TL84r <- RemoveNodes(TL84, c(1,56))
    AssertEqual(NPS(TL84)[-c(1,56),],NPS(TL84r))

    TL84r <- RemoveNodes(TL84, c('Nostoc sp.', 'Umbra limi'))
    AssertEqual(NPS(TL84)[-c(1,56),],NPS(TL84r))

    # Test secondary and cascade
    AssertRaises(RemoveNodes(c3, 'R', method='secondary'))
    AssertRaises(RemoveNodes(c3, 'R', method='cascade'))
    AssertEqual(c(P='P'), NP(RemoveNodes(c4, 'R', method='secondary'), 'node'))
    AssertRaises(RemoveNodes(c4, 'R', method='cascade'))
    AssertEqual(c(O='O'), NP(RemoveNodes(c5, 'R', method='secondary'), 'node'))
    AssertRaises(RemoveNodes(c5, 'R', method='cascade'))
    AssertEqual(c(R='R', O='O'), NP(RemoveNodes(c5, 'C', method='secondary'), 'node'))
    AssertEqual(c(R='R', C='C'), NP(RemoveNodes(c5, 'O', method='secondary'), 'node'))
    AssertEqual(c(C='C', D='D', E='E'), 
                NP(RemoveNodes(c7, 'A', method='secondary'), 'node'))
    AssertEqual(c(D='D', E='E'), 
                NP(RemoveNodes(c7, 'A', method='cascade'), 'node'))
    AssertEqual(c(A='A', D='D', E='E'), 
                NP(RemoveNodes(c7, 'B', method='secondary'), 'node'))
    AssertEqual(c(A='A', B='B', D='D', E='E'), 
                NP(RemoveNodes(c7, 'C', method='secondary'), 'node'))
    AssertEqual(c(A='A', B='B', C='C', E='E'), 
                NP(RemoveNodes(c7, 'D', method='secondary'), 'node'))
    AssertEqual(c(A='A', B='B', C='C', D='D'), 
                NP(RemoveNodes(c7, 'E', method='secondary'), 'node'))

    # 56 - 25 = 31 nodes remain
    AssertEqual(31, NumberOfNodes(RemoveNodes(TL84, BasalNodes(TL84))))
    AssertEqual(14, NumberOfNodes(RemoveNodes(TL84, BasalNodes(TL84), method='secondary')))
    cascaded <- NumberOfNodes(RemoveNodes(TL84, BasalNodes(TL84), method='cascade'))
    AssertEqual(6, cascaded)
    isolated <- c('Asterionella formosa','Chrysosphaerella longispina',
                  'Diceras sp.', 'Rhizosolenia sp.', 'Spinocosmarium sp.',
                  'Staurastrum sp.')
    AssertEqual(isolated, IsolatedNodes(TL84))
}

TestRemoveIsolatedNodes <- function()
{
    AssertRaises(RemoveIsolatedNodes(c1))
    AssertRaises(RemoveIsolatedNodes(c2))
    AssertEqual(c3, RemoveIsolatedNodes(c3, title=c3$title))
    AssertEqual(c4, RemoveIsolatedNodes(c4, title=c4$title))
    AssertEqual(c5, RemoveIsolatedNodes(c5, title=c5$title))
    AssertEqual(c6, RemoveIsolatedNodes(c6, title=c6$title))

    TL84.no.iso <- RemoveIsolatedNodes(TL84, title=TL84$title)
    AssertEqual(50, NumberOfNodes(TL84.no.iso))

    isolated <- c('Asterionella formosa','Chrysosphaerella longispina',
                  'Diceras sp.', 'Rhizosolenia sp.', 'Spinocosmarium sp.',
                  'Staurastrum sp.')
    AssertEqual(isolated, IsolatedNodes(TL84))

    AssertEqual(NPS(TL84)[!NP(TL84,'node') %in% isolated,], NPS(TL84.no.iso))
    AssertEqual(TLPS(TL84), TLPS(TL84.no.iso))
}

TestLumpNodes <- function()
{
    # Don't lump any nodes. The returned community should be identical to the 
    # original.
    lump <- NP(TL84, 'node')
    TL84.lumped <- LumpNodes(TL84, lump, title=CP(TL84, 'title'))
    AssertEqual(TL84, TL84.lumped)

    # Lump all nodes into a single node
    TL84.lumped <- LumpNodes(TL84, paste('S', rep(1, 56)))
    AssertEqual(1, NumberOfNodes(TL84.lumped))
    AssertEqual(weighted.mean(NP(TL84, 'M'), NP(TL84, 'N')), 
                unname(NP(TL84.lumped, 'M')))
    AssertEqual(mean(NP(TL84, 'N')), unname(NP(TL84.lumped, 'N')))
    AssertEqual("producer,invertebrate,vert.ecto", 
                unname(NP(TL84.lumped, 'category')))

    # Lump all nodes into a single node without weighting by N
    TL84.lumped <- LumpNodes(TL84, paste('S', rep(1, 56)), weight.by=NULL)
    AssertEqual(1, NumberOfNodes(TL84.lumped))
    AssertEqual(mean(NP(TL84, 'M')), unname(NP(TL84.lumped, 'M')))
    AssertEqual(mean(NP(TL84, 'N')), unname(NP(TL84.lumped, 'N')))
    AssertEqual("producer,invertebrate,vert.ecto", 
                unname(NP(TL84.lumped, 'category')))

    # Lump some specific nodes
    lump <- NP(TL84, 'node')
    lump[c(1,3)] <- 'Lump 1 and 3'
    lump[c(2,4)] <- 'Lump 2 and 4'
    TL84.lumped <- LumpNodes(TL84, lump)
    AssertEqual(54, NumberOfNodes(TL84.lumped))
    AssertEqual(weighted.mean(NP(TL84, 'M')[c(1,3)], NP(TL84, 'N')[c(1,3)]), 
                unname(NP(TL84.lumped, 'M')['Lump 1 and 3']))
    AssertEqual(mean(NP(TL84, 'N')[c(1,3)]), 
                unname(NP(TL84.lumped, 'N')['Lump 1 and 3']))
    AssertEqual('producer', unname(NP(TL84.lumped, 'category')['Lump 1 and 3']))
    AssertEqual('Bacteria,Chromista', 
                unname(NP(TL84.lumped, 'kingdom')['Lump 1 and 3']))
    AssertEqual(weighted.mean(NP(TL84, 'M')[c(2,4)], NP(TL84, 'N')[c(2,4)]), 
                unname(NP(TL84.lumped, 'M')['Lump 2 and 4']))
    AssertEqual(mean(NP(TL84, 'N')[c(2,4)]), 
                unname(NP(TL84.lumped, 'N')['Lump 2 and 4']))
    AssertEqual('producer', unname(NP(TL84.lumped, 'category')['Lump 2 and 4']))
    AssertEqual('Plantae,Chromista', 
                unname(NP(TL84.lumped, 'kingdom')['Lump 2 and 4']))
    AssertEqual(NPS(TL84)[5:56,], NPS(TL84.lumped)[3:54,])


    # Lump isolated nodes
    AssertEqual(56, NumberOfNodes(TL84))
    AssertEqual(6, length(IsolatedNodes(TL84)))
    lump <- NP(TL84, 'node')
    lump[IsolatedNodes(TL84)] <- 'Isolated'
    TL84.lumped <- LumpNodes(TL84, lump)
    AssertEqual(51, NumberOfNodes(TL84.lumped))
    AssertEqual(weighted.mean(NP(TL84, 'M')[IsolatedNodes(TL84)], 
                              NP(TL84, 'N')[IsolatedNodes(TL84)]), 
                unname(NP(TL84.lumped, 'M')['Isolated']))
    AssertEqual(mean(NP(TL84, 'N')[IsolatedNodes(TL84)]), 
                unname(NP(TL84.lumped, 'N')['Isolated']))
    AssertEqual("producer", unname(NP(TL84.lumped, 'category')['Isolated']))

    # Lump isolated nodes without weighting by N
    AssertEqual(56, NumberOfNodes(TL84))
    AssertEqual(6, length(IsolatedNodes(TL84)))
    lump <- NP(TL84, 'node')
    lump[IsolatedNodes(TL84)] <- 'Isolated'
    TL84.lumped <- LumpNodes(TL84, lump, weight.by=NULL)
    AssertEqual(51, NumberOfNodes(TL84.lumped))
    AssertEqual(mean(NP(TL84, 'M')[IsolatedNodes(TL84)]), 
                unname(NP(TL84.lumped, 'M')['Isolated']))
    AssertEqual(mean(NP(TL84, 'N')[IsolatedNodes(TL84)]), 
                unname(NP(TL84.lumped, 'N')['Isolated']))
    AssertEqual("producer", unname(NP(TL84.lumped, 'category')['Isolated']))


    # Lump Ythan Estuary nodes
    lump <- NP(YthanEstuary, 'node')

    # European flounder:
    # "Platichthys flesus" and "Platichthys flesus (juvenile)"
    # Lump these in to one node
    lump["Platichthys flesus (juvenile)"==lump] <- "Platichthys flesus"

    # Common eider:
    # "Somateria mollissima" and "Somateria mollissima (juvenile)"
    # Lump these in to one node
    lump["Somateria mollissima (juvenile)"==lump] <- "Somateria mollissima"

    # Weight by N
    YthanEstuary.lumped <- LumpNodes(YthanEstuary, lump)
    AssertEqual(unname(NP(YthanEstuary.lumped, 'M')["Somateria mollissima"]), 
                1637.01774691358014024445)
    AssertEqual(unname(NP(YthanEstuary.lumped, 'N')["Somateria mollissima"]), 
                2592)

    # No weighting
    YthanEstuary.lumped2 <- LumpNodes(YthanEstuary, lump, weight.by=NULL)
    AssertEqual(unname(NP(YthanEstuary.lumped2, 'M')["Somateria mollissima"]), 
                1425)
    AssertEqual(unname(NP(YthanEstuary.lumped2, 'N')["Somateria mollissima"]), 
                2592)
}

TestLumpTrophicSpecies <- function()
{
    AssertEqual(1, NumberOfNodes(c1))
    AssertEqual(1, NumberOfNodes(LumpTrophicSpecies(c1)))
    AssertEqual(1, NumberOfNodes(c2))
    AssertEqual(1, NumberOfNodes(LumpTrophicSpecies(c2)))
    AssertEqual(2, NumberOfNodes(c3))
    AssertEqual(2, NumberOfNodes(LumpTrophicSpecies(c3)))
    AssertEqual(3, NumberOfNodes(c4))
    AssertEqual(3, NumberOfNodes(LumpTrophicSpecies(c4)))
    AssertEqual(3, NumberOfNodes(c5))
    AssertEqual(3, NumberOfNodes(LumpTrophicSpecies(c5)))
    AssertEqual(3, NumberOfNodes(c6))
    AssertEqual(3, NumberOfNodes(LumpTrophicSpecies(c6)))
    AssertEqual(c('Inferred','Known'), 
                TLPS(LumpTrophicSpecies(c6))[,'link.evidence'])
    AssertEqual(c(0.5,0.2), 
                TLPS(LumpTrophicSpecies(c6))[,'link.strength'])

    AssertEqual(56, NumberOfNodes(TL84))

    # Exclude isolated species.
    lumped <- LumpTrophicSpecies(TL84, include.isolated=FALSE, weight.by=NULL)
    AssertEqual(21, NumberOfNodes(lumped))

    # From Jonsson et al 2005 AER. Isolated species assigned NA.
    trophic.species <- c(1,2,NA,3,4,3,5,NA,6,1,1,NA,4,7,4,8,6,7,2,3,6,7,
                         4,6,3,NA,3,NA,NA,6,3,9,9,10,11,12,13,14,15,11,
                         15,16,15,9,15,17,18,15,15,15,15,12,19,20,20,21)

    # Are my trophic species the same as those calculated by cheddar?
    AssertEqual(trophic.species, 
                unname(TrophicSpecies(TL84, include.isolated=FALSE)))

    # Check M and N have been averaged correctly
    for(ts in unique(trophic.species[!is.na(trophic.species)]))
    {
        AssertEqual(unname(NP(lumped,'M')[ts]), 
                    mean(NP(TL84,'M')[which(ts==trophic.species)]))
        AssertEqual(unname(NP(lumped,'N')[ts]), 
                    mean(NP(TL84,'N')[which(ts==trophic.species)]))
    }

    # Include isolated species.
    lumped <- LumpTrophicSpecies(TL84, include.isolated=TRUE, weight.by=NULL)
    AssertEqual(22, NumberOfNodes(lumped))

    trophic.species <- c(1,2,3,4,5,4,6,3,7,1,1,3,5,8,5,9,7,8,2,4,7,8,5,7,
                         4,3,4,3,3,7,4,10,10,11,12,13,14,15,16,12,16,17,16,
                         10,16,18,19,16,16,16,16,13,20,21,21,22)

    # Are my trophic species the same as those calculated by cheddar?
    AssertEqual(trophic.species, 
                unname(TrophicSpecies(TL84, include.isolated=TRUE)))

    # Check M and N have been averaged correctly
    for(ts in unique(trophic.species[!is.na(trophic.species)]))
    {
        AssertEqual(unname(NP(lumped,'M')[ts]), 
                    mean(NP(TL84,'M')[which(ts==trophic.species)]))
        AssertEqual(unname(NP(lumped,'N')[ts]), 
                    mean(NP(TL84,'N')[which(ts==trophic.species)]))
    }
}

TestOrderCommunity <- function()
{
    # Order unchanged
    AssertEqual(c1, OrderCommunity(c1, 'node', title='c1'))
    AssertEqual(c2, OrderCommunity(c2, 'node', title='c2'))
    AssertEqual(c3, OrderCommunity(c3, new.order=1:2, title='c3'))
    AssertEqual(c4, OrderCommunity(c4, new.order=1:3, title='c4'))
    AssertEqual(c5, OrderCommunity(c5, new.order=1:3, title='c5'))
    AssertEqual(c6, OrderCommunity(c6, new.order=1:3, title='c6'))

    # Reverse order
    c6r <- OrderCommunity(c6, new.order=3:1)
    AssertEqual(NPS(c6), NPS(c6r)[3:1,])

    target1 <- TLPS(c6)
    target1 <- target1[order(target1$resource, target1$consumer),]
    target2 <- TLPS(c6r)
    target2 <- target2[order(target2$resource, target2$consumer),]
    AssertEqual(target1, target2)

    # Order nodes alphabetically
    SkipwithPondr <- OrderCommunity(SkipwithPond, 'node')
    target1 <- NPS(SkipwithPond)[order(NP(SkipwithPond,'node')),]
    target2 <- NPS(SkipwithPondr)
    AssertEqual(target1, target2)
    target1 <- TLPS(SkipwithPond)
    target1 <- target1[order(target1$resource, target1$consumer),]
    target2 <- TLPS(SkipwithPondr)
    target2 <- target2[order(target2$resource, target2$consumer),]
    AssertEqual(target1, target2)

    # Order by property
    AssertEqual(c("Chromulina sp.","Dactylococcopsis fascicularis", 
                  "Diceras sp.", "Trachelomonas sp.", "Cryptomonas sp. 1", 
                  "Closteriopsis longissimus", "Chroococcus dispersus",
                  "Selenastrum minutum", "Unclassified flagellates",
                  "Dictyosphaerium pulchellum", "Dinobryon sociale",
                  "Rhizosolenia sp.","Nostoc sp.","Mallomonas sp. 1", 
                  "Asterionella formosa","Mallomonas sp. 2",
                  "Cryptomonas sp. 2","Arthrodesmus sp.", 
                  "Dinobryon cylindricum","Peridinium pulsillum", 
                  "Dinobryon bavaricum","Spinocosmarium sp.", 
                  "Staurastrum sp.",  "Glenodinium quadridens", 
                  "Keratella cochlearis","Keratella testudo",
                  "Dinobryon sertularia","Microcystis aeruginosa", 
                  "Kellicottia sp.",  "Conochilus (solitary)", 
                  "Peridinium wisconsinense","Peridinium cinctum", 
                  "Peridinium limbatum","Synedra sp.","Ploesoma sp.", 
                  "Gastropus stylifer","Ascomorpha eucadis", 
                  "Conochiloides dossuarius","Filinia longispina", 
                  "Trichocerca multicrinis", "Trichocerca cylindrica",
                  "Polyarthra vulgaris","Chrysosphaerella longispina", 
                  "Synchaeta sp.","Bosmina longirostris",
                  "Diaphanosoma leuchtenbergianum", "Tropocyclops prasinus",
                  "Leptodiaptomus siciloides","Cyclops varians rubellus",
                  "Orthocyclops modestus","Daphnia pulex",
                  "Holopedium gibberum", "Chaoborus punctipennis",
                  "Phoxinus eos","Phoxinus neogaeus","Umbra limi"), 
                 unname(NP(OrderCommunity(TL84, 'M'), 'node')))
    AssertEqual(c('Umbra limi','Phoxinus neogaeus','Phoxinus eos',
                  'Holopedium gibberum','Daphnia pulex','Filinia longispina',
                  'Leptodiaptomus siciloides','Keratella testudo',
                  'Cyclops varians rubellus','Chaoborus punctipennis',
                  'Ascomorpha eucadis','Diaphanosoma leuchtenbergianum',
                  'Orthocyclops modestus','Synchaeta sp.',
                  'Gastropus stylifer','Conochilus (solitary)',
                  'Trichocerca multicrinis','Tropocyclops prasinus',
                  'Ploesoma sp.','Trichocerca cylindrica',
                  'Bosmina longirostris','Conochiloides dossuarius',
                  'Kellicottia sp.','Polyarthra vulgaris',
                  'Keratella cochlearis','Synedra sp.','Nostoc sp.',
                  'Dinobryon sertularia','Spinocosmarium sp.',
                  'Dinobryon cylindricum','Chrysosphaerella longispina',
                  'Asterionella formosa','Peridinium cinctum',
                  'Staurastrum sp.','Dictyosphaerium pulchellum',
                  'Microcystis aeruginosa','Diceras sp.',
                  'Peridinium wisconsinense','Peridinium limbatum',
                  'Mallomonas sp. 1','Chroococcus dispersus',
                  'Mallomonas sp. 2','Cryptomonas sp. 2','Dinobryon sociale',
                  'Dinobryon bavaricum','Dactylococcopsis fascicularis',
                  'Arthrodesmus sp.','Rhizosolenia sp.','Cryptomonas sp. 1',
                  'Glenodinium quadridens','Closteriopsis longissimus',
                  'Peridinium pulsillum','Chromulina sp.',
                  'Selenastrum minutum','Trachelomonas sp.',
                  'Unclassified flagellates'), 
                unname(NP(OrderCommunity(TL84, 'N'), 'node')))

    # Reorder by function
    AssertEqual(c('Asterionella formosa','Chrysosphaerella longispina',
                  'Diceras sp.','Rhizosolenia sp.','Spinocosmarium sp.',
                  'Staurastrum sp.','Dinobryon bavaricum',
                  'Microcystis aeruginosa','Peridinium limbatum',
                  'Peridinium wisconsinense','Synedra sp.',
                  'Dinobryon sertularia','Mallomonas sp. 1',
                  'Peridinium cinctum','Arthrodesmus sp.',
                  'Closteriopsis longissimus','Mallomonas sp. 2','Nostoc sp.',
                  'Dinobryon cylindricum','Dactylococcopsis fascicularis',
                  'Cryptomonas sp. 2','Dictyosphaerium pulchellum',
                  'Dinobryon sociale','Peridinium pulsillum',
                  'Glenodinium quadridens','Filinia longispina',
                  'Gastropus stylifer','Kellicottia sp.','Keratella testudo',
                  'Ploesoma sp.','Polyarthra vulgaris',
                  'Trichocerca multicrinis','Trichocerca cylindrica',
                  'Phoxinus eos','Phoxinus neogaeus','Ascomorpha eucadis',
                  'Synchaeta sp.','Conochilus (solitary)',
                  'Conochiloides dossuarius','Keratella cochlearis',
                  'Umbra limi','Bosmina longirostris','Cryptomonas sp. 1',
                  'Chroococcus dispersus','Unclassified flagellates',
                  'Chromulina sp.','Selenastrum minutum','Trachelomonas sp.',
                  'Diaphanosoma leuchtenbergianum','Holopedium gibberum',
                  'Orthocyclops modestus','Cyclops varians rubellus',
                  'Leptodiaptomus siciloides','Tropocyclops prasinus',
                  'Chaoborus punctipennis','Daphnia pulex'), 
                unname(NP(OrderCommunity(TL84, 'Degree'), 'node')))
    AssertEqual(c("Daphnia pulex", "Chaoborus punctipennis", 
                  "Cyclops varians rubellus", "Leptodiaptomus siciloides", 
                  "Tropocyclops prasinus", "Holopedium gibberum", 
                  "Orthocyclops modestus", "Cryptomonas sp. 1",           
                  "Chroococcus dispersus", "Unclassified flagellates", 
                  "Chromulina sp.", "Selenastrum minutum", 
                  "Trachelomonas sp.", "Diaphanosoma leuchtenbergianum", 
                  "Bosmina longirostris", "Umbra limi",       
                  "Ascomorpha eucadis", "Synchaeta sp.", 
                  "Conochilus (solitary)", "Conochiloides dossuarius", 
                  "Keratella cochlearis", "Filinia longispina", 
                  "Gastropus stylifer", "Kellicottia sp.", 
                  "Keratella testudo", "Ploesoma sp.", 
                  "Polyarthra vulgaris", "Trichocerca multicrinis", 
                  "Trichocerca cylindrica", "Phoxinus eos", 
                  "Phoxinus neogaeus", "Glenodinium quadridens", 
                  "Cryptomonas sp. 2", "Dictyosphaerium pulchellum", 
                  "Dinobryon sociale", "Peridinium pulsillum", 
                  "Nostoc sp.", "Dinobryon cylindricum", 
                  "Dactylococcopsis fascicularis", "Arthrodesmus sp.", 
                  "Closteriopsis longissimus", "Mallomonas sp. 2", 
                  "Dinobryon sertularia", "Mallomonas sp. 1", 
                  "Peridinium cinctum", "Dinobryon bavaricum", 
                  "Microcystis aeruginosa", "Peridinium limbatum", 
                  "Peridinium wisconsinense", "Synedra sp.", 
                  "Asterionella formosa", "Chrysosphaerella longispina", 
                  "Diceras sp.", "Rhizosolenia sp.", "Spinocosmarium sp.", 
                  "Staurastrum sp."), 
                unname(NP(OrderCommunity(TL84, 'Degree', decreasing=TRUE), 'node')))
}

TestNumberOfNodesByClass <- function()
{
    AssertEqual(c(invertebrate=22, producer=31, vert.ecto=3), 
                NumberOfNodesByClass(TL84))
    AssertEqual(c('<unnamed>'=2, invertebrate=34, producer=1), 
                NumberOfNodesByClass(BroadstoneStream))
}

TestFractionOfNodesByClass <- function()
{
    AssertEqual(c(invertebrate=0.39285714285714284921, 
                  producer=0.55357142857142860315, 
                  vert.ecto=0.05357142857142856845), 
                FractionOfNodesByClass(TL84))
}
