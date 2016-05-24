# Functions requiring trophic links
TestPredationMatrix <- function()
{
    AssertEqual(matrix(0, ncol=1, nrow=1, dimnames=list('S', 'S')), 
               PredationMatrix(c1))
    AssertEqual(matrix(1, ncol=1, nrow=1, dimnames=list('S', 'S')), 
                        PredationMatrix(c2))
    AssertEqual(matrix(c(0,0,1,0), ncol=2, nrow=2, 
                               dimnames=list(c('R','C'), c('R','C'))), 
                        PredationMatrix(c3))
    AssertEqual(matrix(c(0,0,0, 1,0,0, 0,1,0), ncol=3, nrow=3, 
                               dimnames=list(c('R','C','P'), c('R','C','P'))), 
                        PredationMatrix(c4))
    AssertEqual(matrix(c(0,0,0, 1,0,0, 1,1,0), ncol=3, nrow=3, 
                               dimnames=list(c('R','C','O'), c('R','C','O'))), 
                        PredationMatrix(c5))
    AssertEqual(matrix(c(0,0,0, 1,0,0, 0,1,0), ncol=3, nrow=3, 
                               dimnames=list(c('R','C','P'), c('R','C','P'))), 
                        PredationMatrix(c6))
    AssertEqual(unname(NP(TL84, 'node')), rownames(PredationMatrix(TL84)))
    AssertEqual(unname(NP(TL84, 'node')), colnames(PredationMatrix(TL84)))
    AssertEqual(269, sum(PredationMatrix(TL84)))
}

TestInDegree <- function()
{
    AssertEqual(c(S=0), InDegree(c1))
    AssertEqual(c(S=1), InDegree(c2))
    AssertEqual(c(R=0,C=1), InDegree(c3))
    AssertEqual(c(R=0,C=1,P=1), InDegree(c4))
    AssertEqual(c(R=0,C=1,O=2), InDegree(c5))
    AssertEqual(c(R=0,C=1,P=1), InDegree(c6))
    AssertEqual(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,6,6,11,7,16,14,28,6,7,6,19,6,6,6,16,17,6,6,6,6,16,
                  22,9,9,12), unname(InDegree(TL84)))
}

TestNormalisedTrophicGenerality <- function()
{
    AssertEqual(c(S=NaN), NormalisedTrophicGenerality(c1))
    AssertEqual(c(S=1), NormalisedTrophicGenerality(c2))
    AssertEqual(c(R=0,C=2), NormalisedTrophicGenerality(c3))
    AssertEqual(c(R=0,C=1.5,P=1.5), NormalisedTrophicGenerality(c4))
    AssertEqual(c(R=0,C=1,O=2), NormalisedTrophicGenerality(c5))
    AssertEqual(c(R=0,C=1.5,P=1.5), NormalisedTrophicGenerality(c6))
    AssertEqual(56, sum(NormalisedTrophicGenerality(TL84)))
    AssertEqual( 1, mean(NormalisedTrophicGenerality(TL84)))
}

TestOutDegree <- function()
{
    AssertEqual(c(S=0), OutDegree(c1))
    AssertEqual(c(S=1), OutDegree(c2))
    AssertEqual(c(R=1,C=0), OutDegree(c3))
    AssertEqual(c(R=1,C=1,P=0), OutDegree(c4))
    AssertEqual(c(R=2,C=1,O=0), OutDegree(c5))
    AssertEqual(c(R=1,C=1,P=0), OutDegree(c6))
    AssertEqual(c(4,3,0,18,5,18,3,0,1,4,4,0,5,2,5,6,1,2,3,18,1,2,5,1,18,
                  0,18,0,0,1,18,5,5,4,4,7,4,4,4,4,4,3,4,5,4,7,5,4,4,4,4,
                  7,4,1,1,1), unname(OutDegree(TL84)))
}

TestNormalisedTrophicVulnerability <- function()
{
    AssertEqual(c(S=NaN), NormalisedTrophicVulnerability(c1))
    AssertEqual(c(S=1), NormalisedTrophicVulnerability(c2))
    AssertEqual(c(R=2,C=0), NormalisedTrophicVulnerability(c3))
    AssertEqual(c(R=1.5,C=1.5,P=0), NormalisedTrophicVulnerability(c4))
    AssertEqual(c(R=2,C=1,O=0), NormalisedTrophicVulnerability(c5))
    AssertEqual(c(R=1.5,C=1.5,P=0), NormalisedTrophicVulnerability(c6))
    AssertEqual(56, sum(NormalisedTrophicVulnerability(TL84)))
    AssertEqual( 1, mean(NormalisedTrophicVulnerability(TL84)))
}

TestBasalNodes <- function()
{
    AssertEqual(c(S=FALSE), IsBasalNode(c1))
    AssertEqual(0, length(BasalNodes(c1)))
    AssertEqual(c(S=FALSE), IsBasalNode(c2))
    AssertEqual(0, length(BasalNodes(c2)))
    AssertEqual(c(R=TRUE, C=FALSE), IsBasalNode(c3))
    AssertEqual('R', BasalNodes(c3))
    AssertEqual(c(R=TRUE, C=FALSE, P=FALSE), IsBasalNode(c4))
    AssertEqual('R', BasalNodes(c4))
    AssertEqual(c(R=TRUE, C=FALSE, O=FALSE), IsBasalNode(c5))
    AssertEqual('R', BasalNodes(c5))
    AssertEqual(c(R=TRUE, C=FALSE, P=FALSE), IsBasalNode(c6))
    AssertEqual('R', BasalNodes(c6))
    AssertEqual(c(A=TRUE, B=FALSE, C=FALSE, D=FALSE, E=FALSE), IsBasalNode(c7))
    AssertEqual('A', BasalNodes(c7))
    AssertEqual(c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,
                    TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
                    TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,
                    TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE), 
                    unname(IsBasalNode(TL84)))
    AssertEqual(c('Nostoc sp.','Arthrodesmus sp.','Cryptomonas sp. 1',
                  'Cryptomonas sp. 2','Chroococcus dispersus',
                  'Closteriopsis longissimus','Dinobryon bavaricum',
                  'Dinobryon cylindricum','Dactylococcopsis fascicularis',
                  'Dictyosphaerium pulchellum','Dinobryon sertularia',
                  'Dinobryon sociale','Glenodinium quadridens',
                  'Microcystis aeruginosa','Mallomonas sp. 1',
                  'Mallomonas sp. 2','Unclassified flagellates',
                  'Peridinium limbatum','Peridinium cinctum',
                  'Peridinium pulsillum','Peridinium wisconsinense',
                  'Chromulina sp.','Selenastrum minutum','Synedra sp.',
                  'Trachelomonas sp.'), BasalNodes(TL84))
}

TestNonBasalNodes <- function()
{
    AssertEqual(c(S=TRUE), IsNonBasalNode(c1))
    AssertEqual('S', NonBasalNodes(c1))
    AssertEqual(c(S=TRUE), IsNonBasalNode(c2))
    AssertEqual('S', NonBasalNodes(c2))
    AssertEqual(c(R=FALSE,C=TRUE), IsNonBasalNode(c3))
    AssertEqual('C', NonBasalNodes(c3))
    AssertEqual(c(R=FALSE,C=TRUE,P=TRUE), IsNonBasalNode(c4))
    AssertEqual(c('C','P'), NonBasalNodes(c4))
    AssertEqual(c(R=FALSE,C=TRUE,O=TRUE), IsNonBasalNode(c5))
    AssertEqual(c('C','O'), NonBasalNodes(c5))
    AssertEqual(c(R=FALSE,C=TRUE,P=TRUE), IsNonBasalNode(c6))
    AssertEqual(c('C','P'), NonBasalNodes(c6))
    AssertEqual(c(A=FALSE,B=TRUE,C=TRUE,D=TRUE,E=TRUE), IsNonBasalNode(c7))
    AssertEqual(c('B','C','D','E'), NonBasalNodes(c7))
}

TestTopLevelNodes <- function()
{
    AssertEqual(c(S=FALSE), IsTopLevelNode(c1))
    AssertEqual(0, length(TopLevelNodes(c1)))
    AssertEqual(c(S=FALSE), IsTopLevelNode(c2))
    AssertEqual(0, length(TopLevelNodes(c2)))
    AssertEqual(c(R=FALSE,C=TRUE), IsTopLevelNode(c3))
    AssertEqual('C', TopLevelNodes(c3))
    AssertEqual(c(R=FALSE,C=FALSE,P=TRUE), IsTopLevelNode(c4))
    AssertEqual('P', TopLevelNodes(c4))
    AssertEqual(c(R=FALSE,C=FALSE,O=TRUE), IsTopLevelNode(c5))
    AssertEqual('O', TopLevelNodes(c5))
    AssertEqual(c(R=FALSE,C=FALSE,P=TRUE), IsTopLevelNode(c6))
    AssertEqual('P', TopLevelNodes(c6))
    AssertEqual(c(A=FALSE,B=FALSE,C=TRUE,D=FALSE,E=FALSE), IsTopLevelNode(c7))
    AssertEqual('C', TopLevelNodes(c7))
    AssertEqual(c(rep(FALSE,55), TRUE), unname(IsTopLevelNode(TL84)))
    AssertEqual('Umbra limi', TopLevelNodes(TL84))
}

TestNonTopLevelNodes <- function()
{
    AssertEqual(c(S=TRUE), IsNonTopLevelNode(c1))
    AssertEqual('S', NonTopLevelNodes(c1))
    AssertEqual(c(S=TRUE), IsNonTopLevelNode(c2))
    AssertEqual('S', NonTopLevelNodes(c2))
    AssertEqual(c(R=TRUE,C=FALSE), IsNonTopLevelNode(c3))
    AssertEqual('R', NonTopLevelNodes(c3))
    AssertEqual(c(R=TRUE,C=TRUE,P=FALSE), IsNonTopLevelNode(c4))
    AssertEqual(c('R','C'), NonTopLevelNodes(c4))
    AssertEqual(c(R=TRUE,C=TRUE,O=FALSE), IsNonTopLevelNode(c5))
    AssertEqual(c('R','C'), NonTopLevelNodes(c5))
    AssertEqual(c(R=TRUE,C=TRUE,P=FALSE), IsNonTopLevelNode(c6))
    AssertEqual(c('R','C'), NonTopLevelNodes(c6))
    AssertEqual(c(A=TRUE,B=TRUE,C=FALSE,D=TRUE,E=TRUE), IsNonTopLevelNode(c7))
    AssertEqual(c('A','B','D','E'), NonTopLevelNodes(c7))
    AssertEqual(c(rep(TRUE,55), FALSE), unname(IsNonTopLevelNode(TL84)))
    AssertEqual(c('Nostoc sp.','Arthrodesmus sp.','Asterionella formosa',
                  'Cryptomonas sp. 1','Cryptomonas sp. 2',
                  'Chroococcus dispersus','Closteriopsis longissimus',
                  'Chrysosphaerella longispina','Dinobryon bavaricum',
                  'Dinobryon cylindricum','Dactylococcopsis fascicularis',
                  'Diceras sp.','Dictyosphaerium pulchellum',
                  'Dinobryon sertularia','Dinobryon sociale',
                  'Glenodinium quadridens','Microcystis aeruginosa',
                  'Mallomonas sp. 1','Mallomonas sp. 2',
                  'Unclassified flagellates','Peridinium limbatum',
                  'Peridinium cinctum','Peridinium pulsillum',
                  'Peridinium wisconsinense','Chromulina sp.',
                  'Rhizosolenia sp.','Selenastrum minutum',
                  'Spinocosmarium sp.','Staurastrum sp.','Synedra sp.',
                  'Trachelomonas sp.','Ascomorpha eucadis','Synchaeta sp.',
                  'Bosmina longirostris','Conochilus (solitary)',
                  'Cyclops varians rubellus',
                  'Diaphanosoma leuchtenbergianum','Daphnia pulex',
                  'Filinia longispina','Conochiloides dossuarius',
                  'Gastropus stylifer','Holopedium gibberum',
                  'Kellicottia sp.','Keratella cochlearis',
                  'Keratella testudo','Leptodiaptomus siciloides',
                  'Orthocyclops modestus','Ploesoma sp.',
                  'Polyarthra vulgaris','Trichocerca multicrinis',
                  'Trichocerca cylindrica','Tropocyclops prasinus',
                  'Chaoborus punctipennis','Phoxinus eos',
                  'Phoxinus neogaeus'), NonTopLevelNodes(TL84))
}

TestIntermediateNodes <- function()
{
    AssertEqual(c(S=FALSE), IsIntermediateNode(c1))
    AssertEqual(0, length(IntermediateNodes(c1)))
    AssertEqual(c(S=FALSE), IsIntermediateNode(c2))
    AssertEqual(0, length(IntermediateNodes(c2)))
    AssertEqual(c(R=FALSE,C=FALSE), IsIntermediateNode(c3))
    AssertEqual(0, length(IntermediateNodes(c3)))
    AssertEqual(c(R=FALSE,C=TRUE,P=FALSE), IsIntermediateNode(c4))
    AssertEqual('C', IntermediateNodes(c4))
    AssertEqual(c(R=FALSE,C=TRUE,O=FALSE), IsIntermediateNode(c5))
    AssertEqual('C', IntermediateNodes(c5))
    AssertEqual(c(R=FALSE,C=TRUE,P=FALSE), IsIntermediateNode(c6))
    AssertEqual('C', IntermediateNodes(c6))
    AssertEqual(c(A=FALSE,B=TRUE,C=FALSE,D=FALSE,E=FALSE), IsIntermediateNode(c7))
    AssertEqual('B', IntermediateNodes(c7))
    AssertEqual(c(rep(FALSE,31), rep(TRUE,24), FALSE), 
                unname(IsIntermediateNode(TL84)))
    AssertEqual(c('Ascomorpha eucadis','Synchaeta sp.',
                  'Bosmina longirostris','Conochilus (solitary)',
                  'Cyclops varians rubellus',
                  'Diaphanosoma leuchtenbergianum','Daphnia pulex',
                  'Filinia longispina','Conochiloides dossuarius',
                  'Gastropus stylifer','Holopedium gibberum',
                  'Kellicottia sp.','Keratella cochlearis',
                  'Keratella testudo','Leptodiaptomus siciloides',
                  'Orthocyclops modestus','Ploesoma sp.',
                  'Polyarthra vulgaris','Trichocerca multicrinis',
                  'Trichocerca cylindrica','Tropocyclops prasinus',
                  'Chaoborus punctipennis','Phoxinus eos',
                  'Phoxinus neogaeus'), 
                IntermediateNodes(TL84))
}

TestIsolatedNodes <- function()
{
    AssertEqual(c(S=TRUE), IsIsolatedNode(c1))
    AssertEqual('S', IsolatedNodes(c1))
    AssertEqual(c(S=TRUE), IsIsolatedNode(c2))
    AssertEqual('S', IsolatedNodes(c2)) # Cannibalism ignored
    AssertEqual(c(R=FALSE,C=FALSE), IsIsolatedNode(c3))
    AssertEqual(0, length(IsolatedNodes(c3)))
    AssertEqual(c(R=FALSE,C=FALSE,P=FALSE), IsIsolatedNode(c4))
    AssertEqual(0, length(IsolatedNodes(c4)))
    AssertEqual(c(R=FALSE,C=FALSE,O=FALSE), IsIsolatedNode(c5))
    AssertEqual(0, length(IsolatedNodes(c5)))
    AssertEqual(c(R=FALSE,C=FALSE,P=FALSE), IsIsolatedNode(c6))
    AssertEqual(0, length(IsolatedNodes(c6)))
    AssertEqual(c(A=FALSE,B=FALSE,C=FALSE,D=TRUE,E=TRUE), IsIsolatedNode(c7))
    AssertEqual(c('D','E'), IsolatedNodes(c7))
    target <- rep(FALSE, 56)
    target[c(3,8,12,26,28,29)] <- TRUE
    AssertEqual(target, unname(IsIsolatedNode(TL84)))
    AssertEqual(c('Asterionella formosa','Chrysosphaerella longispina',
                  'Diceras sp.', 'Rhizosolenia sp.', 'Spinocosmarium sp.',
                  'Staurastrum sp.'), IsolatedNodes(TL84))
}

TestConnectedNodes <- function()
{
    AssertEqual(c(S=FALSE), IsConnectedNode(c1))
    AssertEqual(0, length(ConnectedNodes(c1)))
    AssertEqual(c(S=FALSE), IsConnectedNode(c2))
    AssertEqual(0, length(ConnectedNodes(c2))) # Cannibalism ignored
    AssertEqual(c(R=TRUE,C=TRUE), IsConnectedNode(c3))
    AssertEqual(c('R','C'), ConnectedNodes(c3))
    AssertEqual(c(R=TRUE,C=TRUE,P=TRUE), IsConnectedNode(c4))
    AssertEqual(c('R','C','P'), ConnectedNodes(c4))
    AssertEqual(c(R=TRUE,C=TRUE,O=TRUE), IsConnectedNode(c5))
    AssertEqual(c('R','C','O'), ConnectedNodes(c5))
    AssertEqual(c(R=TRUE,C=TRUE,P=TRUE), IsConnectedNode(c6))
    AssertEqual(c('R','C','P'), ConnectedNodes(c6))
    AssertEqual(c(A=TRUE,B=TRUE,C=TRUE,D=FALSE,E=FALSE), IsConnectedNode(c7))
    AssertEqual(c('A','B','C'), ConnectedNodes(c7))
    target <- rep(TRUE, 56)
    target[c(3,8,12,26,28,29)] <- FALSE
    AssertEqual(target, unname(IsConnectedNode(TL84)))
    AssertEqual(c('Nostoc sp.','Arthrodesmus sp.','Cryptomonas sp. 1',
                  'Cryptomonas sp. 2','Chroococcus dispersus',
                  'Closteriopsis longissimus','Dinobryon bavaricum',
                  'Dinobryon cylindricum','Dactylococcopsis fascicularis',
                  'Dictyosphaerium pulchellum','Dinobryon sertularia',
                  'Dinobryon sociale','Glenodinium quadridens',
                  'Microcystis aeruginosa','Mallomonas sp. 1',
                  'Mallomonas sp. 2','Unclassified flagellates',
                  'Peridinium limbatum','Peridinium cinctum',
                  'Peridinium pulsillum','Peridinium wisconsinense',
                  'Chromulina sp.','Selenastrum minutum','Synedra sp.',
                  'Trachelomonas sp.','Ascomorpha eucadis','Synchaeta sp.',
                  'Bosmina longirostris','Conochilus (solitary)',
                  'Cyclops varians rubellus',
                  'Diaphanosoma leuchtenbergianum','Daphnia pulex',
                  'Filinia longispina','Conochiloides dossuarius',
                  'Gastropus stylifer','Holopedium gibberum',
                  'Kellicottia sp.','Keratella cochlearis',
                  'Keratella testudo','Leptodiaptomus siciloides',
                  'Orthocyclops modestus','Ploesoma sp.',
                  'Polyarthra vulgaris','Trichocerca multicrinis',
                  'Trichocerca cylindrica','Tropocyclops prasinus',
                  'Chaoborus punctipennis','Phoxinus eos',
                  'Phoxinus neogaeus','Umbra limi'), ConnectedNodes(TL84))
}

TestFractionBasalNodes <- function()
{
    AssertEqual(0, FractionBasalNodes(c1))
    AssertEqual(0, FractionBasalNodes(c2))
    AssertEqual(0.5, FractionBasalNodes(c3))
    AssertEqual(1/3, FractionBasalNodes(c4))
    AssertEqual(1/3, FractionBasalNodes(c5))
    AssertEqual(1/3, FractionBasalNodes(c6))
    AssertEqual(1/5, FractionBasalNodes(c7))
    AssertEqual(25/56, FractionBasalNodes(TL84))
}

TestFractionNonBasalNodes <- function()
{
    AssertEqual(1, FractionNonBasalNodes(c1))
    AssertEqual(1, FractionNonBasalNodes(c2))
    AssertEqual(0.5, FractionNonBasalNodes(c3))
    AssertEqual(2/3, FractionNonBasalNodes(c4))
    AssertEqual(2/3, FractionNonBasalNodes(c5))
    AssertEqual(2/3, FractionNonBasalNodes(c6))
    AssertEqual(4/5, FractionNonBasalNodes(c7))
    AssertEqual(31/56, FractionNonBasalNodes(TL84))
}

TestFractionIntermediateNodes <- function()
{
    AssertEqual(0, FractionIntermediateNodes(c1))
    AssertEqual(0, FractionIntermediateNodes(c2))
    AssertEqual(0, FractionIntermediateNodes(c3))
    AssertEqual(1/3, FractionIntermediateNodes(c4))
    AssertEqual(1/3, FractionIntermediateNodes(c5))
    AssertEqual(1/3, FractionIntermediateNodes(c6))
    AssertEqual(1/5, FractionIntermediateNodes(c7))
    AssertEqual(24/56, FractionIntermediateNodes(TL84))
}

TestFractionTopLevelNodes <- function()
{
    AssertEqual(0, FractionTopLevelNodes(c1))
    AssertEqual(0, FractionTopLevelNodes(c2))
    AssertEqual(0.5, FractionTopLevelNodes(c3))
    AssertEqual(1/3, FractionTopLevelNodes(c4))
    AssertEqual(1/3, FractionTopLevelNodes(c5))
    AssertEqual(1/3, FractionTopLevelNodes(c6))
    AssertEqual(1/5, FractionTopLevelNodes(c7))
    AssertEqual(1/56, FractionTopLevelNodes(TL84))
}

TestFractionNonTopLevelNodes <- function()
{
    AssertEqual(1, FractionNonTopLevelNodes(c1))
    AssertEqual(1, FractionNonTopLevelNodes(c2))
    AssertEqual(0.5, FractionNonTopLevelNodes(c3))
    AssertEqual(2/3, FractionNonTopLevelNodes(c4))
    AssertEqual(2/3, FractionNonTopLevelNodes(c5))
    AssertEqual(2/3, FractionNonTopLevelNodes(c6))
    AssertEqual(4/5, FractionNonTopLevelNodes(c7))
    AssertEqual(55/56, FractionNonTopLevelNodes(TL84))
}

TestFractionIsolatedNodes <- function()
{
    AssertEqual(1, FractionIsolatedNodes(c1))
    AssertEqual(1, FractionIsolatedNodes(c2))
    AssertEqual(0, FractionIsolatedNodes(c3))
    AssertEqual(0, FractionIsolatedNodes(c4))
    AssertEqual(0, FractionIsolatedNodes(c5))
    AssertEqual(0, FractionIsolatedNodes(c6))
    AssertEqual(2/5, FractionIsolatedNodes(c7))
    AssertEqual(6/56, FractionIsolatedNodes(TL84))
}

TestFractionConnectedNodes <- function()
{
    AssertEqual(0, FractionConnectedNodes(c1))
    AssertEqual(0, FractionConnectedNodes(c2))
    AssertEqual(1, FractionConnectedNodes(c3))
    AssertEqual(1, FractionConnectedNodes(c4))
    AssertEqual(1, FractionConnectedNodes(c5))
    AssertEqual(1, FractionConnectedNodes(c6))
    AssertEqual(3/5, FractionConnectedNodes(c7))
    AssertEqual(50/56, FractionConnectedNodes(TL84))
}

TestFractions <- function()
{
    # These should sum to 1
    for(community in list(c1,c2,c3,c4,c5,c6,c7,TL84))
    {
        print(community)
        AssertEqual(1, sum(FractionBasalNodes(community), 
                           FractionIntermediateNodes(community), 
                           FractionTopLevelNodes(community), 
                           FractionIsolatedNodes(community)))

        AssertEqual(1, sum(FractionConnectedNodes(community), 
                           FractionIsolatedNodes(community)))

        AssertEqual(1, sum(FractionBasalNodes(community), 
                           FractionNonBasalNodes(community)))

        AssertEqual(1, sum(FractionTopLevelNodes(community), 
                           FractionNonTopLevelNodes(community)))
    }
}

TestResourcesByNode <- function()
{
    AssertEqual(list(S=vector('character')), ResourcesByNode(c1))
    AssertEqual(list(S='S'), ResourcesByNode(c2))
    AssertEqual(list(R=vector('character'), C='R'), ResourcesByNode(c3))
    AssertEqual(list(R=vector('character'), C='R', P='C'), ResourcesByNode(c4))
    AssertEqual(list(R=vector('character'), C='R', O=c('R','C')), 
                ResourcesByNode(c5))
    AssertEqual(list(R=vector('character'), C='R', P='C'), ResourcesByNode(c6))
}

TestResourcesOfNodes <- function()
{
    AssertEqual(vector('character'), ResourcesOfNodes(c1, 'S'))
    AssertEqual('S', ResourcesOfNodes(c2, 'S'))
    AssertEqual(vector('character'), ResourcesOfNodes(c3, 'R'))
    AssertEqual('R', ResourcesOfNodes(c3, 'C'))
    AssertEqual(vector('character'), ResourcesOfNodes(c4, 'R'))
    AssertEqual('R', ResourcesOfNodes(c4, 'C'))
    AssertEqual('C', ResourcesOfNodes(c4, 'P'))

    # More than one node
    AssertEqual(list(R=vector('character'), C='R'), 
                        ResourcesOfNodes(c4, c('R','C')))
    AssertEqual(list(R=vector('character'), C='R'), 
                        ResourcesOfNodes(c4, 1:2))

    # Bad nodes
    AssertRaises(ResourcesOfNodes(c1, 0))
    AssertRaises(ResourcesOfNodes(c1, 2))
    AssertRaises(ResourcesOfNodes(c1, ''))
    AssertRaises(ResourcesOfNodes(c1, 'x'))
}

TestConsumersByNode <- function()
{
    AssertEqual(list(S=vector('character')), ConsumersByNode(c1))
    AssertEqual(list(S='S'), ConsumersByNode(c2))
    AssertEqual(list(R='C', C=vector('character')), 
                        ConsumersByNode(c3))
    AssertEqual(list(R='C', C='P', P=vector('character')), 
                        ConsumersByNode(c4))
    AssertEqual(list(R=c('C','O'), C='O', O=vector('character')), 
                        ConsumersByNode(c5))
    AssertEqual(list(R='C', C='P', P=vector('character')), 
                        ConsumersByNode(c6))
}

TestConsumersOfNodes <- function()
{
    AssertEqual(vector('character'), ConsumersOfNodes(c1, 'S'))
    AssertEqual('S', ConsumersOfNodes(c2, 'S'))
    AssertEqual('C', ConsumersOfNodes(c3, 'R'))
    AssertEqual(vector('character'), ConsumersOfNodes(c3, 'C'))
    AssertEqual('C', ConsumersOfNodes(c4, 'R'))
    AssertEqual('P', ConsumersOfNodes(c4, 'C'))
    AssertEqual(vector('character'), ConsumersOfNodes(c4, 'P'))

    # More than one node
    AssertEqual(list(R='C', C='P'), 
                        ConsumersOfNodes(c4, c('R','C')))
    AssertEqual(list(R='C', C='P'), 
                        ConsumersOfNodes(c4, 1:2))

    # Bad nodes
    AssertRaises(ConsumersOfNodes(c1, 0))
    AssertRaises(ConsumersOfNodes(c1, 2))
    AssertRaises(ConsumersOfNodes(c1, ''))
    AssertRaises(ConsumersOfNodes(c1, 'x'))
}

TestResourcesAndConsumersByNode <- function()
{
    AssertEqual(list(S=vector('character')), ResourcesAndConsumersByNode(c1))
    AssertEqual(list(S='S'), ResourcesAndConsumersByNode(c2))
    AssertEqual(list(R='C', C='R'), ResourcesAndConsumersByNode(c3))
    AssertEqual(list(R='C', C=c('R','P'), P='C'), 
                ResourcesAndConsumersByNode(c4))
    AssertEqual(list(R=c('C','O'), C=c('R','O'), O=c('R','C')), 
                ResourcesAndConsumersByNode(c5))
    AssertEqual(list(R='C', C=c('R','P'), P='C'), 
                ResourcesAndConsumersByNode(c6))
    AssertEqual(list(A=c('A','B'), B=c('A','C'), C=c('B','C'), D='D', 
                     E=vector(mode='character')), 
                ResourcesAndConsumersByNode(c7))
}

TestNumberOfTrophicLinks <- function()
{
    AssertEqual(0, NumberOfTrophicLinks(c1))
    AssertEqual(1, NumberOfTrophicLinks(c2))
    AssertEqual(1, NumberOfTrophicLinks(c3))
    AssertEqual(2, NumberOfTrophicLinks(c4))
    AssertEqual(3, NumberOfTrophicLinks(c5))
    AssertEqual(2, NumberOfTrophicLinks(c6))
    AssertEqual(5, NumberOfTrophicLinks(c7))
    AssertEqual(5, NumberOfTrophicLinks(c7))
    AssertEqual(269, NumberOfTrophicLinks(TL84))
    AssertEqual(264, NumberOfTrophicLinks(RemoveCannibalisticLinks(TL84)))
}

TestLinkageDensity <- function()
{
    AssertEqual(0, LinkageDensity(c1))
    AssertEqual(1, LinkageDensity(c2))
    AssertEqual(0.5, LinkageDensity(c3))
    AssertEqual(2/3, LinkageDensity(c4))
    AssertEqual(1, LinkageDensity(c5))
    AssertEqual(4.80357142857142882519, LinkageDensity(TL84))
}

TestDegree <- function()
{
    AssertEqual(c(S=0), Degree(c1))
    AssertEqual(c(S=2), Degree(c2))
    AssertEqual(InDegree(c2)+OutDegree(c2), Degree(c2))
    AssertEqual(c(R=1,C=1), Degree(c3))
    AssertEqual(InDegree(c3)+OutDegree(c3), Degree(c3))
    AssertEqual(c(R=1,C=2,P=1), Degree(c4))
    AssertEqual(InDegree(c4)+OutDegree(c4), Degree(c4))
    AssertEqual(c(R=2,C=2,O=2), Degree(c5))
    AssertEqual(InDegree(c5)+OutDegree(c5), Degree(c5))
    AssertEqual(c(4,3,0,18,5,18,3,0,1,4,4,0,5,2,5,6,1,2,3,18,1,2,5,1,
                  18,0,18,0,0,1,18,11,11,15,11,23,18,32,10,11,10,22,10,
                  11,10,23,22,10,10,10,10,23,26,10,10,13), unname(Degree(TL84)))
    AssertEqual(InDegree(TL84)+OutDegree(TL84), Degree(TL84))
}

TestDegreeDistribution <- function()
{
    AssertEqual(c('0'=1), DegreeDistribution(c1))
    AssertEqual(c('0'=0, '1'=0, '2'=1), DegreeDistribution(c2))
    AssertEqual(c('0'=0, '1'=1), DegreeDistribution(c3))
    AssertEqual(c('0'=0, '1'=2/3, '2'=1/3), DegreeDistribution(c4))
    AssertEqual(c('0'=0, '1'=0, '2'=1), DegreeDistribution(c5))
    AssertEqual(c('0'=0, '1'=2/3, '2'=1/3), DegreeDistribution(c6))
}

TestDirectedConnectance <- function()
{
    AssertEqual(0, DirectedConnectance(c1))
    AssertEqual(1, DirectedConnectance(c2))
    AssertEqual(0.25, DirectedConnectance(c3))
    AssertEqual(2/9, DirectedConnectance(c4))
    AssertEqual(1/3, DirectedConnectance(c5))
    AssertEqual(0.08577806122448979054, DirectedConnectance(TL84))
}

TestCannibals <- function()
{
    AssertEqual(c(S=FALSE), IsCannibal(c1))
    AssertEqual(0, length(Cannibals(c1)))
    AssertEqual(0, FractionCannibalistic(c1))
    AssertEqual(c(S=TRUE), IsCannibal(c2))
    AssertEqual('S', Cannibals(c2))
    AssertEqual(1, FractionCannibalistic(c2))
    AssertEqual(c(R=FALSE,C=FALSE), IsCannibal(c3))
    AssertEqual(0, length(Cannibals(c3)))
    AssertEqual(0, FractionCannibalistic(c3))
    AssertEqual(c(R=FALSE,C=FALSE,P=FALSE), IsCannibal(c4))
    AssertEqual(0, length(Cannibals(c4)))
    AssertEqual(0, FractionCannibalistic(c4))
    AssertEqual(c(R=FALSE,C=FALSE,O=FALSE), IsCannibal(c5))
    AssertEqual(0, length(Cannibals(c5)))
    AssertEqual(0, FractionCannibalistic(c5))
    AssertEqual(c(A=TRUE,B=FALSE,C=TRUE,D=TRUE,E=FALSE), IsCannibal(c7))
    AssertEqual(c('A','C','D'), Cannibals(c7))
    AssertEqual(3/5, FractionCannibalistic(c7))
    AssertEqual(c("Cyclops varians rubellus", "Orthocyclops modestus",
                  "Tropocyclops prasinus", "Chaoborus punctipennis", 
                  "Umbra limi"), Cannibals(TL84))
    AssertEqual(5/56, FractionCannibalistic(TL84))
}

TestTrophicChainsStats <- function()
{
    AssertEqual(NULL, TrophicChainsStats(c1))
    AssertEqual(NULL, TrophicChainsStats(c2))
    AssertEqual(1, TrophicChainsStats(c3)$chain.lengths)
    AssertEqual(matrix(c(1,0,0,1), nrow=2, dimnames=list(c('R','C'), NULL)), 
                TrophicChainsStats(c3)$node.pos.counts)
    AssertEqual(2, TrophicChainsStats(c4)$chain.lengths)
    AssertEqual(matrix(c(1,0,0,0,1,0,0,0,1), nrow=3, 
                       dimnames=list(c('R','C','P'), NULL)), 
                TrophicChainsStats(c4)$node.pos.counts)
    AssertEqual(c(1,2), TrophicChainsStats(c5)$chain.lengths)
    AssertEqual(matrix(c(2,0,0,0,1,1,0,0,1), nrow=3, 
                       dimnames=list(c('R','C','O'), NULL)), 
                TrophicChainsStats(c5)$node.pos.counts)
    AssertEqual(2, TrophicChainsStats(c6)$chain.lengths)
    AssertEqual(matrix(c(1,0,0,0,1,0,0,0,1), nrow=3, 
                       dimnames=list(c('R','C','P'), NULL)), 
                TrophicChainsStats(c6)$node.pos.counts)
    AssertEqual(2, TrophicChainsStats(c7)$chain.lengths)
    AssertEqual(matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0), nrow=5, 
                       dimnames=list(c('A','B','C','D','E'), NULL)), 
                TrophicChainsStats(c7)$node.pos.counts)

    # Simpler tests for larger communities                              sum pos 
    #                       Mean chain length       n.chains, longest   count
    to.test <- list(
      list(TL84,             4.83500334001336007361,    5988,   8,      34940), 
      list(TL86,             4.84470882905447730593,    1597,   7,       9334), 
      list(YthanEstuary,     5.08319118225954635903,    7621,  10,      46360), 
      list(SkipwithPond,    10.4850243487305565537,  2538120,  17,   29150370), 
      list(BroadstoneStream, 5.63949610967024828057,    5398,  10,      35840))
    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        print(community)

        s <- TrophicChainsStats(community)
        AssertEqual(to.test[[index]][[2]], mean(s$chain.lengths))
        AssertEqual(to.test[[index]][[3]], length(s$chain.length))
        AssertEqual(to.test[[index]][[4]], ncol(s$node.pos.counts))
        AssertEqual(to.test[[index]][[5]], sum(s$node.pos.counts))
    }
}

TestTrophicChainsStatsOverflow <- function()
{
    community <- TL84
    alist <- cheddar:::.CAdjacencyList(community, ConsumersByNode(community))

    # 1 if basal, 0 otherwise
    is.basal <- as.integer(IsBasalNode(community))

    # Outputs to be filled by the C function
    n.chains <- as.integer(0)
    longest <- as.integer(0)
    status <- as.integer(-1)

    res <- .C('trophic_chains_size', 
              as.integer(alist), 
              as.integer(length(alist)), 
              as.integer(is.basal), 
              as.integer(nrow(alist)), 
              as.integer(1), # test_overflow
              as.integer(0), # max_queue
              n.chains=n.chains,
              longest=longest,
              status=status, 
              PACKAGE='cheddar', NAOK=TRUE, DUP=FALSE)

    AssertEqual(-1, res$status)
}

TestTrophicChains <- function()
{
    CheckChains <- function(community, chains)
    {
        # Helper function
        # Chains of length 1 would indicate isolated node
        AssertEqual(FALSE, any(1==apply(chains, 1, function(r) sum(''!=r))))

        # Should be no duplicated rows
        AssertEqual(FALSE, any(duplicated(chains, MARGIN=1)))

        # No node should appear more than once in a chain
        AssertEqual(FALSE, 
                   any(apply(chains, 1, function(r) any(duplicated(r[''!=r])))))

        # Each chain should start with a basal node
        AssertEqual(TRUE, all(chains[,1] %in% BasalNodes(community)))

        if(any(IsTopLevelNode(community)))
        {
            # Each chain should end with a top-level node
            last <- apply(chains, 1, function(row) row[max(which(""!=row))])
            AssertEqual(TRUE, all(last %in% TopLevelNodes(community)))
        }

        # Make sure that every link in every chain appears in the predation 
        # matrix
        pm <- PredationMatrix(community)
        apply(chains, 1, function(chain)
        {
            last.in.chain <- tail(which(''!=chain), 1)
            from <- chain[1:(last.in.chain-1)]
            to <- chain[2:last.in.chain]
            AssertEqual(rep(1, length(to)), as.vector(diag(pm[from,to])))
        })

        # How could this not be the case?
        AssertEqual(TRUE, ncol(chains)<=NumberOfNodes(community))
    }

    # Some extra test cases
    test1 <- Community(nodes=data.frame(node=c('A', 'B', 'C')), 
                       trophic.links=cbind(resource=c('A', 'B', 'C'), 
                                           consumer=c('B', 'C', 'B')), 
                       properties=list(title='test1'))

    test2 <- Community(nodes=data.frame(node=c('A', 'B', 'C')), 
                       trophic.links=cbind(resource=c('A', 'B', 'A', 'C'), 
                                           consumer=c('B', 'C', 'C', 'B')), 
                       properties=list(title='test2'))

    test3 <- Community(nodes=data.frame(node=c('A', 'B', 'C')), 
                       trophic.links=cbind(resource=c('A', 'B', 'A'), 
                                           consumer=c('B', 'C', 'C')), 
                       properties=list(title='test2'))


    to.test <- list(list(c1,NULL), 
      list(c2,NULL), 
      list(c3,data.frame(Node.1='R',Node.2='C', stringsAsFactors=FALSE)), 
      list(c4,data.frame(Node.1='R',Node.2='C',Node.3='P', stringsAsFactors=FALSE)), 
      list(c5,data.frame(Node.1=c('R','R'), 
                         Node.2=c('O','C'), 
                         Node.3=c('', 'O'), stringsAsFactors=FALSE)), 
      list(c6,data.frame(Node.1='R',Node.2='C',Node.3='P', stringsAsFactors=FALSE)), 
      list(c7,data.frame(Node.1='A',Node.2='B',Node.3='C', stringsAsFactors=FALSE)),
      list(test1,data.frame(Node.1='A',Node.2='B',Node.3='C', stringsAsFactors=FALSE)), 
      list(test2,data.frame(Node.1=c('A','A'),
                            Node.2=c('B','C'), 
                            Node.3=c('C','B'), stringsAsFactors=FALSE)),
      list(test3,data.frame(Node.1=c('A','A'),
                            Node.2=c('C','B'), 
                            Node.3=c( '','C'), stringsAsFactors=FALSE)))
    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        expected <- to.test[[index]][[2]]
        print(community)

        chains <- TrophicChains(community)
        AssertEqual(expected, chains)
        if(!is.null(chains))
        {
            CheckChains(community, chains)
        }
    }

    # Webs with many nodes - don't check each link, just dims and first chain
    to.test <- list(
      list(TL84, c(5988, 8), list('Nostoc sp.','Diaphanosoma leuchtenbergianum',
                                  'Umbra limi')),
      list(TL86, c(1597, 7), list('Ankyra judayi','Daphnia rosea',
                                  'Chaoborus punctipennis',
                                  'Micropterus salmoides')),
      list(YthanEstuary, c(7621, 10), list('Diatoms','Crangon crangon',
                                           'Somateria mollissima')),
      list(BroadstoneStream, c(5398, 10), list('Potamophylax cingulatus', 
                                               'Cordulegaster boltonii')), 
      list(Benguela, c(15430, 14), list('Phytoplankton', 'Anchovy', 'Seals', 
                                        'Sharks')))
    to.test <- head(to.test, -1)

    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        expected.dim <- to.test[[index]][[2]]
        expected.first.chain <- to.test[[index]][[3]]
        print(community)

        chains <- TrophicChains(community)
        AssertEqual(expected.dim, dim(chains))
        AssertEqual(expected.first.chain, 
                    unname(chains[1,1:length(expected.first.chain),drop=TRUE]))
        CheckChains(community, chains)
    }
}

TestThreeNodeChains <- function()
{
    CheckChains <- function(community, chains)
    {
        # Helper function
        # Sanity checks
        # Should be no duplicated rows
        AssertEqual(FALSE, any(duplicated(chains, MARGIN=1)))

        # No node should appear more than once in a chain except as bottom and 
        # top.
        AssertEqual(FALSE, 
                    any(apply(chains, 1, function(r) r[1]==r[2] | r[2]==r[3])))

        # Make sure that every link in every chain appears in the predation 
        # matrix
        pm <- PredationMatrix(community)
        apply(chains, 1, function(chain)
        {
            from <- head(chain, -1)
            to <- tail(chain, -1)
            AssertEqual(rep(1, length(to)), diag(pm[from,to]))
        })
    }

    to.test <- list(
      list(c1, NULL), 
      list(c2, NULL), 
      list(c3, NULL), 
      list(c4, data.frame(bottom='R',intermediate='C',top='P')), 
      list(c5, data.frame(bottom='R',intermediate='C',top='O')), 
      list(c6, data.frame(bottom='R',intermediate='C',top='P')),
      list(c7, data.frame(bottom='A',intermediate='B',top='C')))
    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        expected <- to.test[[index]][[2]]
        print(community)

        chains <- ThreeNodeChains(community)
        AssertEqual(expected, chains)
        if(!is.null(chains))
        {
            CheckChains(community, chains)
        }

        chains <- ThreeNodeChains(community, exclude.loops=TRUE)
        AssertEqual(expected, chains)
        if(!is.null(chains))
        {
            CheckChains(community, chains)
        }
    }

    # Webs with many nodes - don't check each link, just nrows
    to.test <- list(list(TL84, 1044, 1042), 
                    list(TL86, 651, 645),
                    list(YthanEstuary, 1616, 1616),
                    list(SkipwithPond, 2433, 2401))
    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        expected.rows <- to.test[[index]][[2]]
        expected.rows.no.loops <- to.test[[index]][[3]]
        print(community)

        chains <- ThreeNodeChains(community)
        AssertEqual(expected.rows, nrow(chains))
        CheckChains(community, chains)

        chains <- ThreeNodeChains(community, exclude.loops=TRUE)
        AssertEqual(expected.rows.no.loops, nrow(chains))
        CheckChains(community, chains)
    }
}

TestTrophicHeight <- function()
{
    AssertEqual(c(S=NA), TrophicHeight(c1, include.isolated=TRUE))
    AssertEqual(c(S=NA), TrophicHeight(c1, include.isolated=FALSE))
    AssertEqual(c(S=NA), TrophicHeight(c2, include.isolated=TRUE))
    AssertEqual(c(S=NA), TrophicHeight(c2, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2), TrophicHeight(c3, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2), TrophicHeight(c3, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,P=3), TrophicHeight(c4, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2,P=3), TrophicHeight(c4, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,O=2.5), TrophicHeight(c5, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2,O=2.5), TrophicHeight(c5, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,P=3), TrophicHeight(c6, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2,P=3), TrophicHeight(c6, include.isolated=FALSE))
    AssertEqual(c(A=1,B=2,C=3,D=1,E=1), TrophicHeight(c7, include.isolated=TRUE))
    AssertEqual(c(A=1,B=2,C=3,D=NA,E=NA), 
                TrophicHeight(c7, include.isolated=FALSE))

    # Lower tolerance required for Mac Pro in GIS lab.
    AssertEqual(c(rep(1,31),2,2,2,2,10/3,2,2.41860465116279055309,2,2,
                  2,2,2,2,2,2,4.2,2,2,2,2,10/3,
                  4.60252672497570447518, 5.16833667334669311089, 
                  5.16833667334669311089, 5.83500334001336007361), 
                unname(TrophicHeight(TL84)), 
                tolerance=1e-6)
    AssertEqual(c(1,1,NA,1,1,1,1,NA,1,1,1,NA,rep(1,13),NA,1,NA,NA,1,1,
                  2,2,2,2,10/3,2,2.41860465116279055309,2,2,
                  2,2,2,2,2,2,4.2,2,2,2,2,10/3,
                  4.60252672497570447518, 5.16833667334669311089, 
                  5.16833667334669311089, 5.83500334001336007361), 
                unname(TrophicHeight(TL84, include.isolated=FALSE)), 
                tolerance=1e-6)
}

TestTrophicLevel <- function()
{
    AssertEqual(c(S=NA), PreyAveragedTrophicLevel(c1))
    AssertEqual(c(S=NA), PreyAveragedTrophicLevel(c2))
    AssertEqual(c(R=1,C=2), PreyAveragedTrophicLevel(c3))
    AssertEqual(c(R=1,C=2,P=3), PreyAveragedTrophicLevel(c4))
    AssertEqual(c(R=1,C=2,O=2.5), PreyAveragedTrophicLevel(c5))
    AssertEqual(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,2,2,2,2,3.14285714285714,2,
                          2.10714285714286,2,2,2,2,2,2,2,2,3.20535714285714,
                          2,2,2,2,3.14285714285714,3.17134353741497,
                          3.52995086923658,3.52995086923658,
                          3.80267814196386), 
                unname(PreyAveragedTrophicLevel(TL84)))

    AssertEqual(c(S=NA), PreyAveragedTrophicLevel(c1, include.isolated=FALSE))
    AssertEqual(c(S=NA), PreyAveragedTrophicLevel(c2, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2), PreyAveragedTrophicLevel(c3, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,P=3), PreyAveragedTrophicLevel(c4, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,O=2.5), PreyAveragedTrophicLevel(c5, include.isolated=FALSE))

    # c7 is not energetically feasible so all trophic levels are NA
    AssertEqual(c(A=NA,B=NA,C=NA,D=NA,E=NA), 
                PreyAveragedTrophicLevel(c7, include.isolated=FALSE))
    AssertEqual(c(A=NA,B=NA,C=NA,D=NA,E=NA), 
                PreyAveragedTrophicLevel(c7, include.isolated=TRUE))

    AssertEqual(c(1,1,NA,1,1,1,1,NA,1,1,1,NA,1,1,1,1,1,1,1,1,1,1,1,1,
                          1,NA,1,NA,NA,1,1,2,2,2,2,3.14285714285714,2,
                          2.10714285714286,2,2,2,2,2,2,2,2,3.20535714285714,
                          2,2,2,2,3.14285714285714,3.17134353741497,
                          3.52995086923658,3.52995086923658,
                          3.80267814196386), 
                unname(PreyAveragedTrophicLevel(TL84, include.isolated=FALSE)))
}

TestTrophicSpecies <- function()
{
    AssertEqual(c(S=NA), TrophicSpecies(c1, include.isolated=FALSE))
    AssertEqual(c(S=1), TrophicSpecies(c1, include.isolated=TRUE))
    AssertEqual(c(S=NA), TrophicSpecies(c2, include.isolated=FALSE))
    AssertEqual(c(S=1), TrophicSpecies(c2, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2), TrophicSpecies(c3, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2), TrophicSpecies(c3, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2,P=3), TrophicSpecies(c4, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,P=3), TrophicSpecies(c4, include.isolated=TRUE))
    AssertEqual(c(R=1,C=2,O=3), TrophicSpecies(c5, include.isolated=FALSE))
    AssertEqual(c(R=1,C=2,O=3), TrophicSpecies(c5, include.isolated=TRUE))

    # From Jonsson et al 2005 AER. Isolated species assigned NA.
    target <- c(1,2,NA,3,4,3,5,NA,6,1,1,NA,4,7,4,8,6,7,2,3,6,7,
                4,6,3,NA,3,NA,NA,6,3,9,9,10,11,12,13,14,15,11,
                15,16,15,9,15,17,18,15,15,15,15,12,19,20,20,21)
    AssertEqual(target, unname(TrophicSpecies(TL84, include.isolated=FALSE)))

    # Isolated species included - these are given a trophic level of 1.
    target <- c(1,2,3,4,5,4,6,3,7,1,1,3,5,8,5,9,7,8,2,4,7,8,5,7,
                4,3,4,3,3,7,4,10,10,11,12,13,14,15,16,12,16,17,16,
                10,16,18,19,16,16,16,16,13,20,21,21,22)
    AssertEqual(target, unname(TrophicSpecies(TL84, include.isolated=TRUE)))
}

TestTrophicLinksForNodes <- function()
{
    AssertEqual(NULL, TrophicLinksForNodes(c1, 'S'))
    AssertEqual(data.frame(resource='S', consumer='S', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c2, 'S'))

    AssertEqual(data.frame(resource='R', consumer='C', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c3, 'R'))
    AssertEqual(data.frame(resource='R', consumer='C', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c3, 'C'))
    AssertEqual(data.frame(resource='R', consumer='C', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c3, c('R','C')))

    AssertEqual(data.frame(resource='R', consumer='C', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c4, 'R'))
    AssertEqual(data.frame(resource=c('R','C'), consumer=c('C','P'), 
                           stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c4, 'C'))
    AssertEqual(data.frame(resource='C', consumer='P', stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c4, 'P'))
    AssertEqual(data.frame(resource=c('R','C'), consumer=c('C','P'), 
                           stringsAsFactors=FALSE), 
                TrophicLinksForNodes(c4, c('R','C','P')))
}

TestResolveTrophicLinksToNodeNames <- function()
{
    f <- cheddar:::.ResolveTrophicLinksToNodeNames
    expected = matrix(c('S', 'S'), nrow=1, 
                      dimnames=list(NULL, c('resource', 'consumer')))
    AssertEqual(expected, f(c1, data.frame(resource='S', consumer='S')))
    AssertEqual(expected, f(c1, data.frame(resource=1, consumer='S')))
    AssertEqual(expected, f(c1, data.frame(resource=1, consumer=1)))
    AssertEqual(expected, f(c1, data.frame(resource='S', consumer=1)))
    AssertRaises(f(c3, data.frame(resource='x', consumer='R')))
    AssertRaises(f(c3, data.frame(resource='R', consumer='')))
    AssertRaises(f(c3, matrix(c(0, 1), ncol=2)))
    AssertRaises(f(c3, matrix(c(0, 1), ncol=2)))
}

TestResolveTrophicLinksToRowIndices <- function()
{
    f <- cheddar:::.ResolveTrophicLinksToRowIndices
    AssertEqual(NA, f(c1, data.frame(resource='S', consumer='S')))
    AssertEqual(NA, f(c1, data.frame(resource=1, consumer=1)))
    AssertEqual(1, f(c2, data.frame(resource='S', consumer='S')))
    AssertEqual(1, f(c2, data.frame(resource=1, consumer=1)))
    AssertEqual(1, f(c3, data.frame(resource='R', consumer='C')))
    AssertEqual(NA, f(c3, data.frame(resource='R', consumer='R')))
    AssertEqual(NA, f(c3, data.frame(resource='C', consumer='C')))
    AssertEqual(NA, f(c3, data.frame(resource='C', consumer='R')))
    AssertEqual(c(1,NA,NA,NA), 
                f(c3, data.frame(resource=c('R','R','C','C'), 
                                 consumer=c('C','R','C','R'))))
    AssertRaises(f(c3, data.frame(resource='x', consumer='R')))
    AssertRaises(f(c3, data.frame(resource='R', consumer='')))
    AssertRaises(f(c3, matrix(resource=0, consumer=1)))
    AssertRaises(f(c3, matrix(resource=1, consumer=0)))
}

TestDoLinksExist <- function()
{
    f <- cheddar:::.DoLinksExist
    AssertEqual(FALSE, f(c1, data.frame(resource='S', consumer='S')))
    AssertEqual(FALSE, f(c1, data.frame(resource=1, consumer=1)))
    AssertEqual(TRUE, f(c2, data.frame(resource='S', consumer='S')))
    AssertEqual(TRUE, f(c2, data.frame(resource=1, consumer=1)))
    AssertEqual(TRUE, f(c3, data.frame(resource='R', consumer='C')))
    AssertEqual(FALSE, f(c3, data.frame(resource='R', consumer='R')))
    AssertEqual(FALSE, f(c3, data.frame(resource='C', consumer='C')))
    AssertEqual(FALSE, f(c3, data.frame(resource='C', consumer='R')))
    AssertEqual(c(TRUE,FALSE,FALSE,FALSE), 
                f(c3, data.frame(resource=c('R','R','C','C'), 
                                 consumer=c('C','R','C','R'))))
}

TestRemoveCannibalisticLinks <- function()
{
    AssertEqual(0, NumberOfTrophicLinks(c1))
    AssertEqual(0, NumberOfTrophicLinks(RemoveCannibalisticLinks(c1)))
    AssertEqual(1, NumberOfTrophicLinks(c2))
    AssertEqual(0, NumberOfTrophicLinks(RemoveCannibalisticLinks(c2)))
    AssertEqual(1, NumberOfTrophicLinks(RemoveCannibalisticLinks(c3)))
    AssertEqual(2, NumberOfTrophicLinks(RemoveCannibalisticLinks(c4)))
    AssertEqual(3, NumberOfTrophicLinks(RemoveCannibalisticLinks(c5)))
    AssertEqual(2, NumberOfTrophicLinks(RemoveCannibalisticLinks(c7)))
    AssertEqual(269, NumberOfTrophicLinks(TL84))
    AssertEqual(5, length(Cannibals(TL84)))
    AssertEqual(264, NumberOfTrophicLinks(RemoveCannibalisticLinks(TL84)))
    AssertEqual(0, length(Cannibals(RemoveCannibalisticLinks(TL84))))
}

TestShortestPathLengths <- function()
{
    # Helper function that creates expected matrices
    F <- function(data, nodes)
    {
        m <- matrix(data, ncol=length(nodes), dimnames=list(nodes, nodes))
        mode(m) <- 'numeric' # Necessary for a 1 x 1 matrix that contains only NA
        return (m)
    }

    AssertEqual(F(0, 'S'), ShortestPaths(c1))
    AssertEqual(F(0, 'S'), ShortestPaths(c2))
    AssertEqual(F(c(0,1, 1,0), c('R','C')), ShortestPaths(c3))
    AssertEqual(F(c(0,1,2, 1,0,1, 2,1,0), c('R','C','P')), ShortestPaths(c4))
    AssertEqual(F(c(0,1,1, 1,0,1, 1,1,0), c('R','C','O')), ShortestPaths(c5))
    AssertEqual(F(c(0,1,2, 1,0,1, 2,1,0), c('R','C','P')), ShortestPaths(c6))
}

TestSumDietGap <- function()
{
    # Test the C function directly
    f <- function(pm, expected)
    {
        n <- ncol(pm)
        sum_diet_gaps <- as.integer(0)
        status <- as.integer(-1)

        res <- .C('sum_diet_gaps', 
                  as.integer(pm), 
                  as.integer(n), 
                  as.integer((1:n) - 1),   # 0-indexed
                  sum_diet_gaps=sum_diet_gaps, 
                  status=status)
        AssertEqual(res$sum_diet_gaps, expected)
        return (res)
    }

    f(matrix(0, ncol=1), 0)
    f(matrix(1, ncol=1), 0)
    f(matrix(c(0,0,0,0), ncol=2), 0)
    f(matrix(c(1,1,1,1), ncol=2), 0)
    f(matrix(c(1,0,1,0), ncol=2), 0)
    f(matrix(c(0,1,0,1), ncol=2), 0)

    f(matrix(c(1, 1, 0, 
               0, 1, 1,
               1, 0, 1), ncol=3, byrow=TRUE), 1)
    f(matrix(c(1, 1, 0, 
               0, 0, 1,
               1, 1, 1), ncol=3, byrow=TRUE), 2)
    f(matrix(c(1, 1, 1, 
               0, 0, 0,
               1, 1, 1), ncol=3, byrow=TRUE), 3)
    f(matrix(c(1, 1, 1, 
               0, 0, 0,
               0, 0, 0), ncol=3, byrow=TRUE), 0)
    f(matrix(c(1, 1, 1, 
               0, 1, 1,
               0, 0, 1), ncol=3, byrow=TRUE), 0)
    f(matrix(c(0, 0, 1, 
               0, 1, 1,
               1, 1, 1), ncol=3, byrow=TRUE), 0)

    f(matrix(c(1, 0, 0, 1, 
               1, 0, 0, 0,
               0, 0, 0, 1,
               1, 0, 0, 1), ncol=4, byrow=TRUE), 2)

    f(matrix(c(1, 1, 0, 0, 0, 0, 
               1, 1, 1, 1, 0, 0, 
               1, 1, 1, 1, 0, 0, 
               0, 0, 0, 0, 0, 0, 
               1, 1, 1, 1, 0, 0, 
               0, 1, 0, 1, 0, 0), ncol=6, byrow=TRUE), 4)

    f(matrix(c(1, 1, 0, 1, 1, 1, 
               0, 1, 1, 1, 0, 0, 
               0, 1, 1, 1, 1, 1, 
               0, 1, 1, 0, 1, 1, 
               0, 0, 0, 1, 1, 1, 
               1, 1, 1, 1, 0, 1), ncol=6, byrow=TRUE), 9)

    f(matrix(c(0, 0, 0, 0, 0, 0, 
               0, 0, 0, 1, 0, 1, 
               0, 0, 1, 1, 1, 1, 
               0, 1, 1, 1, 1, 1, 
               1, 1, 1, 1, 1, 1, 
               0, 1, 0, 0, 0, 1), ncol=6, byrow=TRUE), 0)

    f(matrix(c(0, 0, 0, 0, 0, 0, 
               0, 0, 0, 1, 0, 1, 
               0, 0, 1, 0, 1, 0, 
               0, 1, 1, 0, 0, 0, 
               1, 1, 1, 1, 1, 0, 
               0, 1, 0, 0, 0, 1), ncol=6, byrow=TRUE), 6)

    f(matrix(c(1, 0, 1, 0, 1, 0, 
               0, 1, 0, 1, 0, 1, 
               1, 0, 1, 0, 1, 0, 
               0, 1, 0, 1, 0, 1, 
               1, 0, 1, 0, 1, 0, 
               0, 1, 0, 1, 0, 1), ncol=6, byrow=TRUE), 12)

    f(matrix(c(0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 1), ncol=6, byrow=TRUE), 0)

    f(PredationMatrix(TL84), 424)

    # Test the R function
    AssertEqual(SumDietGaps(TL84), 424)
}

TestMinimiseSumDietGaps <- function()
{
    # Can't test for exact values because simulated annealing learning is 
    # stochastic and not guaranteed to reach a global minimum
    to.test <- list(list(c2, 0), 
                    list(c3, 0),
                    list(c4, 0),
                    list(c5, 0),
                    list(c6, 0),
                    list(TL84, 10), 
                    list(TL86, 10), 
                    list(YthanEstuary, 277),
                    list(SkipwithPond, 31),
                    list(BroadstoneStream, 7),
                    list(Benguela, 27))
    for(index in 1:length(to.test))
    {
        community <- to.test[[index]][[1]]
        expected <- to.test[[index]][[2]]
        print(community)

        res <- MinimiseSumDietGaps(community, T.start=1, T.stop=0.1, c=0.99, 
                                   swaps.per.T=1000, n=3)

        if(0==expected)
        {
            AssertEqual(0, res$sum.gaps)
        }
        else
        {
            if( (res$sum.gaps / expected) > 3)
            {
                stop(paste('Best diet gap for', community, 'of', res$sum.gaps,
                           'is much higher than the expected value of', 
                            expected))
            }
        }
    }
}

TestTrophicSimilarity <- function()
{
    # Helper function that creates expected matrices
    F <- function(data, nodes)
    {
        m <- matrix(data, ncol=length(nodes), dimnames=list(nodes, nodes))
        mode(m) <- 'numeric' # Necessary for a 1 x 1 matrix that contains only NA
        return (m)
    }
    AssertEqual(F(NA, 'S'), TrophicSimilarity(c1))
    AssertEqual(F(NA, 'S'), TrophicSimilarity(c2))
    AssertEqual(F(c(1,0,0,1), c('R','C')), TrophicSimilarity(c3))
    AssertEqual(F(c(1,0,1,0,1,0,1,0,1), c('R','C','P')), TrophicSimilarity(c4))
    AssertEqual(F(c(1,1/3,1/3,1/3,1,1/3,1/3,1/3,1), c('R','C','O')), 
                TrophicSimilarity(c5))
    AssertEqual(F(c(1,0,1,0,1,0,1,0,1), c('R','C','P')), TrophicSimilarity(c6))
    AssertEqual(F(c(1,1/3,1/3,0,0, 1/3,1,1/3,0,0, 1/3,1/3,1,0,0, 0,0,0,NA,NA, 0,0,0,NA,NA), 
                  LETTERS[1:5]), TrophicSimilarity(c7))
}

TestOmnivores <- function()
{
    AssertEqual(vector(mode='character'), Omnivores(c1))
    AssertEqual(vector(mode='character'), Omnivores(c2))
    AssertEqual(vector(mode='character'), Omnivores(c3))
    AssertEqual(vector(mode='character'), Omnivores(c4))
    AssertEqual('O', Omnivores(c5))
    AssertEqual(vector(mode='character'), Omnivores(c6))
    AssertEqual(vector(mode='character'), Omnivores(c7))

    AssertEqual(c('Microzooplankton', 'Zooplankton', 'Ctenophores', 
        'Chrysaora quinquecirrha', 'Other suspension feeders', 'Mya', 
        'Crassostrea virginica', 'Callinectes sapidus', 'Anchoa mitchilli', 
        'Brevoortia tyrannus', 'Trinectes maculatus', 'Leiostomus xanthurus', 
        'Morone americana', 'Pomatomus saltatrix', 'Paralichthys dentatus', 
        'Morone saxatilis'), Omnivores(ChesapeakeBay))
    AssertEqual(c('Microzooplankton', 'Zooplankton', 'Ctenophores', 
        'Chrysaora quinquecirrha', 'Other suspension feeders', 'Mya', 
        'Crassostrea virginica', 'Callinectes sapidus', 'Anchoa mitchilli', 
        'Brevoortia tyrannus', 'Trinectes maculatus', 'Leiostomus xanthurus', 
        'Pomatomus saltatrix', 'Paralichthys dentatus', 'Morone saxatilis'), 
      Omnivores(ChesapeakeBay, level=ChainAveragedTrophicLevel))
}

TestOmnivory <- function()
{
    AssertEqual(0, Omnivory(c1))
    AssertEqual(0, Omnivory(c2))
    AssertEqual(0, Omnivory(c3))
    AssertEqual(0, Omnivory(c4))
    AssertEqual(1/3, Omnivory(c5))
    AssertEqual(0, Omnivory(c6))
    AssertEqual(0, Omnivory(c7))
    # Value of 0.455 for Chesapeake Bay presented by Bersier et al 2002 Ecology
    res <- Omnivory(ChesapeakeBay, level=ChainAveragedTrophicLevel)
    AssertEqual(0.455, round(res, 3))
}

TestNodeQuantitativeDescriptors <- function()
{
    res <- NodeQuantitativeDescriptors(ChesapeakeBay, 'biomass.flow')

    # The values in Bersier et al 2002, Table 1:
    expected <- matrix(c(
0, 7,      0, 80051, 0,    3.19, 0,    0,    0,    0,    0,    0,    2.62, 5.47,
0, 8,      0,  2977, 0,    3.21, 0,    0,    0,    0,    0,    0,    2.64, 0.21,
0, 6,      0,294955, 0,    3.63, 0,    0,    0,    0,    0,    0,    2.99,22.96,
0, 1,      0, 18086, 0,    1,    0,    0,    0,    0,    0,    0,    0.82, 0.39,
0, 1,      0, 88721, 0,    1,    0,    0,    0,    0,    0,    0,    0.82, 1.9,
1, 1,  88721, 31638, 1,    1,    0.5,  0.74, 0,    0,    0.65, 3.9,  0.82, 0.68,
3, 5,  64224, 11742, 2.13, 2.42, 0.47, 0.83, 0.89, 1,    1.39, 6,    2,    0.61,
3, 7,  46389,  9855, 1.81, 2.48, 0.42, 0.77, 0.89, 0.56, 1.19, 3.69, 2.04, 0.52,
3, 1,  10447,   552, 2.01, 1,    0.67, 0.97, 2,    1.01, 1.31, 0.92, 0.82, 0.01,
2, 0,   1711,     0, 1.88, 0,    1,    1,    1,    0.88, 1.23, 0.14, 0,    0,
3, 1,   4594,   538, 1.41, 1,    0.58, 0.92, 0.89, 0.27, 0.92, 0.28, 0.82, 0.01,
3, 2,   2488,   224, 1.41, 1.18, 0.54, 0.93, 0.89, 0.26, 0.92, 0.15, 0.97, 0.01,
3, 0,   4830,     0, 1.41, 0,    1,    1,    0.89, 0.26, 0.92, 0.3,  0,    0,
1, 5, 160831,   609, 1,    3.41, 0.23, 0.99, 0,    0,    0.65, 7.06, 2.8,  0.04,
1, 6,  25062,   547, 1,    2.94, 0.25, 0.94, 0,    0,    0.65, 1.1,  2.42, 0.03,
1, 2,  57335,  4593, 1,    1.07, 0.48, 0.92, 0,    0,    0.65, 2.52, 0.88, 0.11,
2, 0,  54048,     0, 1.89, 0,    1,    1,    0,    0,    1.24, 4.49, 0,    0,
1, 6,  14075,  1027, 1,    1.31, 0.43, 0.91, 0,    0,    0.65, 0.62, 1.08, 0.03,
7, 2,   8593,   295, 4.12, 1.05, 0.8,  0.99, 2.59, 1.48, 2.69, 1.55, 0.86, 0.03,
1, 0,      5,     0, 1,    0,    1,    1,    0,    0,    0.65, 0,    0,    0,
1, 1,     26,     0, 1,   1,    0.5,  0.99, 0,    0,    0.65, 0,    0.82, 0,
3, 6,   1823,   153, 1.59, 3.58, 0.31, 0.84, 0.89, 0.55, 1.04, 0.13, 2.94, 0.1,
3, 3,    273,    21, 1.42, 2.64, 0.35, 0.87, 0.89, 0.36, 0.93, 0.02, 2.17, 0,
1, 0,      5,     0, 1,    0,    1,    1,    0,    0,    0.65, 0,    0,    0,
3, 0,     10,     0, 1.93, 0,    1,    1,    0,    0,    1.26, 0,    0,    0,
4, 0,     96,     0, 2.95, 0,    1,    1,    0.75, 0.36, 1.93, 0.01, 0,    0,
5, 1,    483,    10, 2.66, 1,    0.72, 0.99, 0.65, 0.14, 1.74, 0.06, 0.82, 0,
3, 0,    150,     0, 2.56, 0,    1,    1,    0.89, 0.38, 1.67, 0.02, 0,    0,
3, 0,    219,     0, 2.24, 0,    1,    1,    0,    0,    1.46, 0.02, 0,    0,
3, 0,     16,     0, 2.41, 0,    1,    1,    0.89, 0.90, 1.58, 0,    0,    0,
1, 1,     91,     4, 1,    1,    0.5,  0.96, 0,    0,    0.65, 0,    0.82, 0,
4, 0,     26,     0, 3.11, 0,    1,    1,    1.83, 0.82, 2.03, 0,    0,    0,
4, 0,     30,     0, 2.52, 0,    1,    1,    0,    0, 1.65, 0,       0,    0
), ncol=14, byrow=TRUE)
    colnames(expected) <- colnames(res)
    rownames(expected) <- rownames(res)

    # Some differences to the published values
    # TODO Contact Bersier to get to the bottom of these differences
    expected['Ctenophores','bIn'] <- 10446          # Change from 10447
    expected['Macoma spp','bIn'] <- 57334           # Change from 57335
    expected['Leiostomus xanthurus','bIn'] <- 482   # Change from 483
    expected['Arius felis','bIn'] <- 218            # Change from 219
    expected['Microzooplankton','nP'] <- 2.43       # Change from 2.42
    expected['Leiostomus xanthurus','d.prime'] <- 0.73 # Change from 0.72
    expected['Callinectes sapidus','v'] <- 0.01     # Change from 0.03
    expected['Anchoa mitchilli','v'] <- 0.01        # Change from 0.1

    # bIn and bOut are rounded to zero dp in the paper
    res[,'bIn'] <- round(res[,'bIn'], 0)
    res[,'bOut'] <- round(res[,'bOut'], 0)
    AssertEqual(round(res, 2), expected)
}

TestQuantitativeDescriptors <- function()
{
    res <- QuantitativeDescriptors(ChesapeakeBay, 'biomass.flow')
    expected <- matrix(c(
        0.364,   0.364,   0.455,
        0.485,   0.485,   0.394,
        0.152,   0.152,   0.152,
        0.75,    0.711,   0.317,    # Col 3 shown as 0.316 in the paper
        2.212,   1.372,   2.096,
        0.067,   0.042,   0.064,
        0.397,   0.451,   0.011,
        0.055,   0.055,   0.195,
        0.288,   0.287,   0.112,
        0.260,   0.207,   0.682,
        3.157,   1.605,   1.517,
        3.000,   1.618,   0.675,
        1.096,   0.488,   1.910,
        6.000,   2.970,   6.408,
        0.455,   0.510,   0.280, 
        2.607,   1.801,   1.375,
        3.476,   1.911,   2.818,    # Col 3 shown as 2.812 in the paper
        0.739,   0.645,   1.882,
        1.195,   1.032,   4.066
        ), ncol=3, byrow=TRUE)
    colnames(expected) <- colnames(res)
    rownames(expected) <- rownames(res)
    AssertEqual(round(res, 3), expected)
}
