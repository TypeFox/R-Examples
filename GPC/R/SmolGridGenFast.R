SmolGridGenFast <- function(d,ThetaNodes,SparseGridTmp,CoordDummy,JListDummy) { # this function changes only [NodeCoord counter JList] 
    if (d == 1){
        for (i in 1:(sum(!is.na(ThetaNodes[d,])))){
            JListDummy[d] = i
            CoordDummy[d] = ThetaNodes[d,i]
            SparseGridTmp$n = SparseGridTmp$n + 1
            # cat(JListDummy)
            SparseGridTmp$JList[SparseGridTmp$n,] = JListDummy            
            SparseGridTmp$NodeCoord[SparseGridTmp$n,] = CoordDummy
        }
    } else {
        for (i in 1:(sum(!is.na(ThetaNodes[d,])))){
            JListDummy[d] = i;
            CoordDummy[d] = ThetaNodes[d,i];
            SparseGridTmp = SmolGridGenFast(d-1,ThetaNodes,SparseGridTmp,CoordDummy,JListDummy)      
        }
    }
    return(SparseGridTmp)
}