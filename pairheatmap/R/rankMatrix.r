############################################
##  R package pairheatmap
##  fucntion: rankMatrix.r
############################################


## rank functions
rankMatrix <- function(m1, m2, group=NULL,
                       clustMethod="complete", clustMembers=NULL,
                        orderRowGroup, clusterRow, clusterCol, clusterColTogether
                       )
{
    # check same row
    if (nrow(m1) != nrow(m2)) stop("matrix 1 must have same row as matrix 2")

    rOrder <- c(1:nrow(m1))
    cOrder <- c(1:ncol(m1))
    cOrder2 <- c(1:ncol(m2))

    if (is.null(group))
    {

        clusterResult.m1 <- clustGroup(m1,
                         clusterRow, clusterCol, clustMethod, clustMembers,
                         rOrder, cOrder)
        clusterResult.m2 <- clustGroup(m2,
                         clusterRow, clusterCol, clustMethod, clustMembers,
                         rOrder, cOrder2)
        rowInd=clusterResult.m1$rowInd
        colInd.m1 =clusterResult.m1$colInd
        colInd.m2 =clusterResult.m2$colInd

    }
    else 
    {
        if (!is.vector(group) || length(group)!=nrow(m1))
            stop("group length must be equal to the row of matrix 1 and 2!")
        
        group.unique <- unique(group)

        if (!is.null(orderRowGroup)) group <- factor(group, levels=orderRowGroup, ordered=TRUE)
        else  group <- factor(group, levels=group.unique, ordered=TRUE)

        if (clusterRow)
        {
            rowResult <- rankGroup(m1, m2, group,
                              clusterRow, clusterCol, clustMethod, clustMembers)
            rowInd <- rowResult$rowInd
        }
        else rowInd <- order(group)

        if (clusterCol)
        {
            clusterResult.m1 <- clustGroup(m1, "FALSE", clusterCol,
                         clustMethod, clustMembers,
                         rowInd, cOrder)
            clusterResult.m2 <- clustGroup(m2, "FALSE", clusterCol,
                         clustMethod, clustMembers,
                         rowInd, cOrder2)
            colInd.m1 <- clusterResult.m1$colInd
            colInd.m2 <- clusterResult.m2$colInd
        }
        else
        {
            colInd.m1 <- cOrder
            colInd.m2 <- cOrder2
        }

        # control col relation between two heatmaps
        if (clusterColTogether)
        {
            if (ncol(m1) == ncol(m2))
            {
               colInd.m2 <- colInd.m1
            }
            else if  (ncol(m1) > ncol(m2))
            {
                done.ind <- which(colInd.m1 <= ncol(m2))
                colInd.m2 <- colInd.m1[done.ind]
            }
            else
            {
                colInd.m2 <- c(colInd.m1, (ncol(m1)+1) : ncol(m2))
            }
        }
        
    }

   return(list(rowInd=rowInd, colInd.m1=colInd.m1, colInd.m2=colInd.m2))
}


# rankGroupSimilarity(m1, m2, group, FALSE)
#rankGroupSimilarity <- function(m1, m2, group, rank)
#{
#   group.unique <- unique(group)
#   score <- NULL
#   for (i in group.unique)
#   {
#      ind <- which(group==i)
#      score <- c(score, getGroupSimScore(m1[ind,], m2[ind,]))
#   }
   
#   if (rank) group.unique <- group.unique[order(score, decreasing=FALSE)]
#   group <- factor(group, levels=group.unique, ordered=TRUE)
   
#   m1 <- m1[order(group),]
#   m2 <- m2[order(group),]
#   indOrder <- order(group)
#   groupOrder <- group[order(group)]
   
#   return(list(m1=m1, m2=m2, indOrder=indOrder, groupOrder=groupOrder))
#}

clustGroup <- function(m1,
                         clusterRow, clusterCol, method, members,
                         rOrder,
                         cOrder
                         )
{
        if (clusterCol)
        {
            hc.col = hclust(dist(t(m1)), method=method, members=members)
            cOrder <- cOrder[hc.col$order]
        }

        if (clusterRow)
        {
            hc.row = hclust(dist(m1), method=method, members=members)
            rOrder <- rOrder[hc.row$order]
        }
        
        return(list(rowInd=rOrder, colInd =cOrder ))

} 

## cluster row
rankGroup <- function(m1, m2, rawGroup,
                         clusterRow, clusterCol, clustMethod, clustMembers)
{

    rowgroup <- rawGroup[order(rawGroup)]

    ## cluster sparately
    rowgroup.unique <- rev(unique(rowgroup))
    new.rowInd <- NULL

    for (i in 1:length(rowgroup.unique))
    {
        select.ind <- which(rawGroup==rowgroup.unique[i])
        m1.g <- m1[select.ind,]
        m2.g <- m2[select.ind,]
        
        if (length(select.ind)>2)
        {
            clusterResult <- clustGroup(m1.g,
                             clusterRow, clusterCol, clustMethod, clustMembers,
                             select.ind, 1:ncol(m1.g))
            new.rowInd <- c(new.rowInd, clusterResult$rowInd)

        }
        else
        {
            new.rowInd <- c(new.rowInd, select.ind)
        }
        
    }

    return(list(rowInd=new.rowInd, colInd=1:ncol(m1)))
    
}

