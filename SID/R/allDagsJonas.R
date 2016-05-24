allDagsJonas <- 
function (adj,row.names)
{
    # Input: adj. mat of a DAG with row.names, containing the undirected component that 
    #        should be extended
    # !!!! the function can probably be faster if we use partial orderings 
    
    a <- adj[row.names,row.names]
    
    if(any((a + t(a))==1))
    {
        #warning("The matrix is not entirely undirected.")
        return(-1)
    }
    return(allDagsIntern(adj, a, row.names,NULL))
}

