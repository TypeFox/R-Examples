## Create potential list (rip, universe)
##
.createPotList <- function(rip.order, universe)
{
    cliques     <- rip.order$cliques
    APlist  <- as.list(rep(NA,length(cliques)))
    for ( i  in 1:length(cliques))
        {
            cq    <- cliques[[ i ]]
            vlab  <- universe$levels[cq]
            APlist[[ i ]] <- parray(cq, vlab)
        }
    APlist
}

## Insert cpt's into potential list (cptlist, APlist)
##
.insertCPT <- function(cptlist, APlist, details=0)
{
    if (details>=1) cat(".Inserting cpt's in potential list [.insertCPT]\n")

    APnames <- lapply(APlist, function(x) names(dimnames(x)))
    CPnames <- unname(lapply(cptlist, function(x) varNames(x)))

    ## FIXME .findHosts can be replaced by gRbase::get_superset_
    hosts    <- .findHosts( CPnames, APnames )

    for ( i  in 1:length(cptlist))
        {
            cptc <- cptlist[[ i ]]
            j    <- hosts[ i ]
            APlist[[ j ]] <- tableOp( APlist[[ j ]], cptc, "*" )
        }
    .infoPrint(details, 4, {cat("....APlist (after insertion):\n"); print(APlist) })
    APlist
}

.insertNA <- function(list.of.tables)
{
    lapply(list.of.tables,
           function(xxx)
           {
               xxx[] <- NA
               xxx
           } )
}


## FIXME .findHosts can be replaced by gRbase::get_superset_
.findHosts <- function( xx, yy ){
  unlist(lapply(1:length(xx), function( i ) which(isin(yy, xx[[  i  ]], index=T)>0)[1]))
}














