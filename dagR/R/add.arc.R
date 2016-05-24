add.arc <-
function(dag, arc, type=0)
{ # this adds arcs;
  # type=0 (directed; default); type=1 (association);
  # it removes those parts of the DAG that need to be re-obtained;
  # CAVE: it uses the node numbering of the DAG, not those used in dag.init;
    dag$arc <- rbind(dag$arc, arc)
    dag$arc.type <- c(dag$arc.type, type)
    dag$paths <- NULL
    dag$pathsN <- NULL
    dag$path.status <- NULL
    dag$searchType <- NULL
    dag$searchRes <- NULL
    return(dag)
}
