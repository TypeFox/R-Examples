`create.ID.col` <-
function(x,sample.id=c("sample_name","sample_treat")){

    # get the columns needed to create the id
    colIndex <- colnames(x[[4]]) %in% sample.id

    # paste the appropriate columns linewise together
    # us as.matrix for the accession via colIndex, in case there is only one column of the 
    # sample description it would return a vector rather than a matrix, so we have to coerce it
    identifier <- as.character(apply(as.matrix(x[[4]][,colIndex], nrow=nrow(x[[4]])), 1, paste, collapse=""))

    # add the new identifier vector to the sampledescription
    x[[4]] <- cbind(x[[4]],data.frame(identifier=identifier))

    return(x)
}

