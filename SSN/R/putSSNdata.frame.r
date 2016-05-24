putSSNdata.frame <- function(DataFrame, x, Name = "Obs")
{
    if(!is.data.frame(DataFrame)) {
        cat("Incoming DataFrame is not actually a data.frame\n")
        return(x)
    }

    if (!"pid" %in% names(DataFrame)| !all(row.names(DataFrame) == DataFrame$pid)) {
        stop("data.frame 'pid' column is missing or does not match the data.frame row.names\n")
        return(x)
    }

    if(Name == "Obs") {
        if(class(x) == "SpatialStreamNetwork"){
            if(!all(DataFrame$pid == x@obspoints@SSNPoints[[1]]@point.data$pid)) {
                stop("Input data.frame is incompatible with existing SSN structure\n")
                return(x)
            }
            x@obspoints@SSNPoints[[1]]@point.data <- DataFrame
            return(x)
        }
        if((class(x) == "glmssn.predict") |
           (class(x) == "influenceSSN") |
           (class(x) == "glmssn"))
        {
            if(!all(DataFrame$pid == x$ssn.object@obspoints@SSNPoints[[1]]@point.data$pid)) {
                stop("Input data.frame is incompatible with existing SSN structure\n")
                return(x)
            }
            x$ssn.object@obspoints@SSNPoints[[1]]@point.data <- DataFrame
            return(x)
        }
    }
    if(Name != "Obs") {
        if(class(x) == "SpatialStreamNetwork"){
            if(!all(DataFrame$pid == x@predpoints@SSNPoints[x@predpoints@ID == Name][[1]]@point.data$pid)) {
                stop("Input data.frame is incompatible with existing SSN structure\n")
                return(x)
            }
            x@predpoints@SSNPoints[x@predpoints@ID == Name][[1]]@point.data <- DataFrame
            return(x)
        }
        if((class(x) == "glmssn.predict") |
           (class(x) == "influenceSSN") |
           (class(x) == "glmssn"))
        {
            if(!all(DataFrame$pid == x$ssn.object@predpoints@SSNPoints[x$ssn.object@predpoints@ID == Name][[1]]@point.data$pid)) {
                stop("Input data.frame is incompatible with existing SSN structure\n")
                return(x)
            }
            x$ssn.object@predpoints@SSNPoints[x$ssn.object@predpoints@ID == Name][[1]]@point.data <- DataFrame
            return(x)
        }
    }

}
