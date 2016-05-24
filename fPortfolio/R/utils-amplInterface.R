
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  amplModelOpen         Opens a writes to an AMPL model file
#  amplModelAdd          Adds model specs to an existing AMPL model file
#  amplModelShow         Shows the content of an AMPL .mod file
# FUNCTION:             DESCRIPTION:
#  amplDataOpen          Opens and writes the header to an AMPL data file
#   amplDataAddValue      Adds a numeric value to an AMPL data file
#   amplDataAddVector     Adds a numeric vector to an AMPL data file
#   amplDataAddMatrix     Adds a numeric matrix to an AMPL data file
#   amplDataSemicolon     Adds a semicolon on the end of a data input line
#  amplDataShow           Shows the content of an AMPL data file
# FUNCTION:             DESCRIPTION:
#  amplRunOpen           Opens a run file
#  amplRunAdd            Adds run specs to an existing AMPL run file
#  amplRunShow           Shows the content of an AMPL run file   
# FUNCTION:             DESCRIPTION:
#  amplOutShow           Shows the content of an AMPL output txt file     
################################################################################


amplModelOpen <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Writes an AMPL model to a .mod file
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "mod", sep = ".")
    
    # Write Model to File:
    write("", file=file, ncolumns=1, append=FALSE)   

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplModelAdd <-
    function(model, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Writes an AMPL model to a .mod file
    
    # Arguments:
    #   model - a character vector with the lines making the model file
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "mod", sep = ".")
    
    # Write Model to File:
    write(model, file=file, ncolumns=1, append=TRUE)   

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplModelShow <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Prints an AMPL .mod file
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "mod", sep = ".")
    
    # Trace File:
    cat(readLines(file), sep = "\n")

    # Return Value:
    invisible()
}


################################################################################


amplDataOpen <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Writes the header to an AMPL data file
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "dat", sep = ".")
    
    # Write to the file:
    write("", file=file, ncolumns=1, append=FALSE)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplDataAdd <- 
    function(name, data, type, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Check data type and adds to AMPL data file
    
    # Arguments:
    #   name - AMPL data name
    #   data - the data object, a numeric value, vector or matrix
    #   type - eiher "value", "vector" or "matrix"
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:

    # Add Data:
    if (type == "value") {
        amplDataAddValue(data=name, value=data, project=project) }
    if (type == "vector") {
        amplDataAddVector(data=name, vector=data, project=project) }
    if (type == "matrix") {
        amplDataAddMatrix(data=name, matrix=data, project=project) }
        
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplDataShow <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Prints an AMPL data file
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "dat", sep = ".")
    
    # Trace File:
    cat(readLines(file), sep = "\n")

    # Return Value:
    invisible()
}


# -----------------------------------------------------------------------------


amplDataSemicolon <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Adds a semicolon on the end of a data input line
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "dat", sep = ".")
    
    # Add a Semicolon and an Empty Line:
    write(";", file=file, ncolumns=1, append=TRUE)
    write(" ", file=file, ncolumns=1, append=TRUE)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplDataAddValue <-
    function(data, value, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Adds a numeric value to an AMPL data file
    
    # Arguments:
    #   data - a character string, the name of the value
    #   value - a numeric value, the value of the numeric input variable
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    name <- data
    
    # Compose File Name:
    file <- paste(project, "dat", sep = ".")
    
    # Write Name:
    x <- paste("param", name, ":=", collapse=" ")
    write(x, file=file, ncolumns=1, append=TRUE)
    
    # Write Numeric Value:
    write(value, file=file, ncolumns=1, append=TRUE)
    amplDataSemicolon(project) 
    
    # Return Value:
    invisible()
}    


# ------------------------------------------------------------------------------


amplDataAddVector <-
    function(data, vector, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Adds a numeric vector to an AMPL data file
    
    # Arguments:
    #   data - a character string, the name of the vector
    #   value - a numeric vector, the values of the numeric input vector
    #   project -  a project name, gives the root name of the model file
    
    # FUNCTION:
    
    name <- data
    
    # Compose File Name:
    file = paste(project, "dat", sep = ".")
    
    # Write Name:
    x <- paste("param", name, ":=", collapse=" ") 
    write(x, file=file, ncolumns=1, append=TRUE)
    
    # Write Vector:
    vector <- as.vector(vector)
    N <- length(vector)
    write(t(cbind(1:N, vector)), file=file, ncolumns=2, append=TRUE)
    amplDataSemicolon(project) 
    
    # Return Value:
    invisible()
}   


# ------------------------------------------------------------------------------


amplDataAddMatrix <-
    function(data, matrix, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Adds a numeric matrix to an AMPL data file
    
    # Arguments:
    #   data - a character string, the name of the matrix
    #   value - a numeric matrix, the values of the numeric input matrix
    #   project -  a project name, gives the root name of the model file
    
    # FUNCTION:
    
    name <- data
    
    # Compose File Name:
    file <- paste(project, "dat", sep = ".")
    
    # Write Name:
    x <- paste("param", name, ":", collapse=" ")
    write(x, file=file, ncolumns=1, append=TRUE)
    
    # Write Matrix:
    N <- ncol(matrix)
    x <- paste(paste(1:N, collapse = " "), ":=")
    write(x, file = file, ncolumns = 1, append = TRUE)
    X <- cbind(1:nrow(matrix), matrix)
    colnames(X) = rownames(X) = NULL
    write(t(X), file = file, ncolumns = 1, append = TRUE, sep = " ")
    amplDataSemicolon(project) 
    
    # Return Value:
    invisible()
}


################################################################################


amplRunOpen <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Writes an AMPL run to a .run file   
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "run", sep = ".")
    
    # Open:
    write("", file=file, ncolumns=1, append=FALSE) 
    
    # Return Value:
    invisible()
}


# -----------------------------------------------------------------------------


amplRunAdd <-
    function(run, project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Writes an AMPL run to a .run file   
    
    # Arguments:
    #   run - a character vector with the lines making the model file
    #   project -  a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "run", sep = ".")
    
    # Write Run Commands to File:
    write(run, file=file, ncolumns=1, append=FALSE) 
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


amplRunShow <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Prints an AMPL .run file
    
    # Arguments:
    #   project -  a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file = paste(project, "run", sep = ".")
    
    # Trace File:
    cat(readLines(file), sep = "\n")

    # Return Value:
    invisible()
}


################################################################################


amplOutShow <-
    function(project)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Shows the content of an AMPL output txt file     
    
    # Arguments:
    #   project - a project name, gives the root name of the model file
    
    # FUNCTION:
    
    # Compose File Name:
    file <- paste(project, "txt", sep = ".")
    
    # Trace File:
    cat(readLines(file), sep = "\n")

    # Return Value:
    invisible()
}


################################################################################

