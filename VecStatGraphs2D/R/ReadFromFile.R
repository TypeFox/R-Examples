ReadFromFile <- function (input) 
{
    vector = try(read.table(input, header = FALSE, dec = ".", 
        ))
    if (class(vector) == "try-error") {
        print("Error in file!")
        return(0)
    }
    return(vector)
}
