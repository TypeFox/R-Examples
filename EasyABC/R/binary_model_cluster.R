## model wrapper
binary_model_cluster <- function(command) {
    invoke <- function(param) {
        numclust = param[1]
        nparam = length(param)
        input_file_name = paste("input", numclust, sep = "")
        output_file_name = paste("output", numclust, sep = "")
        write.table(param, file = input_file_name, row.names = F, col.names = F, 
            quote = F)
        system2(command, args = as.character(numclust))
        file.remove(input_file_name)
        output = read.table(output_file_name, header = F)
        file.remove(output_file_name)
        output
    }
} 
