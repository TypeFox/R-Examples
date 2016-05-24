## model wrapper
binary_model <- function(command) {
    invoke <- function(param) {
        write.table(param, file = "input", row.names = F, col.names = F, quote = F)
        system(command)
        file.remove("input")
        output = read.table("output", header = F)
        file.remove("output")
        output
    }
} 
