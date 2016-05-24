copy.magpie <- function(input_file,output_file) {
  write.magpie(read.magpie(input_file),output_file)
}