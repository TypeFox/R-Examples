ReadBNFFile <- function(filename) {
  # reads a bnf grammar file and returns a list structure
  
  # read the file line by line
  con=file(filename, open="r")
  lines=readLines(con) 
  close(con)
  
  # parse the lines
  rule_list = list()
  for (l in lines) {
    l = trim_space(l)
    gram_line = strsplit(l, "::=")[[1]]
    if (length(gram_line) > 1) {
      # split and trim rules
      rules = strsplit(trim_space(gram_line[2]), "|", fixed=TRUE)
      for (j in seq_along(rules[[1]])) {
        rules[[1]][[j]] = trim_space(rules[[1]][[j]])
      }
      
      # add rules to list
      i = length(rule_list) + 1
      rule_list[[i]] = list()
      rule_list[[i]][[1]] = trim_space(gram_line[[1]])
      rule_list[[i]][[2]] = as.list(rules[[1]])
    }
  }
  
  return (rule_list)
}

