ReadInstr = function(fn){
  OUT = list()
  for(i in 1:length(fn)){ # loop through files
    f = file(fn, 'r')
    RESP = list()
    PZ = 0 # specifies whether we're in the poles/zeros part of the file
    while(length({line = readLines(f, 1)}) > 0){ # loop until we read an empty line (end of file)

      # break line into character vector with whitespace as separator
      line = strsplit(line, '[[:space:]]')[[1]]
      line = line[line != ''] # remove empty elements of character vector
      
      # check for sensitivity line
      if(!is.na(charmatch('SENSITIVITY', line))){
        RESP$Sense = as.numeric(line[4])
      }
      # check for constant line
      if(!is.na(charmatch('A0', line))){
        RESP$Knorm = as.numeric(line[4])
      }
      # check for start of zeros
      if(!is.na(charmatch('ZEROS', line))){
        RESP$nz = as.numeric(line[2])
        PZ = 1
        next
      }
      # check for start of poles
      if(!is.na(charmatch('POLES', line))){
        RESP$np = as.numeric(line[2])
        PZ = 2
        next
      }
      # check for end of poles--a line saying CONSTANT is after the end of the poles section
      if(!is.na(charmatch('CONSTANT', line))){
        PZ = 0
      }
      # read in zeros, if in zeros section
      if(PZ == 1){
        line = as.numeric(line)
        RESP$zeros = c(RESP$zeros, line[1] + 1i * line[2])
      }
      # read in poles, if in poles section
      if(PZ == 2){
        line = as.numeric(line)
        RESP$poles = c(RESP$poles, line[1] + 1i * line[2])
      }
    }

    # We currently have a response of counts per m displacement.  To fit other
    # functions in RSEIS, we want this as counts for velocity.  So, we remove
    # a zero at the origin.
    w = which(RESP$zeros == 0)[1]
    RESP$zeros = RESP$zeros[-w]
    OUT[[i]] = RESP
    close(f)
  }
  return(OUT)
}
