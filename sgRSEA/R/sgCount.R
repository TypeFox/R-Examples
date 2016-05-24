sgCount <-
function(fqInfo, annotation, seqlen=20, output_prefix='sgCount'){
    path <-paste(system.file(package="sgRSEA"), "sgCount_R.py", sep="/")
    command <- paste("python", path, fqInfo, annotation, seqlen, output_prefix)
    response <- system(command,intern=T)
    print(response)
    }
