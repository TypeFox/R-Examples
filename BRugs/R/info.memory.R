"infoMemory" <-
    function(){
        command <- "BugsEmbed.AllocatedMemory"
        res <- .CmdInterpreter(command)
        buffer <- file.path(tempdir(), "buffer.txt")
        res <- readLines(buffer)
        mem <- as.numeric(gsub("^([0-9]+).+", "\\1", res))
        mem
    }
