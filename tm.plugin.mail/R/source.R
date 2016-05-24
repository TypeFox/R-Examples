MBoxSource <-
function(mbox, encoding = "")
{
    f <- file(mbox)
    open(f)
    on.exit(close(f))
    message.nr <- 0
    offsets <- integer(0)
    lines <- integer(0)
    while (length(line <- readLines(f, 1)) == 1) {
        if (grepl("^From ", line, useBytes = TRUE)) {
            message.nr <- message.nr + 1
            offsets[message.nr] <- seek(f)
            lines[message.nr] <- 0
        } else
            lines[message.nr] <- lines[message.nr] + 1
    }
    SimpleSource(encoding = encoding, length = length(lines), reader = readMail,
                 mbox = mbox, msgOffsets = offsets, msgLines = lines,
                 class = "MBoxSource")
}

getElem.MBoxSource <-
function(x)
{
    f <- file(x$mbox)
    open(f)
    on.exit(close(f))
    seek(f, x$msgOffsets[x$position])
    list(content = iconv(readLines(f, x$msgLines[x$position]),
                         x$encoding, "UTF-8", "byte"),
         uri = x$mbox)
}
