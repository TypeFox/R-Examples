library(gdata)
saveto <- tempfile(pattern = "test.txt", tmpdir = tempdir())

write.fwf(x = data.frame(a=1:length(LETTERS), b=LETTERS),
          file=saveto, eol="\r\n")
