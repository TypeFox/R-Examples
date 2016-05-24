msgpack.writeResult <-
function(filename, result) {
	fl <- file(filename, "wb")
	writeBin(result, fl)
	close(fl)
}
