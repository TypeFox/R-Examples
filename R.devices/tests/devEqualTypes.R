library("R.devices")
devEqualTypes <- R.devices:::.devEqualTypes
png <- grDevices::png
postscript <- grDevices::postscript

message("*** devEqualTypes() ...")

message("*** devEqualTypes('png', 'png') ...")
res <- devEqualTypes("png", "png")
stopifnot(res)

message("*** devEqualTypes('png', png) ...")
res <- devEqualTypes("png", png)
stopifnot(res)

message("*** devEqualTypes(png, 'png') ...")
res <- devEqualTypes(png, "png")
stopifnot(res)

message("*** devEqualTypes(foo, png) ...")
foo <- png
res <- devEqualTypes(foo, png)
stopifnot(res)

message("*** devEqualTypes('png', postscript) ...")
res <- devEqualTypes("png", postscript)
stopifnot(!res)

message("*** devEqualTypes(postscript, 'png') ...")
res <- devEqualTypes(postscript, "png")
stopifnot(!res)

message("*** devEqualTypes('non-existing', png) ...")
res <- devEqualTypes("non-existing", png)
stopifnot(!res)

message("*** devEqualTypes(png, 'non-existing') ...")
res <- devEqualTypes(png, "non-existing")
stopifnot(!res)

message("*** devEqualTypes() ... DONE")
