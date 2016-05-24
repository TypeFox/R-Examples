## Get data from an unattached package:
g.data.get <- function(item, dir) get(load(g.data.mash(dir, item)))
