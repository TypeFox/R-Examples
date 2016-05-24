## Put data into an unattached package:
g.data.put <- function(item, value, dir) {
    assign(item, value)
    save(list=item, file=g.data.mash(dir, item))
}
