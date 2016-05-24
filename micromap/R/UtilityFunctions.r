if(getRversion() >= "2.15.1")
   utils::globalVariables(c("color", "coordsx", "coordsy", "fill.color",
      "hole", "IDpoly", "pGrp", "pGrpOrd", "plot.height", "plot.width",
      "print.file", "region", "tmp.adj", "tmp.data", "tmp.data1", "tmp.data2",
      "tmp.data3", "tmp.data4", "tmp.data5", "tmp.data6", "tmp.data7",
      "tmp.label", "tmp.labels", "tmp.y", "textx", "texty", "xmax", "xmin",
      "ymax", "ymin"))


all_atts <- function(a, att.name)
   unlist(unlist(a)[which(att.name==names(unlist(a)))])


all_attsb <- function(a, att.name)
   as.logical(unlist(a)[which(att.name==names(unlist(a)))])


all_equal2 <- function(v)
   all(unlist(sapply(1:length(v), function(x) v[x]==v[x:length(v)])))


make.string <- function(vct) {
	rtn <- vct[1]
	if (length(vct) > 1) for(i in 2:length(vct)) rtn <- paste(rtn, vct[i],sep=", ")
	
	rtn
}


right <- function(txt, i)
   substring(txt, nchar(txt)-i+1)


subplot <- function(x, y)
   viewport(layout.pos.row = x, layout.pos.col = y)





