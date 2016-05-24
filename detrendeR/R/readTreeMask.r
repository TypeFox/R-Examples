readTreeMask =function (rwl, stc = c(5, 2, 1)) 
{
    if (sum(stc) != 8) 
        stop("Site-Tree-Core mask does not sum to 8")
    ids = colnames(rwl)							     #series names
  #  ids = format(ids,width=8, align="l")

test = function (x, site.chars=stc) {
out=c(NA,NA,NA)
	  out[1] = substring(x, 1, stc[1])
        out[2] = substring(x, stc[1]+1, sum(stc[1:2]))
        out[3] = substring(x, sum(stc[1:2])+1, sum(stc[1:2])+stc[3])
return(out)
}

out = t(sapply(ids, test,site.chars=stc))
out = data.frame(out)

tree.series = ids
tree.vec = as.numeric(out[, 2])
tree.ids = unique(out[, 2])

    core.vec = rep(NA, length(tree.vec))
    n.trees = length(tree.ids)
    for (i in 1:n.trees) {
        n.cores = length(core.vec[tree.vec == i])
        core.vec[tree.vec == i] = seq(1, n.cores)
    }
 
	out =data.frame(out,tree = tree.vec, core = core.vec)
	out<-out[order(out[,4]),]					         #sort by tree 
colnames(out) = c("Site", "Tree", "Core","tree", "core")
return(out)
}

#read.tree.mask(rwl)
