formatRESULT = function(dat, triplicate="Triplicate", score="score", t="Time"){
dat1 = dat[dat[[triplicate]]==1,]

urep = sort(unique(dat[[triplicate]]))
dat.scores = sapply(urep, function(x){dat[dat[[triplicate]]==x,][[score]]})
colnames(dat.scores) = paste(score, urep, sep="")
dat = data.frame(dat1[,colnames(dat1)!=triplicate & colnames(dat1)!=score], dat.scores)

datbefore = dat[as.character(dat[[t]])=="Before",]
colnames(datbefore)[grep(score, colnames(datbefore))] = paste(score, "before", urep, sep="")

datafter = dat[as.character(dat[[t]])=="After",]
colnames(datafter)[grep(score, colnames(datafter))] = paste(score, "after", urep, sep="")

res = data.frame(datbefore[,c("MainPlate","Plate", "Norm", "well", "row", "col", "welltype")], datbefore[,paste(score, "before", urep, sep="")], datafter[,paste(score, "after", urep, sep="")])

res = data.frame(ID=paste(res[["MainPlate"]], res[["Plate"]], res[["well"]], sep="_"), res)
return(res)
}
