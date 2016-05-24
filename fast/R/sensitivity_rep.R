`sensitivity_rep` <-
function(data.zoo, xval=index(data.zoo), direction=1, data=coredata(data.zoo), numberf, order=4, legend = paste("P",1:order, sep=""), cukier=TRUE, reorder=1:numberf, ...){
sen <- apply(data,direction,sensitivity, numberf=numberf, order=order, cukier=cukier, make.plot=FALSE, reorder=reorder)
p.sensitivity(sen=sen, xval=xval, legend=legend, ...)

return(sen)

}

