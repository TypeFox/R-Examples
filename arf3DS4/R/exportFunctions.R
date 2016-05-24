#############################################
# arf3DS4 S4 DISPLAY FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#arf2Excel

arf2Excel <- 
function(mod,filename='estimates.txt')
{
	modmat <- cbind(1:.model.regions(mod),matrix(.model.estimates(mod),,.model.params(mod),byrow=T))
	dimnames(modmat)=list(1:.model.regions(mod),c('reg','x','y','z','sx','sy','sz','rxy','ryz','rxz','a'))
	write.table(modmat,file=filename,quote=F,sep='\t',col.names=T,row.names=F)
	
}
