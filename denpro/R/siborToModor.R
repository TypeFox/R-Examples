siborToModor<-function(tree){
#
# From ordering in siblings to ordering of modes
# We have the right ordering in profile
#
data<-plotprof(tree,plot=F,data=T,cutlev=NULL,ptext=0,info=NULL,
infolift=0,infopos=0)
vecs<-data$vecs
#
parent<-tree$parent
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc
moodinum<-mlkm$lkm    #length(modloc)
#
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i]
    xcor[i]<-vecs[loc,1]    
}
modloc<-omaord2(modloc,xcor)       
#
return(modloc)
}
