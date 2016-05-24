library(intamap)
library(automap)
library(gstat)
set.seed(13531)

plotFigs = FALSE
npts = 1000
pts=SpatialPoints(cbind(runif(npts),runif(npts)))
d = SpatialPointsDataFrame(cbind(0,0),data.frame(z=1))
observations = krige(z~1,d,pts,vgm(1, "Sph", .5,anis=c(90,0.5)),nsim=1,nmax=50,beta=0)

#spplot(sim)

predictionLocations = spsample(observations, 1000, "regular")

# We dont know the projection of the data at this stage, assume it is
# somehow metric



formulaString = as.formula(sim1~1)
krigingObject = createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
  formulaString = formulaString
)
class(krigingObject) = c("automap")

checkSetup(krigingObject)
object = preProcess(krigingObject)
objTemp = estimateAnisotropy(object)
#rotate Data
objTemp$observations = rotateAnisotropicData(objTemp$observations, objTemp$anisPar)
#Estimate Variogram Model
vario = autofitVariogram(objTemp$formulaString, objTemp$observations, model = "Sph")$var_model
objTemp$anisPar
if (plotFigs) {
  spplot(object$observations, "sim1", col.regions=bpy.colors())
  spplot(objTemp$observations, "sim1", col.regions=bpy.colors())
  plot(variogram(sim1~1, object$observations, alpha=c(0, 90)), vario)
}
vario


vmod = vgm(1, "Sph", 1, anis = c(90, 0.5))
krigingObject$observations = krige(z~1, d, pts, vmod, nsim = 1, nmax = 50, beta = 0)
object = preProcess(krigingObject)
objTemp = estimateAnisotropy(object)
objTemp$anisPar
objTemp$observations = rotateAnisotropicData(objTemp$observations, objTemp$anisPar)
vario = autofitVariogram(objTemp$formulaString, objTemp$observations, model="Sph")$var_model
vario
vmod

vmod = vgm(1, "Sph", 2, anis=c(45, 0.2))
krigingObject$observations = krige(z~1, d, pts, vmod, nsim = 1, nmax = 50, beta = 0)
object = preProcess(krigingObject)
objTemp = estimateAnisotropy(object)
objTemp$anisPar
objTemp$observations = rotateAnisotropicData(objTemp$observations, objTemp$anisPar)
vario = autofitVariogram(objTemp$formulaString, objTemp$observations, model = "Sph")$var_model
vario
vmod


vmod = vgm(1,"Sph",1,anis=c(135,0.5))
krigingObject$observations = krige(z~1,d,pts,vmod,nsim=1,nmax=50,beta=0)
object = preProcess(krigingObject)
objTemp=estimateAnisotropy(object)
objTemp$anisPar
objTemp$observations=intamap:::rotateAnisotropicData(objTemp$observations,objTemp$anisPar)
vario = autofitVariogram(objTemp$formulaString,objTemp$observations,model="Sph")$var_model
vario
vmod



vmod = vgm(1,"Sph",.3,anis=c(90,0.5))
krigingObject$observations = krige(z~1,d,pts,vmod,nsim=1,nmax=100,beta=0)
object = preProcess(krigingObject)
objTemp=estimateAnisotropy(object)
objTemp$anisPar
objTemp$observations=intamap:::rotateAnisotropicData(objTemp$observations,objTemp$anisPar)
vario = autofitVariogram(objTemp$formulaString,objTemp$observations,model="Sph")$var_model
vario
vmod
if (plotFigs) {
  p1 = plot(variogram(sim1~1,object$observations,alpha=c(0,90)),vmod,ylim = c(0,1.2),xlim=c(0,0.6),main="orig,orig")
  p2 = plot(variogram(sim1~1,object$observations,alpha=c(0,90)),vario,ylim = c(0,1.2),xlim=c(0,0.6),main="orig,fitted")
  p3 = plot(variogram(sim1~1,objTemp$observations,alpha=c(0,90)),vmod,ylim = c(0,1.2),xlim=c(0,0.6),main="rot,orig")
  p4 = plot(variogram(sim1~1,objTemp$observations,alpha=c(0,90)),vario,ylim = c(0,1.2),xlim=c(0,0.6),main = "rot,fitted")

  print(p1,position = c(0,0.5,0.5,1),more = TRUE)
  print(p2,position = c(0.5,0.5,1,1),more = TRUE)
  print(p3,position = c(0,0,0.5,0.5),more = TRUE)
  print(p4,position = c(0.5,0,1,0.5))

  plot(variogram(sim1~1,objTemp$observations,alpha=c(seq(0,150,30))),vmod,ylim = c(0,1.2),xlim=c(0,0.6),main="rot,orig")

  spplot(object$observations,"sim1",col.regions=bpy.colors())
  spplot(objTemp$observations,"sim1",col.regions=bpy.colors())
}




data(sic2004)
coordinates(sic.val)=~x+y
sic.val$value=sic.val$dayx
x = sic.test$x
y = sic.test$y

coordinates(sic.test)=~x+y

stest = sic.test[(x > -10000 & x < 140000 & y > 100000 & y < 240000),]
 
obj<-createIntamapObject(formulaString = "joker~1",
   observations=sic.val,
   predictionLocations = stest,
   params = list(doAnisotropy = TRUE),
   class = "automap" )
obj = preProcess(obj)
obj = estimateParameters(obj)
obj$anisPar
obj = spatialPredict(obj)
obj = postProcess(obj)
summary(as.data.frame(obj$outputTable))
