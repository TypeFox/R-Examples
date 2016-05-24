if(getRversion() >= "2.15.1")  utils::globalVariables("Goldmann.demogr")

kf.sort <-
function(inf=NULL, is.octopus=FALSE,
                   range.sex=NULL, range.age=NULL, range.qual=NULL,   
                   plot.isopter='III4e', CI.or.Quant='CI', 
                   force23=TRUE)
{
### This function completes the data sorting process
### It compiles all individual files into a matrix  specified by
### sex, age and quality variables. This allows the user to analyse data

### range.sex must be NULL (include all subjects), 'Female', or 'Male'
### range.age must be NULL (all ages), or a (numeric) a single value or a vector of length 2 specifying a closed range
### range.qual must be NULL (all subjects), 
###                  or a single value from {"Good witness", "Fair witness", "Poor witness"}
###                  or a vector of length 2 from that set
### plot.isopter must be one of 'III4e', 'I4e', 'I2e'
### CI.or.Quant is either 'CI' or 'Quant' for 95% CI's or 95% quantile envelope
### force23 is a binary indicator - FALSE to define the closure
###         of the bands at sectors 23 and 1; TRUE to define it at sectors 23 and 2
###
### The graphic output of this function is just to visualise data
### To get estimates and proper CI's use function quantile.fit


if(is.null(inf)) {goldmann.demogr<- Goldmann.demogr} 
 else            {goldmann.demogr <- inf}


names(goldmann.demogr)[1]<- 'ID'


b0<- (1:dim(goldmann.demogr)[1])[!is.na(goldmann.demogr$Eye)] ### extracts positions individuals with data available
b00<-  !is.na(goldmann.demogr$Eye) ### logical - individuals with meaningful data



if(!is.null(range.sex)) b0a<- b00 & is.element(goldmann.demogr$Sex, range.sex)   else b0a<- rep(TRUE, length(b00))
if(!is.null(range.age)) b0b<- b00 & is.element(goldmann.demogr$Age, range.age)  else b0b<- rep(TRUE, length(b00))
if(!is.null(range.qual)) b0c<- b00 & is.element(goldmann.demogr$Quality.of.test, range.qual) else b0c<- rep(TRUE, length(b00))

b00<- b0a & b0b & b0c



goldmann.demogr<- goldmann.demogr[b00,]  ### select filtered individuals



list.of.files<- paste(goldmann.demogr$ID, substring(goldmann.demogr$Eye,1,1),sep='') 
### generates list of files or objects, if inf is NULL

num.files<- length(list.of.files)                            

### declare space for distances and frequencies over sectors
mat.output<- matrix(0, nrow=num.files*24, ncol=8) ### ncolo = id, sector, dists x 3, freqs x 3

name.aux<- rep( list.of.files, rep( 24, num.files))
sec.aux<- rep(paste('S',1:24,sep=''), num.files)
names.aux<- paste( name.aux,sec.aux, sep='-')
names1<- c('III4e','I4e','I2e')
dimnames(mat.output)<- list(names.aux, 
c('ID', 'Sector',paste('d',names1,sep=''), paste('f',names1,sep='')))


mat.output[,2]<- rep(1:24, num.files)


k0<-0

for (i in list.of.files)
{
k0<-k0+1
s0<- seq(1 + (k0-1)*24, k0*24)
mat.output[s0, 1]<-rep(k0, 24)

aux<- kf.sector(file.name=i, is.octopus=is.octopus)

mat.output[s0, 3:5]<-aux$dists
mat.output[s0, 6:8]<-aux$freqs
}### close loop on list.of.files

### close is.goldmann

mat.output<- as.data.frame(mat.output)

mat.output$Sex<- rep( goldmann.demogr$Sex, rep( 24,num.files))
mat.output$Age<- rep( goldmann.demogr$Age, rep(24, num.files))
mat.output$Quality.of.test <- rep( goldmann.demogr$Quality.of.test, rep( 24,num.files))

##### 

arcsM<- rad(seq( 0,  345, by=15)) ### angles for central lines of arcs


ind.isopter<- switch(plot.isopter,  'III4e' = 3, 'I4e' = 4,'I2e' = 5)
col.isopter<- switch(plot.isopter,  'III4e' = 'black', 'I4e' = 'royalblue3',
                                    'I2e' = 'maroon')


b1<- mat.output[, ind.isopter]>0 

### get confidence intervals

if(CI.or.Quant == 'CI')
{
   test1<-tapply( mat.output[b1,ind.isopter], mat.output[b1, 'Sector'], mean)
   test2<-sqrt( tapply( mat.output[b1,ind.isopter], mat.output[b1,'Sector'], var)) #sd

 c1<- min( mat.output[b1, 'Sector'])
 
  if(c1>1)
   {
      test1<- c( NA, test1)
      test2<- c( NA, test2)
   }
   inner<- cbind((test1 -1.96*test2)* cos(arcsM), (test1 -1.96*test2)* sin(arcsM))
   outer<-  cbind((test1 +1.96*test2)* cos(arcsM), (test1 +1.96*test2)* sin(arcsM))
   middle<- cbind(test1*cos(arcsM), test1*sin(arcsM))
      
          
   ind.closure<-1
   if (is.na(inner[1,1]) | force23) ind.closure <- 2
  
    inner<- rbind( inner, inner[ind.closure,])
    middle<- rbind( middle, middle[ind.closure,])
    outer<- rbind( outer, outer[ind.closure,])
   
  

   set.template()
   if(force23)
   {
      lines(inner[-1,], col='firebrick4', lwd=2)
      lines(middle[-1,], col=col.isopter, lwd=2)
      lines(outer[-1,], col='firebrick4', lwd=2)
   }
  if(!force23)
   {
      lines(inner, col=2, lwd=2)
      lines(outer, col=2, lwd=2)
      lines(middle, col=2, lwd=2)
       }

   x1<- c(outer[,1], rev(inner[,1]))
   y1<- c(outer[,2], rev(inner[,2]))

 ### end of CI procedure
}

if(CI.or.Quant == 'Quant')
{
#### get quantiles

   q1<- tapply( mat.output[b1,ind.isopter], mat.output[b1, 'Sector'], quantile, 
            probs=c(0.025, 0.5, 0.975))
                       
   q1<- matrix( unlist(q1), ncol=3, byrow=TRUE) ### ok

   
 c1<- min( mat.output[b1, 'Sector'])
   
   if(c1>1)  q1<- rbind( rep(NA,3), q1) ### add 1st sector

   inner<- cbind( q1[,1]*cos(arcsM), q1[,1]*sin(arcsM))
   middle<- cbind( q1[,2]*cos(arcsM), q1[,2]*sin(arcsM))
   outer<- cbind( q1[,3]*cos(arcsM), q1[,3]*sin(arcsM))

   ind.closure<-1
   if (is.na(inner[1,1]) | force23) ind.closure <- 2
   
   
    inner<- rbind( inner, inner[ind.closure,])
    middle<- rbind( middle, middle[ind.closure,])
    outer<- rbind( outer, outer[ind.closure,])
  
   set.template()
   
   if(force23)
   {
      lines(inner[-1,], col='firebrick4', lwd=2)
      lines(middle[-1,], col=col.isopter, lwd=2)
      lines(outer[-1,], col='firebrick4', lwd=2)
   }
   if(!force23)
   {
      lines(inner, col='firebrick4', lwd=2)
      lines(middle, col=col.isopter, lwd=2)
      lines(outer, col='firebrick4', lwd=2)
   }

   x1<- c(outer[,1], rev(inner[,1]))
   y1<- c(outer[,2], rev(inner[,2]))
} ### end of quantiles procedure

test<- polygon(x1, y1, col='snow3', density=36, border='grey')


if(force23)
   {
      lines(inner[-1,], col='firebrick4', lwd=2)
      lines(middle[-1,], col=col.isopter, lwd=2)
      lines(outer[-1,], col='firebrick4', lwd=2)
   }
if(!force23)
   {
      lines(inner, col='firebrick4', lwd=2)
      lines(middle, col=col.isopter, lwd=2)
      lines(outer, col='firebrick4', lwd=2)
   }

       
areas<- c(areapl( na.omit(inner)), 
          areapl( na.omit(middle)), 
                 areapl( na.omit(outer))) ## from library splancs
names( areas)<- c('inner','middle','outer')            

####  change distances == 0 to NA's

for (i in 3:5) mat.output[,i]<- ifelse( mat.output[,i]==0, NA, mat.output[,i])

### compute rose diagram of frequencies 


mat.output<- as.data.frame(mat.output)


invisible(list(mat.output=mat.output, regions=list(inner=inner,middle=middle,outer=outer), areas=areas))

}
