kFquant <-
function(inf=NULL, is.octopus=FALSE,
                       range.sex=NULL, range.age=NULL, range.qual=NULL,   
                       plot.iso='III4e', show.raw=FALSE, 
                       tau=c(0.025, 0.25, 0.5, 0.75, 0.975))
                      
{
### range.sex must be NULL (include all subjects), 'Female', or 'Male'
### range.age must be NULL (all ages), or a (numeric) a single value or a vector of length 2 specifying a closed range
### range.qual must be NULL (all subjects), 
###                  or a single value from {"Good witness", "Fair witness", "Poor witness"}
###                  or a vector of length 2 from that set
### plot.iso must be one of 'III4e', 'I4e', 'I2e'
###

if(is.null(inf)) goldmann.demogr<- Goldmann.demogr
else  goldmann.demogr <- inf


names(goldmann.demogr)[1]<- 'ID'


b0<- (1:dim(goldmann.demogr)[1])[!is.na(goldmann.demogr$Eye)] ### extracts positions individuals with data available
b00<-  !is.na(goldmann.demogr$Eye) ### logical - individuals with meaningful data



if(!is.null(range.sex)) b0a<- b00 & is.element(goldmann.demogr$Sex, range.sex)   else b0a<- rep(TRUE, length(b00))
if(!is.null(range.age)) b0b<- b00 & is.element(goldmann.demogr$Age, range.age)  else b0b<- rep(TRUE, length(b00))
if(!is.null(range.qual)) b0c<- b00 & is.element(goldmann.demogr$Quality.of.test, range.qual) else b0c<- rep(TRUE, length(b00))

b00<- b0a & b0b & b0c



goldmann.demogr<- goldmann.demogr[b00,]  ### select filtered individuals



list.of.files<- paste(goldmann.demogr$ID, substring(goldmann.demogr$Eye,1,1),sep='') ### generates list of files
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

                     aux<- kf.sector(i, is.octopus=is.octopus)

                     mat.output[s0, 3:5]<-aux$dists
                     mat.output[s0, 6:8]<-aux$freqs
              }





mat.output<- as.data.frame(mat.output)

mat.output$Sex<- rep( goldmann.demogr$Sex, rep( 24,num.files))
mat.output$Age<- rep( goldmann.demogr$Age, rep(24, num.files))
mat.output$Quality.of.test <- rep( goldmann.demogr$Quality.of.test, rep( 24,num.files))
#print(mat.output)
##### 

for (i in 3:5) mat.output[,i]<- ifelse( mat.output[,i]==0, NA, mat.output[,i])


ind.isopter<- switch(plot.iso,  'III4e' = 3, 'I4e' = 4,'I2e' = 5)
col.isopter<- switch(plot.iso,  'III4e' = 'black', 'I4e' = 'royalblue3',
                                    'I2e' = 'maroon')

iso.band.col<- switch(plot.iso,  'III4e' = 'gray45', 'I4e' = 'gray45',
                                    'I2e' = 'gray45')


arcsM<- rad(seq( 0,  345, by=15)) ### angles for centra lines of arcs

## To super-impose raw data points:                                                             
if(show.raw) kFcheck(mat.output, name.iso=plot.iso,
                           plot.lines=FALSE)
else(set.template())

if(plot.iso =='III4e')
{

test<- mat.output[, c(1,2,3)]
test<-data.frame( ID=test$ID, theta=arcsM[test$Sector],
                  r=test$dIII4e)
}



if(plot.iso =='I4e')
{

test<- mat.output[, c('ID','Sector', 'dI4e')]
test<-data.frame( ID=test$ID, theta=arcsM[test$Sector],
                  r=test$dI4e)
}


if(plot.iso =='I2e')
{

test<- mat.output[, c('ID','Sector', 'dI2e')]
test<-data.frame( ID=test$ID, theta=arcsM[test$Sector],
                  r=test$dI2e)
}



arcsX<-seq( 0, 2*pi, length=1000) ## fine seq to predict
arcsX<- data.frame(theta=arcsX)
arcsX<- unlist(arcsX)
c1<- cos(arcsX)
s1<- sin(arcsX)




mod3<- lqmm(  r ~ cos( theta) + sin(theta) + cos(2*theta) + sin(2*theta)
               + I(sin(theta)*cos(2*theta)), tau=tau,
               random= ~1, group=ID,  na.action=na.omit,
               data=test)
                           
pred3<- matrix(NA, ncol=length(tau), nrow=length( arcsX))

for (i in 1:length(tau))
{
beta<- mod3$theta_x
pred3[,i] <- beta[1,i] + beta[2,i]*cos( arcsX) + beta[3,i]*sin(arcsX) +
             beta[4,i]*cos(2*arcsX) + beta[5,i]*sin(2*arcsX) +
             beta[6,i]*sin(arcsX) * cos(2*arcsX)
}


innerx<-  pred3[,1]*c1
innery<-  pred3[,1]*s1

outerx<- pred3[,5]*c1
outery<- pred3[,5]*s1

inner<- cbind(innerx, innery)
outer<- cbind(outerx, outery)

inner<- rbind( inner, inner)
outer<- rbind( outer, outer)

x1<- c(outer[,1], rev(inner[,1]))
y1<- c(outer[,2], rev(inner[,2]))

test<- polygon(x1, y1, col='snow3', density=36, border='grey')



for ( i in seq(10,70,by=10)) 
{
draw.circle( 0,0, radius=i,border='gray', nv=500)
text( 0, -i,  i, cex=2/3, col='slategray')
text( 0, i, i, cex=2/3, col='slategray')
text( -i,2, i, cex=2/3, col='slategray')
text(i,2,i,cex=2/3, col='slategray')
}



text(c(-80,-90,80,90),rep(2,4), rep(c(80,90),2),
cex=2/3, col='slategray')

text(par()$usr[1] -4, 0, "180", col='darkslategray',cex=3/5, family='sans')
text(par()$usr[2] + 2, 0, "0", col='darkslategray',cex=3/5, family='sans')



lines( pred3[,1]*c1, pred3[,1]*s1,lwd=2,col=iso.band.col,lty=1)
lines( pred3[,3]*c1, pred3[,3]*s1,lwd=3,col= col.isopter,lty=1)
lines( pred3[,5]*c1, pred3[,5]*s1,lwd=2,col=iso.band.col,lty=1)
lines( pred3[,2]*c1, pred3[,2]*s1,lwd=1,col='blueviolet',lty=2)
lines( pred3[,4]*c1, pred3[,4]*s1,lwd=1,col='blueviolet',lty=2)

}
