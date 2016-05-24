require(adimpro)
if(Sys.getenv("ImageMagick")==""){
cat("Please install ImageMagick\n")
} else {
#
#     smoothing images
#
X11(width=8,height=4)
par(mfrow=c(1,2),mar=c(2,2,2,.2),mgp=c(2,1,0))
cat("Read img1\n")
img1 <- read.image(system.file("img/wias.ppm",package="adimpro"))
show.image(img1)
cat("Generate img2 containing oure noise\n")
img2 <- make.image(array(rnorm(30000,0,c(rep(3000,10000),rep(5000,10000),
rep(8000,10000))),c(100,100,3)))
show.image(img2)
readline("Press enter to add the noisy image to img1")

addrandom <- function(x,y){
means <- apply(y,3,mean)
x+sweep(y,3,means)
}
cat("Using fun=''+'' would not work correctly since images can not contain negative grey/color values \n")
img3 <- combine(img1,img2,fun=addrandom,rescale = FALSE)
show.image(img3)
cat("Generate img4 containing colored noise\n")
img4 <- make.image(array(rnorm(30000,0,c(rep(6000,10000),rep(10000,10000),
rep(16000,10000))),c(100,100,3)))
img4 <- awsimage(img4,hmax=1.5,aws=FALSE)
cat("add img4 to img1\n")
addrandom2 <- function(x,y,scale){
means <- apply(y,3,mean)
x+scale*sweep(y,3,means)
}
img5 <- combine(img1,img4,fun=addrandom2,rescale = FALSE,scale=2)
show.image(img5)

readline("Press enter for an example using masks")
simg1 <- segment(img3,hmax=20,select=TRUE,connect=TRUE)
simg1data <- extract.image(simg1)
levels <- as.integer(names(table(simg1data)))
mask <- simg1data==levels[length(levels)]
mask1 <- make.image(mask)
mask2 <- make.image(!mask)
prod <- function(x,y) x*(y/65535)
# image values are in 0:65535
img6 <- combine(img3,mask1,prod)
img7 <- combine(img3,mask2,prod)
add <- function(x,y,a,b) a*x+b*y
img8 <- combine(img6,img7,a=1.2,b=.6)
show.image(img8)
}

