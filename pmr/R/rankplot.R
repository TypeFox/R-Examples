rankplot <- function(dset, trans=FALSE, label.type="item", circle.col="black", circle.bg=FALSE){
##require(graphics)
nitem <- ncol(dset)-1
test <- matrix(data = 0, nrow = factorial(nitem), ncol = nitem, byrow = TRUE)
temp1 <- 1:nitem
i <- 1

if (nitem<3 || nitem>4){
message("Only ranking data with 3 or 4 items is supported.")
}

## generate a list of all possible rankings
for (j in 1:(nitem^nitem-1)){
temp1[1] <- nitem - j%%nitem
temp2 <- j - j%%nitem
for (k in nitem:2){
temp1[k] <- nitem - temp2%/%(nitem^(k-1))
temp2 <- temp2 - (nitem-temp1[k])*(nitem^(k-1))
}
temp2 <- 0
for (l in 1:nitem){
for (m in 1:nitem){
if (temp1[l] == temp1[m] && l != m){
temp2 <- 1
}
}
}
if (temp2 == 0){
for (p in 1:nitem){
test[i,p] = temp1[p]
}
i <- i+1
}
}


n <- rep(0,factorial(nitem))
for (j in 1:factorial(nitem)){
for (k in 1:nrow(dset)){
temp_ind <- 0
for (l in 1:nitem){
if (test[j,l] != dset[k,l]) {temp_ind <- temp_ind + 1}
}
if (temp_ind == 0) {n[j] <- dset[k,nitem+1]}
}
}
test2 <- cbind(test, n)

if (nitem==3){
label <- c(123,132,312,321,231,213)
freq <- rep(0,6)
for (i in (1:6)){
for (j in (1:6)){
if(test2[i,1]*100+test2[i,2]*10+test2[i,3] == label[j]){
freq[j] <- test2[i,4]
}
}
}

## A|B|C as label
if (label.type=="item"){
label[1] <- "ABC"
label[2] <- "ACB"
label[3] <- "BCA"
label[4] <- "CBA"
label[5] <- "CAB"
label[6] <- "BAC"
}

## rank as label
##if (label.type=="ordering"){
##label[1] <- "123"
##label[2] <- "213"
##label[3] <- "132"
##label[4] <- "231"
##label[5] <- "312"
##label[6] <- "321"
##}

x = c(0,1,1.5,1,0,-0.5)
y = c(0,0,-0.866,-1.732,-1.732,-0.866)
graphics::plot(x,y,axes=FALSE,ann=FALSE,xaxt="n",yaxt="n", xlim=c(-1.482,2.982), ylim=c(-2.232,2.232))
graphics::symbols(x,y,circles=freq^0.5/sum(freq^0.5)/2,inches=FALSE, add=TRUE, fg=circle.col, bg=circle.bg)
graphics::text(x,y,labels=label)
graphics::polygon(x,y)
}

if (nitem==4){
plot(c(0,3), c(0,-3), xlab=" ", ylab=" ", type="n", axes=FALSE,ann=FALSE,xaxt="n",yaxt="n")
# draw the labels and the circles
label <- rep(0,24)
freq <- rep(0,24)
for (i in (1:24)){
label[i] <- test2[i,1]*1000+test2[i,2]*100+test2[i,3]*10+test2[i,4]
}
for (i in (1:24)){
for (j in (1:24)){
if(test2[i,1]*1000+test2[i,2]*100+test2[i,3]*10+test2[i,4] == label[j]){
freq[j] <- test2[i,5]
}
}
}

# ordering as label
if (label.type=="ordering"){
label[1] <- "1234"
label[2] <- "2134"
label[3] <- "1324"
label[4] <- "2314"
label[5] <- "3124"
label[6] <- "3214"
label[7] <- "1243"
label[8] <- "2143"
label[9] <- "1342"
label[10] <- "2341"
label[11] <- "3142"
label[12] <- "3241"
label[13] <- "1423"
label[14] <- "2413"
label[15] <- "1432"
label[16] <- "2431"
label[17] <- "3412"
label[18] <- "3421"
label[19] <- "4123"
label[20] <- "4213"
label[21] <- "4132"
label[22] <- "4231"
label[23] <- "4312"
label[24] <- "4321"
}
# A|B|C|D as label
if (label.type=="item"){
label[1] <- "ABCD"
label[2] <- "BACD"
label[3] <- "ACBD"
label[4] <- "BCAD"
label[5] <- "CABD"
label[6] <- "CBAD"
label[7] <- "ABDC"
label[8] <- "BADC"
label[9] <- "ACDB"
label[10] <- "BCDA"
label[11] <- "CADB"
label[12] <- "CBDA"
label[13] <- "ADBC"
label[14] <- "BDAC"
label[15] <- "ADCB"
label[16] <- "BDCA"
label[17] <- "CDAB"
label[18] <- "CDBA"
label[19] <- "DABC"
label[20] <- "DBAC"
label[21] <- "DACB"
label[22] <- "DBCA"
label[23] <- "DCAB"
label[24] <- "DCBA"
}
if (trans==FALSE){
x <- c(0.76,1.26,0.44,1.51,0.64,1.24,1.15,1.61,0.55,2.2,0.77,1.98,1.27,2.23,1,2.55,1.53,2.15,1.88,2.38,1.66,2.72,1.98,2.56)
y <- c(-1.02,-0.56,-1.22,-0.24,-0.95,-0.4,-1.43,-1,-1.93,-0.4,-1.7,-0.58,-2.04,-1.15,-2.33,-0.88,-1.93,-1.36,-2.23,-1.77,-2.55,-1.55,-2.37,-1.83)
}
else {
x <- c(1,1.66,0.55,1.98,0.77,1.53,1.27,1.88,0.44,2.56,0.64,2.15,1.15,2.38,0.76,2.72,1.24,1.98,1.61,2.23,1.26,2.55,1.51,2.2)
y <- c(-2.33,-2.55,-1.93,-2.37,-1.70,-1.93,-2.04,-2.23,-1.22,-1.83,-0.95,-1.36,-1.43,-1.77,-1.02,-1.55,-0.4,-0.58,-1,-1.15,-0.56,-0.88,-0.24,-0.4)
}
graphics::symbols(x,y,circles=freq^0.5/sum(freq^0.5)/2,inches=FALSE,add=TRUE, fg=circle.col, bg=circle.bg)
graphics::text(x,y+0.1,labels=label)

# draw the truncated octahedron
Ax = c(0.64,0.77,1.53,2.15,1.98,1.24)
Ay = c(-0.95,-1.7,-1.93,-1.36,-0.58,-0.4)
graphics::polygon(Ax,Ay)
Dx = c(1.15,1.27,1.88,2.38,2.23,1.61)
Dy = c(-1.43,-2.04,-2.23,-1.77,-1.15,-1)
graphics::polygon(Dx,Dy,lty=3)
Ex = c(1.98,2.15,2.56,2.72,2.55,2.2)
Ey = c(-0.58,-1.36,-1.83,-1.55,-0.88,-0.4)
graphics::polygon(Ex,Ey)
Fx = c(0.55,1,1.66,1.98,1.53,0.77)
Fy = c(-1.93,-2.33,-2.55,-2.37,-1.93,-1.7)
graphics::polygon(Fx,Fy)
Gx = c(1.66,1.98,2.56,2.72)
Gy = c(-2.55,-2.37,-1.83,-1.55)
graphics::lines(Gx,Gy)
Hx = c(0.55,0.44,0.64,1.24,1.51,2.2)
Hy = c(-1.93,-1.22,-0.95,-0.4,-0.24,-0.4)
graphics::lines(Hx,Hy)
Ix = c(0.44,0.76,1.26,1.51)
Iy = c(-1.22,-1.02,-0.56,-0.24)
graphics::lines(Ix,Iy,lty=3)
Jx = c(1.26,1.61)
Jy = c(-0.56,-1)
graphics::lines(Jx,Jy,lty=3)
Kx = c(0.76,1.15)
Ky = c(-1.02,-1.43)
graphics::lines(Kx,Ky,lty=3)
Lx = c(1,1.27)
Ly = c(-2.33,-2.04)
graphics::lines(Lx,Ly,lty=3)
Mx = c(1.66,1.88)
My = c(-2.55,-2.23)
graphics::lines(Mx,My,lty=3)
Nx = c(2.23,2.55)
Ny = c(-1.15,-0.88)
graphics::lines(Nx,Ny,lty=3)
Ox = c(2.38,2.72)
Oy = c(-1.77,-1.55)
graphics::lines(Ox,Oy,lty=3)

}

}