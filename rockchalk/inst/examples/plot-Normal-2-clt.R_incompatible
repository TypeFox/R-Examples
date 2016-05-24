### Filename: Normal2_2008.R
### Paul Johnson  March 28, 2008

### This code should be available somewhere in
### http://pj.freefaculty.org/R.  If it is not email me
### <pauljohn@ku.edu>

sampleSize <- 500

mymean <- 0

mystd <- 1.5

myx <- seq( mymean - 3*mystd, mymean + 3*mystd, length.out = sampleSize)

## Creates an empty container to hold samples:
mySamples <- list()

## for loop creates 100 samples and stores them
for (i in 1:100) {
  mySamples[[i]] <- rnorm(sampleSize,  mean = mymean, sd = mystd)
}


par(mfcol = c(5,2))

for (i in 1:10) {
  hist(mySamples[[i]],main = paste("Sample", i) )
}


### force ranges to match
for (i in 1:10) {
  mytitle <- paste("Sample", i, "mean", round(mean(mySamples[[i]]),2),"std.dev.", round(sd(mySamples[[i]]),2) )
  myhist <- hist(mySamples[[i]],main= mytitle, freq = FALSE, xlim = c(-6,6), xlab = "Normally Distributed Variable" )
}





par(mfcol = c(1,1))



### sapply applies the mean function to each sample. The "s" in
### sapply stands for "simplify".
myMeans <- sapply(mySamples, mean)

hist(myMeans,freq = FALSE, main = paste("Sampling Distribution of ",100," Samples"))

mtext(text = paste("True Mean = ", mymean,  "True Standard Deviation = ", mystd) , side = 3)

### Note the range of the means:
summary(myMeans)

### Compare against the range in one sample, say the 5th
summary(mySamples[[5]])

text( -0.2, 2.5, paste("Mean of Means =", round(mean(myMeans),3),"\n Std. of Means ", round(sd(myMeans),3)))


par(mfcol = c(2,1))

 myhist <- hist(mySamples[[1]],main= mytitle, freq = FALSE, xlim = c(-6,6), xlab = "Normally Distributed Variable" )


hist(myMeans,freq = FALSE, main = paste("Sampling Distribution of ",100," Samples"), xlim = c(-6,6))


### Another way. write a function to create the data and draw it.
### Use lapply to call it over and over




sampleSize <- 500

par(mfrow = c(5,2))


createDist <-function(i){
  z <- rnorm(sampleSize,  mean = mymean, sd = mystd)
  title <- paste("Run", i, "Normal mean=",mymean," std.dev=",mystd)
  hist(z,breaks = 20,main = title)
}

createDist(1)
createDist(2)
createDist(3)
createDist(4)
createDist(5)
createDist(6)
createDist(7)
createDist(8)
createDist(9)
createDist(10)


lapply(1:10,createDist);

par(mfrow = c(1,1))









### Another way to show a sampling distribution. This uses R's built in "Replicate"


sampleSize <- 10000

repls <- 500

mytitle <- paste("Replications=",repls,"SampleSize=",sampleSize,"Mean=",mymean,"Std.Dev=",mystd)

means <- replicate(repls,mean(rnorm(sampleSize,mean = mymean,sd = mystd)))

myhist <- hist(means, main = mytitle, breaks = 15)

# Display some results in text inside the histo

# I know there should be a better way, but don't know it

# Note that in text, \n has a "next line" effect

mylabel <- paste("E(x)=",mymean,"\nobserved

mean=",round(mean(means),2),"\nobserved std.dev \nof means=",round(sd(means),2))

# This puts text into the histogram, guessing

# the positions from the 6'th bin's position

text(myhist$breaks[6],.8*max(myhist$counts),mylabel ,pos = 2)







sampleSize <- 10000

repls <- 500

mytitle <- paste("uniforms")

means <- replicate(repls,mean(runif(sampleSize,min = 0,max = 222)))

par(mfcol = c(5,2))

for (i in 1:10)
hist(runif(sampleSize, min = 0, max = 222))

x11()

myhist <- hist(means, main = mytitle, breaks = 15)

# Display some results in text inside the histo

# I know there should be a better way, but don't know it

# Note that in text, \n has a "next line" effect

mylabel <- paste("E(x)=",mymean,"\nobserved

mean=",round(mean(means),2),"\nobserved std.dev \nof means=",round(sd(means),2))

# This puts text into the histogram, guessing

# the positions from the 6'th bin's position

text(myhist$breaks[6],.8*max(myhist$counts),mylabel ,pos = 2n)
