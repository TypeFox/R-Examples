 simslave <- function (){
    request <-1
    result <-2
    chunksize <- 1000 #job size for slaves
    anytag <- mpi.any.tag()
    while (1) { #slaves do their jobs first and waiting instructions from master
	x <- runif(chunksize)
	y <- runif(chunksize)
	incirvec <- (x-.5)^2 + (y-.5)^2 < .25
	xy <- c(x[incirvec],y[incirvec])
	#send the result to master
	mpi.send(xy, 2, dest=0, tag=result,comm=.comm)	
	#receive instructions from master
	mpi.recv(integer(1),type=1,source=0,tag=anytag,comm=.comm)
	tag <- mpi.get.sourcetag()[2]
	if (tag!=request)
	    break
    }
}

simPI <- function (n,epsilon=1e-4,comm=1) 
{
    tsize <- mpi.comm.size(comm)
    if (tsize < 2)
	stop("It seems no slaves running")
   
    #send the function simslave to all slaves
    mpi.bcast.Robj2slave(simslave, comm=comm)

    #let slaves run the function simslave
    mpi.bcast.cmd(simslave(), comm=comm)

    chunksize <- 1000
    request <-1
    result <- 2
    count <- 0
    totalincir <- 0
    anysrc <- mpi.any.source()

   #prepare an empty plot
    plot(c(0,0),xlim=c(0,1),ylim=c(0,1), ylab="",type="n")

    while (1) {
	#receive done job from slaves
   	xy<-mpi.recv(double(2*chunksize),type=2,source=anysrc,
			tag=result,comm=comm)
	#receive buffer is biger than actual data. So need to get real length
	incir <-mpi.get.count(2)/2
	src <- mpi.get.sourcetag()[1]
	totalincir <- incir + totalincir
	count <- count+chunksize
	mypi <- 4*totalincir/count
	x <- xy[1:incir]
	y<- xy[(incir+1):(2*incir)]
	#add incircle points to the plot
	points(x,y,pch=".",cex=.2,col=src)
	err <- abs(mypi-pi)    
	if (err > epsilon && count < n)
	    mpi.send(integer(1),type=1,dest=src,tag=request,comm=comm)
	else { 
	    #tell slaves to stop with tag=0
	    mpi.send(integer(1),type=1,dest=src,tag=0,comm=comm)
	    break
        }
    }
    #only one slave is stopped. So have to others to stop as well
    if (tsize > 2){
    for (i in 1:(tsize-2)){
	#continue receiving other slaves jobs
   	xy<-mpi.recv(double(2*chunksize),type=2,source=anysrc,
			tag=result,comm=comm)
	incir <-mpi.get.count(2)/2
	src <- mpi.get.sourcetag()[1]
	totalincir <- incir + totalincir
	count <- count+chunksize
	mypi <- 4*totalincir/count
	x <- xy[1:incir]
	y<- xy[(incir+1):(2*incir)]
	points(x,y,pch=".",cex=.2,col=src)
	#tell other slaves to stop
	mpi.send(integer(1),type=1,dest=src,tag=0,comm=comm)
    }
    }
    mypi
}
