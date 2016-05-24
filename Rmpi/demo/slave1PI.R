slave1 <- function (){
    n <- mpi.bcast(integer(1),1,comm=.comm)
    request <-1
    job <-2
    anysrc <- mpi.any.source()
    anytag <- mpi.any.tag()
    id <- mpi.comm.rank(.comm)
    tsize <- mpi.comm.size(.comm)
    if (id == tsize-1){ 	#the server that collects results
	mypi <- 0
	done <- 0
	while (done < tsize-2){
	   mypi <- mpi.recv(double(1), type=2, source=anysrc,
			tag=anytag,comm=.comm)+mypi
	   tag <- mpi.get.sourcetag()[2]
	   done <- done+(tag==0)
	} 
	#send the final result to master
	mpi.send(mypi,2,dest=0,tag=3,comm=.comm)
    }
    else {	# the workers that do the real jobs
    	while (1) {
    	    #send master a request
    	    mpi.send(integer(1),type=1,dest=0,tag=request,comm=.comm)
    	    jobrange<-mpi.recv(integer(2),type=1,source=0,
			tag=anytag,comm=.comm)
	    tag <- mpi.get.sourcetag()[2]
	    if (tag==job){
	    	mypi <- 4*sum(1/(1+((seq(jobrange[1],jobrange[2])-.5)/n)^2))/n
	    	mpi.send(mypi, 2, dest=tsize-1, tag=job,comm=.comm)	
	    }
	    else { #tag=0 means stop
	    	mpi.send(0, 2, dest=tsize-1, tag=0,comm=.comm)	
		break
	    }
        }
    }
}

master1PI <- function (n,maxjoblen,comm=1) 
{
    tsize <- mpi.comm.size(comm)
    if (tsize <3)
	stop("Need at least 2 slaves")
    #send the function slave1 to all slaves
     mpi.bcast.Robj2slave(slave1, comm=comm)

    #let slave run the function slave1
    mpi.bcast.cmd(slave1(), comm=comm)

    #send n to all slaves	
    mpi.bcast(as.integer(n),type=1,comm=comm)
    
    count <- 0
    request <- 1
    job <- 2
    anysrc <- mpi.any.source()
    while (count < n){
	mpi.recv(integer(1), type=1, source=anysrc, tag=request, comm=comm)
    	src <- mpi.get.sourcetag()[1]
	jobrange <- c(count+1, min(count+maxjoblen, n))
	mpi.send(as.integer(jobrange),type=1,dest=src,tag=job,comm=comm)
	count <- count+maxjoblen
    }
    #tell workers to stop with tag=0
    for (i in 1:(tsize-2)){
	mpi.recv(integer(1),type=1,source=anysrc,tag=request,comm=comm)
    	src <- mpi.get.sourcetag()[1]
	mpi.send(integer(1),type=1,dest=src,tag=0,comm=comm)
    }
    mpi.recv(double(1), type=2, source=tsize-1, tag=3, comm=comm)
}
