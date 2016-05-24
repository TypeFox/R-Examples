BclimLayer <-
function(pollen.loc,required.data3D,nsamples=1000){

  pollen.data <- read.table(pollen.loc,header=TRUE)
  
  nslices <- nrow(pollen.data)
  
  if(is.null(required.data3D)) stop("You need to supply a response surfaces file. See help(Bclim) for more details")
  
  ##### Format pollen data
	#getting rowSums ~ 1000
  myrowsums <- rowSums(pollen.data)
  myrowsums[myrowsums==0] <- 1
	Pollen<-as.matrix(round(pollen.data/myrowsums*1000))
	Mu<-V<-matrix(0,175616,28)

  ##### loading predictors & parameters
	M<-required.data3D$Mu;Var<-required.data3D$V;quadpts<-required.data3D$quadpts
	Id<-required.data3D$Id;Mu[Id,]<-M;V[Id,]<-Var;Mu<-c(Mu);V<-c(V)
	quadprobs<-required.data3D$quadprobs;nquadpts<-13;alpha<-required.data3D$alpha
	delta<-required.data3D$delta;Buffer<-required.data3D$Buffer
	Resultslices<-rep(0,length=(175616*nslices))
	Climate<-required.data3D$Climate

  ##### Running C function
	result<-.C("PalaeoRecon3D",as.integer(nslices),as.double(Resultslices),as.integer(Pollen),as.double(Mu),as.double(V), as.double(alpha),as.double(delta),as.double(Buffer), as.integer(nquadpts), as.double(quadpts), as.double(quadprobs)
  ,PACKAGE="Bclim")
		res<-result[[2]];res[res==0]<--Inf

  ##### return posterior dists in a matrix of dimension nslices x 2500
	post<-matrix(exp(res),nslices,175616,byrow=TRUE);post<-post/rowSums(post)

  #####provide n samples from MTCO/GDD5 space including jittering
  samples<-rep(0,nsamples*nslices*3)	

	for(i in 1:nslices){
		Firstsamp<-sample(1:175616, nsamples, replace = TRUE, prob = post[i,])
		jitter<-cbind(runif(nsamples,-1,1),runif(nsamples,-1,1),runif(nsamples,-1,1))*cbind(rep(72.53,nsamples),rep(.7,nsamples),rep(10.206,nsamples))

		jit_locations<-Climate[Firstsamp,]+jitter
		samples[(i-1)*nsamples + (1:nsamples)]<-jit_locations[,1]		
		samples[nslices*nsamples+(i-1)*nsamples + (1:nsamples)]<-jit_locations[,2]
		samples[nslices*nsamples*2+(i-1)*nsamples + (1:nsamples)]<-jit_locations[,3]
	}

	#restructuring samples array
	dim(samples)=c(nsamples,nslices,3)
	
	#naming dimensions
	dimnames(samples)=list(NULL,NULL,c("GDD5","MTCO","AET/PET"))

  return(samples)
}
