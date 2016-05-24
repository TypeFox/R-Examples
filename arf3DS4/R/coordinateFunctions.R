#############################################
# arf3DS4 S4 COORDINATE FUNCTIONS			#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#cropVolume					
#cropVolumeAuto
#makeROImask						[user]
#readFSLmat
#createRegs							[user]
#setRegFiles						[user]
#setRegParams						[user]
#arf2MNI
#MNI2arf
#arf2high
#MNI2atlas
#reqFlip
#flipAxis
#euclidDist							[user]
#makeLowResStruct					[user]
#makeLowResStructAvg				[user]
#getTalairach
#getHarvardOxford
#getAtlasLabels						
#calcHarvardOxfordProbs
#calcTalairach
#getModelAtlas						[user]
#highresBlobs
#loadReg							[user]
#saveReg							[user]



cropVolume <- 
function(filename,resizeToDim,quiet=F) 
#resize nifti images (crops last slices to destination dimensions)
{
	
	dat <- readData(filename)
	data_array <- .fmri.data.datavec(dat)
	
	x=.fmri.data.dims(dat)[2]
	y=.fmri.data.dims(dat)[3]
	z=.fmri.data.dims(dat)[4]
	t=.fmri.data.dims(dat)[5]
	
	dim(data_array) <- c(x,y,z,t)
	
	nx=resizeToDim[2]
	ny=resizeToDim[3]
	nz=resizeToDim[4]
	
	doCrop=TRUE
	
	if(nx>x | ny>y | nz>z) {warning('can only crop images!');doCrop=FALSE}
	if(nx==x & ny==y & nz==z) {doCrop=FALSE}
	
	if(doCrop) {
		
		if(!quiet) cat('Dims',.fmri.data.dims(dat),'>> ')
			
		if(nx<x) {
			remx <- seq(nx+1,x,1)
			data_array <- data_array[-remx,,,]
			.fmri.data.dims(dat)[2]=nx
		} 
		
		if(ny<y) {
			remy <- seq(ny+1,y,1)
			data_array <- data_array[,-remy,,]
			.fmri.data.dims(dat)[3]=ny
		} 
		
		if(nz<z) {
			remz <- seq(nz+1,z,1)
			data_array <- data_array[,,-remz,]
			.fmri.data.dims(dat)[4]=nz
		} 
		
	
		if(!quiet) cat(.fmri.data.dims(dat),'\n')
		if(!quiet) cat('filename',.fmri.data.filename(dat),'\n')	
		if(!quiet) cat('datavec should be',nx*ny*nz*t,'long and is now',length(as.vector(data_array)),'\n')
		cat(paste('[crop] cropped file ',.fmri.data.filename(dat),' (',x,'x',y,'x',z,'x',t,') to: ',nx,'x',ny,'x',nz,'x',t,sep=''),'\n')
		if(file.exists(.fmri.data.filename(dat))) file.remove(.fmri.data.filename(dat))
		writeData(dat,as.vector(data_array))
	}
	
		
}


cropVolumeAuto <- 
function(betadir,quiet=T) 
#resize nifti images (crops last slices to destination dimensions)
{
	
	
	#determinesmallest
	filelist <- list.files(betadir,full.names=T)
	dims=matrix(NA,length(filelist),8)
	for(i in 1:length(filelist)) {
		
		dat = readData(filelist[i])
		dims[i,]=.fmri.data.dims(dat)
		
	}
	
	minx=min(dims[,2])
	miny=min(dims[,3])
	minz=min(dims[,4])
	
	for(i in 1:length(filelist)) cropVolume(filelist[i],c(1,minx,miny,minz),quiet)
	
	
}

readFSLmat <- 
function(filename)
#read FSL affine matrix file
{
	#set separator
	sp <- .Platform$file.sep
	
	#check if only file else pre-append working directory
	if(length(grep(sp,filename))==0) filename <- paste(getwd(),sp,filename,sep='')	
	
	#read in matrix
	mat = as.matrix(read.table(filename,header=F))
	
	return(invisible(mat))
}

makeROImask <-
function(fmridata,maskdata) 
#function to mask fmridata with non-zero elements in maskdata
{

	if(class(maskdata)=='fmri.data') whichzero = which(.fmri.data.datavec(maskdata)==0) else whichzero = which(maskdata==0)
	
	if(length(whichzero)>0) .fmri.data.datavec(fmridata)[whichzero]=0

	return(fmridata)

} 
 
createRegs <- 
function(arfdata) 
#createRegs creates registration files for each
{
	
	#set separator
	sp <- .Platform$file.sep
	
	#make new registration object
	registration <- new('registration')
	
	#get betafiles plus path of registration
	filelist <- .data.betafiles(arfdata)
	path <- .data.regDir(arfdata)
	
	#check betafile integrity and make paths in regDir
	for(filename in filelist) {
		if(file.exists(filename)) {
			
			#get info from betafile and set linkedfile
			info <- getFileInfo(filename)
			dirname <- .nifti.fileinfo.filename(info)
			
			#create dir and create regfilename
			if(!file.exists(paste(path,sp,dirname,sep=''))) dir.create(paste(path,sp,dirname,sep=''))
			
			.registration.linkedfile(registration) <- filename
			.registration.fullpath(registration) <- paste(path,sp,.nifti.fileinfo.filename(info),sep='')
			.registration.filename(registration) <- .data.regRda(arfdata)
			
			#save objects
			save(registration,file=paste(.registration.fullpath(registration),sp,.registration.filename(registration),sep=''))
						
		} else warning('No betafile found to match reg to')
	}
	#return(invisble(registration))
}


setRegFiles <- 
function(registration,examp2stand='example_func2standard.mat',examp2high='example_func2highres.mat',high2stand='highres2standard.mat',example_func='example_func.nii.gz',highres='highres.nii.gz',standard='standard.nii.gz') 
#setRegfiles fills the registration object with the correct matrices and nifti images (uses FSL standard as default)
{
	
	#set and check regFiles
	.registration.examp2stand(registration) <- examp2stand
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.examp2stand(registration),sep=''))) stop('ex2stand does not exist')
	.registration.examp2high(registration) <- examp2high
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.examp2high(registration),sep=''))) stop('ex2high does not exist')
	.registration.high2stand(registration) <- high2stand
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.high2stand(registration),sep=''))) stop('high2stand does not exist')
	.registration.example(registration) <- example_func
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.example(registration),sep=''))) stop('standard does not exist')
	.registration.highres(registration) <- highres
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.highres(registration),sep=''))) stop('highres does not exist')
	.registration.standard(registration) <- standard
	if(!file.exists(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.standard(registration),sep=''))) stop('standard does not exist')
	
	
	#save regfile object
	save(registration,file=paste(.registration.fullpath(registration),.Platform$file.sep,.registration.filename(registration),sep=''))
	
	#return registration object
	return(invisible(registration))
}


setRegParams <- 
function(registration) 
#setRegparams reads in registration parameters from the files in registration object and sets appropriate matrices
{
	#load registration volumes
	examp = readData(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.example(registration),sep=''))
	highres = readData(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.highres(registration),sep=''))
	standard = readData(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.standard(registration),sep=''))
	
	#Set pixdim matrices
	.registration.Dex(registration) = diag(c(.fmri.data.pixdim(examp)[2:4],1))
	.registration.Dhi(registration) = diag(c(.fmri.data.pixdim(highres)[2:4],1))
	.registration.Dst(registration) = diag(c(.fmri.data.pixdim(standard)[2:4],1))
	
	#set x-axis swap matrix
	.registration.SXhi(registration) = diag(c(-1,1,1,1))
	.registration.SXhi(registration)[1,4] = .fmri.data.dims(highres)[2]-1
	
	#set affine transformation matrices
	.registration.Aex2hi(registration) = readFSLmat(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.examp2high(registration),sep=''))
	.registration.Ahi2st(registration) = readFSLmat(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.high2stand(registration),sep=''))
	
	#set origin offset matrix
	.registration.OXst(registration) = rbind(.fmri.data.srow_x(standard),.fmri.data.srow_y(standard),.fmri.data.srow_z(standard),c(0,0,0,1))	
	
	#save object
	save(registration,file=paste(.registration.fullpath(registration),.Platform$file.sep,.registration.filename(registration),sep=''))
	
	return(invisible(registration))
	
}



arf2MNI <- 
#arfToMNI converts arf-native voxel locations to MNI_152 standard coordinates		
function(xyz_coor,registration) 
{
	#minus one to correct for vector in FSL starting at zero	
	xyz = c(xyz_coor-1,1)
	
	#examp_vox to mm
	examp_mm = .registration.Dex(registration)%*%xyz
	#cat('examp_mm',examp_mm,'\n')
	
	#examp_mm to high_mm
	high_mm = .registration.Aex2hi(registration)%*%examp_mm
	#cat('hiigh_mm',high_mm,'\n')
	
	#FLIP
	#high_mm to high_vox
	high_vox = solve(.registration.Dhi(registration))%*%high_mm
	#cat('high_vox',high_vox,'\n')	
	
	#x-axis flipped
	high_vox_flipped = .registration.SXhi(registration)%*%high_vox
	#cat('high_vox_flipped',high_vox_flipped,'\n')
	
	#high_vox to high_mm
	high_mm_flipped = .registration.Dhi(registration)%*%high_vox_flipped
	#cat('high_mm_flipped',high_mm_flipped,'\n')
	#END FLIP
	
	#high_mm to standard_mm
	stand_mm_flipped = .registration.Ahi2st(registration)%*%high_mm_flipped
	#cat('stand_mm_flipped',stand_mm_flipped,'\n')
	
	#standard_mm to standard_vox 
	stand_vox_flipped = solve(.registration.Dst(registration))%*%stand_mm_flipped
	#cat('stand_vox_flipped',stand_vox_flipped,'\n')
	
	#standard_vox to MNI (origin offset)
	stand_mni_flipped = .registration.OXst(registration)%*%stand_vox_flipped
	#cat('stand_mni_flipped',stand_mni_flipped,'\n')
	
	#return MNI
	return(stand_mni_flipped[-length(stand_mni_flipped)])
		
}

MNI2arf <- 
#arfToMNI converts MNI_152 standard coordinates	to arf-native voxel locations	
function(xyz_coor,registration) 
{
	#add one (at end to corrrect for arf coordinates starting at 1 in R)
	xyz = c(xyz_coor,1)

	#inverse standard_vox to MNI (origin offset)
	stand_mni_flipped = solve(.registration.OXst(registration))%*%xyz
	
	#inverse standard_mm to standard_vox 
	stand_vox_flipped = (.registration.Dst(registration))%*%stand_mni_flipped
	
	#inverse high_mm to standard_mm
	stand_mm_flipped = solve(.registration.Ahi2st(registration))%*%stand_vox_flipped
	
	#inverse high_vox to high_mm
	high_mm_flipped = solve(.registration.Dhi(registration))%*%stand_mm_flipped
	
	#FLIP
	#inverse x-axis flipped
	high_vox_flipped = solve(.registration.SXhi(registration))%*%high_mm_flipped
		
	#inverse high_mm to high_vox
	high_vox = (.registration.Dhi(registration))%*%high_vox_flipped
		
	#inverse examp_mm to high_mm
	high_mm = solve(.registration.Aex2hi(registration))%*%high_vox
	#END FLIP
	
	#inverse examp_vox to mm
	examp_mm = solve(.registration.Dex(registration))%*%high_mm
	
	#return arf
	return(examp_mm[-length(examp_mm)]+1)
	
}


arf2high <- 
#arfToMNI converts arf-native mm locations to highres subject scan.		
function(xyz_coor,registration) 
{
	#minus one to correct for vector in FSL starting at zero	
	xyz = c(xyz_coor-1,1)
	
	#examp_vox to mm
	examp_mm = .registration.Dex(registration)%*%xyz
	
	#examp_mm to high_mm
	high_mm = .registration.Aex2hi(registration)%*%examp_mm
	
	#return highres
	return(high_mm[-length(high_mm)])
	
}



MNI2atlas <- 
#MNI converts MNI_152 standard coordinates to atlas index coordinates 		
function(xyz_coor,registration) 
{
	#add one (at end to corrrect for arf coordinates starting at 1 in R)
	xyz = c(xyz_coor,1)
	
	#inverse standard_vox to MNI (origin offset)
	stand_mni_flipped = solve(.registration.OXst(registration))%*%xyz
			
	#return arf
	return(stand_mni_flipped [-length(stand_mni_flipped)]+1)
}

reqFlip <-
function(fmridata)
#check if a flip is required 
{
	flip = c(T,F,F)
	
	if(.fmri.data.sform_code(fmridata)>0) {
		
		if(.fmri.data.srow_x(fmridata)[1]<0) flip[1]=F
		if(.fmri.data.srow_y(fmridata)[2]<0) flip[2]=T
		if(.fmri.data.srow_z(fmridata)[3]<0) flip[3]=T
		
		
	} #else warning('sform_code not set, assuming radiological orientation')
	
	return(flip)
	
}

flipAxis <-
function(data_array,axis=c('x','y','z'))
#flips the axis of a data_array
{
	#which axis
	axis = match.arg(axis)
	
	#create standard affine matrix
	new_dat = data_array
	affine = diag(c(1,1,1,1))
	
	#flip x-axis
	if(axis=='x') {
		affine[1,1] = -1
		affine[1,4] = dim(data_array)[1]+1
		for(x in 1:dim(data_array)[1]) {
			nc = affine%*%c(x,1,1,1)
			new_dat[nc[1],,] = data_array[x,,]
		}
	} 
	
	#flip y-axis
	if(axis=='y') {
		affine[2,2] = -1
		affine[2,4] = dim(data_array)[2]+1
		for(y in 1:dim(data_array)[2]) {
			nc = affine%*%c(1,y,1,1)
			new_dat[,nc[2],] = data_array[,y,]
		}
	} 
	
	#flip z-axis
	if(axis=='z') {
		affine[3,3] = -1
		affine[3,4] = dim(data_array)[3]+1
		for(z in 1:dim(data_array)[3]) {
			nc = affine%*%c(1,1,z,1)
			new_dat[,,nc[3]] = data_array[,,z]
		}
	} 
	
	return(new_dat)	
}




euclidDist <-
function(arfmodel,thres=5,quiet=T) 
#calculate euclidian distances between estimates
{
	
	theta = matrix(.model.estimates(arfmodel),10)
	distmat = matrix(0,ncol(theta),ncol(theta))
	
	p=0
	#for all off-diagonal elements calculate Euclidian distance
	for(i in 1:(ncol(theta)-1)) {
		for(j in (i+1):ncol(theta)) {
			
			dist=0
			for(k in 1:3) dist = dist + (theta[k,i]-theta[k,j])^2
			distmat[i,j]=distmat[j,i]=c(sqrt(dist))
			
			opp=FALSE
			if(theta[10,j]<0 & theta[10,i]>0) opp=TRUE
			if(theta[10,j]>0 & theta[10,i]<0) opp=TRUE
			
			if(i!=j & distmat[i,j]<thres & opp) {
				if(!quiet) cat('region',i,'[',round(theta[c(1,2,3,4,5,6,10),i]),'] -',j,'[',round(theta[c(1,2,3,4,5,6,10),j]),'] distance:',sqrt(dist),'\n')
			}
			p=p+1	
		}
	}
	return(distmat)
}


makeLowResStruct <-
function(arfdata,experiment=NULL)
#make low resolution structural image from high_res T1 image
{
	
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#get runs from dataDir
	runs = list.files(.data.regDir(arfdata),full.names=F)
	
	if(length(runs)==0) stop('No runs directories found found in',.data.regDir(arfdata),',please run registration.')
	
	for(rundir in runs) {
		
		#load registration parameters
		registration = loadRda(paste(.data.regDir(arfdata),.Platform$file.sep,rundir,.Platform$file.sep,.data.regRda(arfdata),sep=''))
		
		#read in lowres + highresfiles
		examp = readData(.registration.linkedfile(registration))
		highres = readData(paste(.registration.fullpath(registration),.Platform$file.sep,.registration.highres(registration),sep=''))
		
		
		#set dimensions of highres image
		highdat = .fmri.data.datavec(highres)
		dim(highdat) = c(.fmri.data.dims(highres)[2],.fmri.data.dims(highres)[3],.fmri.data.dims(highres)[4])
		
		#set dimensions of lowres image
		dimx = .fmri.data.dims(examp)[2]
		dimy = .fmri.data.dims(examp)[3]
		dimz = .fmri.data.dims(examp)[4]
		
		#make newfile and fill with zeroes
		newdat = examp
		.fmri.data.filename(newdat) = .experiment.lowresFile(experiment)
		.fmri.data.fullpath(newdat) = paste(.data.regDir(arfdata),.Platform$file.sep,rundir,sep='')
		newdatavec = rep(0,dimx*dimy*dimz)
		dim(newdatavec) = c(dimx,dimy,dimz)		
		
		#transform highresimage to lowres
		for(z in 1:dimz) {
			for(y in 1:dimy) {
				for(x in 1:dimx) {
							
					xyz = c(x,y,z,1)
					examp_mm = .registration.Dex(registration)%*%xyz
					high_mm = .registration.Aex2hi(registration)%*%examp_mm
					high_vox = solve(.registration.Dhi(registration))%*%high_mm
					newdatavec[x,y,z] = highdat[high_vox[1],high_vox[2],high_vox[3]] 

				}
			}
		}
		
		.fmri.data.cal_min(newdat) = min(newdatavec)
		.fmri.data.cal_max(newdat) = max(newdatavec)
		
		#write data to new image
		invisible(writeData(newdat,as.vector(newdatavec)))
		
	}
	
}

makeLowResStructAvg <-
function(arfmodel,experiment=NULL)
#make average of low_resolution structural image from multiple lowres images (and save in modeldir)
{
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	sp = .Platform$file.sep
	
	#get run list
	runs = list.files(.model.regDir(arfmodel),full.names=F)
	
	#set avg to zero
	avgdat = 0
		
	#load and sum images
	for(rundir in runs) {
		registration = loadRda(paste(.model.regDir(arfmodel),sp,rundir,sp,.model.regRda(arfmodel),sep=''))
		fn = list.files(path=.registration.fullpath(registration),pattern=.experiment.lowresFile(experiment),full.names=T)
		lrdat = readData(fn[1])
		avgdat = avgdat + .fmri.data.datavec(lrdat)
	}
	avgdat = avgdat / length(runs)
	
	#set parameters for newfile and save to modelpath
	avgstruct = lrdat
	.fmri.data.filename(avgstruct) = .experiment.lowresAvg(experiment)
	.fmri.data.fullpath(avgstruct) = .model.modeldatapath(arfmodel)
	.fmri.data.cal_min(avgstruct) = min(avgdat)
	.fmri.data.cal_max(avgstruct) = max(avgdat)
	invisible(writeData(avgstruct,avgdat))

}


getTalairach <-
function(FSLDIR='/usr/local/fsl',which='2mm')
#get talairach indices from FSL talairach file
{
	sp <- .Platform$file.sep
	
	
	#set atlaspath of FSL
	atlaspath = paste(FSLDIR,sp,'data/atlases',sep='')
	taldat = read.table(paste(atlaspath,sp,'talairach.xml',sep=''),fill=T,skip=15,sep='\n',stringsAsFactors=F)
	
	talmat = matrix(NA,1106,6)
	for(i in 1:1106) {
	
		#get index and names
		splits = strsplit(strsplit(taldat[i,1],'>')[[1]],'<')
		index = as.numeric(strsplit(strsplit(splits[[1]],' ')[[2]],'=')[[2]][2])
		namevec = strsplit(splits[[2]][1],"\\.")[[1]]
		
		talmat[i,1]=index
		talmat[i,2:6]=namevec
		
	}
	
	#set 2mm or 1mm file
	which = match.arg(which)
	if(which=='2mm') niiname = 'talairach-labels-2mm.nii.gz'
	if(which=='1mm') niiname = 'talairach-labels-1mm.nii.gz'
	
	#read in file and make array (MNI voxel space)
	talname = paste(FSLDIR,sp,'data/atlases/talairach',sp,niiname,sep='')
	talfile = readData(talname)
	talmap = .fmri.data.datavec(talfile)
	dim(talmap) = c(.fmri.data.dims(talfile)[2],.fmri.data.dims(talfile)[3],.fmri.data.dims(talfile)[4])
	
	#make atlas class
	talairach <- list(data=talmat,map=talmap)
	attr(talairach,'class') <- 'TalairachAtlas'
	
	return(talairach)
		
}


getHarvardOxford <-
function(FSLDIR='/usr/local/fsl',which='2mm')
#get harvard-oxford probability map (cortical and subcortical)
{
	sp <- .Platform$file.sep
		
	#get cortical and subcortical data
	atlaspath = paste(FSLDIR,sp,'data/atlases',sep='')
	cortdat  = read.table(paste(atlaspath,sp,'HarvardOxford-Cortical.xml',sep=''),fill=T,skip=16,sep='\n',stringsAsFactors=F,quote="")
	cortdat = apply(cortdat,2,function(x) gsub("\"",'',x))
	
	cortical = matrix(NA,48,2)
	for(i in 1:48) {
		
		#get index and names
		splits = strsplit(strsplit(cortdat[i,1],'>')[[1]],'<')
		index = as.numeric(strsplit(strsplit(splits[[1]],' ')[[2]],'=')[[2]][2])
		namevec = strsplit(splits[[2]][1],"\\.")[[1]]
		
		cortical[i,1]=index
		cortical[i,2]=namevec
		
	}
	
	subdat = read.table(paste(atlaspath,sp,'HarvardOxford-Subcortical.xml',sep=''),fill=T,skip=16,sep='\n',stringsAsFactors=F,quote="")
	subdat = apply(subdat,2,function(x) gsub("\"",'',x))
	
	subcort = matrix(NA,21,2)
	for(i in 1:21) {
		
		#get index and names
		splits = strsplit(strsplit(subdat[i,1],'>')[[1]],'<')
		index = as.numeric(strsplit(strsplit(splits[[1]],' ')[[2]],'=')[[2]][2])
		namevec = strsplit(splits[[2]][1],"\\.")[[1]]
		
		subcort[i,1]=index
		subcort[i,2]=namevec
		
	}
	
	which = match.arg(which)
	if(which=='2mm') {
		cortnii = 'HarvardOxford-cort-prob-2mm.nii.gz'
		subnii = 'HarvardOxford-sub-prob-2mm.nii.gz'
	}
	
	if(which=='1mm') {
		cortnii = 'HarvardOxford-cort-prob-1mm.nii.gz'
		subnii = 'HarvardOxford-sub-prob-1mm.nii.gz'
	}
	
	#load cortical probability volume
	cortfile = paste(FSLDIR,sp,'data/atlases/HarvardOxford',sp,cortnii,sep='')
	corticalvol = readData(cortfile)
	corticalmap = .fmri.data.datavec(corticalvol)
	dim(corticalmap) = c(.fmri.data.dims(corticalvol)[2],.fmri.data.dims(corticalvol)[3],.fmri.data.dims(corticalvol)[4],.fmri.data.dims(corticalvol)[5])
	
	#load subcortical probability volume
	subfile = paste(FSLDIR,sp,'data/atlases/HarvardOxford',sp,subnii,sep='')
	subcortvol = readData(subfile)
	subcortmap = .fmri.data.datavec(subcortvol)
	dim(subcortmap) = c(.fmri.data.dims(subcortvol)[2],.fmri.data.dims(subcortvol)[3],.fmri.data.dims(subcortvol)[4],.fmri.data.dims(subcortvol)[5])
	
	#make atlas class
	harvardoxford = list(cortical=list(data=cortical,map=corticalmap),subcortical=list(data=subcort,map=subcortmap))
	attr(harvardoxford,'class') <- 'HarvardOxfordAtlas'
	
	return(harvardoxford)
	
}


getAtlasLabels <-
function(coordinates,registration,coortype=c('arf','mni'),atlas=c('both','Talairach','HarvardOxford'),...)
#get atlas labels for a set of coordinates using a registration class and atlas type
{
	if(!is.matrix(coordinates)) stop('input coordinates must be in a matrix')
	if(ncol(coordinates)!=3) stop('input coordinate matrix must have three columns (x,y,z)')
	
	#get atlases
	atlas = match.arg(atlas)
	if(atlas=='both') {
		talairach = getTalairach(...)
		harvard = getHarvardOxford(...)
	}
	if(atlas=='Talairach') talairach = getTalairach(...)
	if(atlas=='HarvardOxford') harvard = getHarvardOxford(...)
	
	#convert coordinates to MNI voxel space
	coortype = match.arg(coortype)
	if(coortype=='arf') coordinates = t(apply(coordinates,1,function(coordinates,registration) arf2MNI(coordinates,registration),registration=registration))
	coordinates = t(round(apply(coordinates,1,function(coordinates,registration) MNI2atlas(coordinates,registration),registration=registration)))
	
	#get labels for the coordinates
	atlasvec = vector(nrow(coordinates),mode='list')
	for(i in 1:nrow(coordinates)) {
		if(atlas=='both') atlasvec[[i]] = list(harvard=calcHarvardOxfordProbs(coordinates[i,],harvard),talairach=calcTalairach(coordinates[i,],talairach))
		if(atlas=='Talairach') atlasvec[[i]] = list(talairach=calcTalairach(coordinates[i,],talairach))
		if(atlas=='HarvardOxford') atlasvec[[i]] = list(harvard=calcHarvardOxfordProbs(coordinates[i,],harvard))
		
	}
	
	return(atlasvec)
}

calcHarvardOxfordProbs <-
function(xyz,harvard)
#calculate harvardOxford probabilities, return vector with percentage labels
{
	
	#check validity of coordinates
	valid = TRUE
	if(xyz[1]<1 | xyz[1]>dim(harvard$cortical$map)[1]) valid=FALSE
	if(xyz[2]<1 | xyz[2]>dim(harvard$cortical$map)[2]) valid=FALSE
	if(xyz[3]<1 | xyz[3]>dim(harvard$cortical$map)[3]) valid=FALSE
	
	labelvec = character(0)
	
	if(valid) {
		#get cortical and subcortical data
		cortical_probvec = harvard$cortical$map[xyz[1],xyz[2],xyz[3],] 
		subcort_probvec = harvard$subcortical$map[xyz[1],xyz[2],xyz[3],] 
		
		cort_large_probs = which(cortical_probvec>=1)
		sub_large_probs = which(subcort_probvec>=1)
	
		for(i in 1:length(cort_large_probs)) {
			index = cort_large_probs[i]
			labelvec = c(labelvec,paste(cortical_probvec[cort_large_probs[i]],'% ',harvard$cortical$data[index,2],sep=''))
		}
		
		for(i in 1:length(sub_large_probs)) {
			index = sub_large_probs[i]
			labelvec = c(labelvec,paste(subcort_probvec[sub_large_probs[i]],'% ',harvard$subcortical$data[index,2],sep=''))
		}
		
	}
	
	return(labelvec)
	
}

calcTalairach <-
function(xyz,talairach)
#calculate talairach index labels, return vector with matching labels
{
	#check validity of coordinates
	valid = TRUE
	if(xyz[1]<1 | xyz[1]>dim(talairach$map)[1]) valid=FALSE
	if(xyz[2]<1 | xyz[2]>dim(talairach$map)[2]) valid=FALSE
	if(xyz[3]<1 | xyz[3]>dim(talairach$map)[3]) valid=FALSE
	
	labelvec = character(0)
	
	if(valid) {
		#get talairach labels
		index = talairach$map[xyz[1],xyz[2],xyz[3]]
		labelvec = c(labelvec,talairach$data[index+1,2:6])
	}
	
	return(labelvec)
}

getModelAtlas <-
function(arfmodel,regrun=1,saveastext=F) 
#wrapper for Atlas coordinates of model estimates
{
	sp <- .Platform$file.sep
	runs = list.files(.model.regDir(arfmodel),full.names=F)[regrun]
	registration = loadRda(paste(.model.regDir(arfmodel),sp,runs,sp,.model.regRda(arfmodel),sep=''))
	
	estimates = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
	coors = t(round(estimates[c(1,2,3),]))
	
	atlas = getAtlasLabels(coors,registration,coortype='arf')
	
	if(saveastext) sink(paste(.model.modelname(arfmodel),'-atlas.txt',sep=''),split=T)
	
	for(i in 1:nrow(coors)) {
		cat('Region ',i,' [',coors[i,1],' ',coors[i,2],' ',coors[i,3],'] @ ',round(estimates[10,i]),'\n',sep='')
		cat('HarvardOxford:\n')
		for(j in 1:length(atlas[[i]]$harvard)) cat('  ',atlas[[i]]$harvard[j],'\n')
		cat('Talairach:\n')
		for(j in 1:length(atlas[[i]]$talairach)) cat('  ',atlas[[i]]$talairach[j],'\n')
		cat('\n')	
	}
	
	if(saveastext) sink()

	return(invisible(atlas))
	
}


highresBlobs <-
function(arfmodel,registration) 
#make highresolution blobs (with outlines at 95% spat.ext)
{
	
	ests = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
	
	centers = matrix(NA,.model.regions(arfmodel),3)
	
	for(reg in 1:ncol(ests)) {
		#minus one to correct for vector in FSL starting at zero	
		xyz = c(ests[c(1,2,3),reg]-1,1)
		
		#examp_vox to mm
		examp_mm = .registration.Dex(registration)%*%xyz
		
		#examp_mm to high_mm
		high_mm = .registration.Aex2hi(registration)%*%examp_mm
		
		#high_mm to high_vox
		high_vox = solve(.registration.Dhi(registration))%*%high_mm
	
		#set center locations in highres vox coordinates
		centers[reg,] = round(high_vox[-length(high_vox)])-1	
		
		#calculate sigma matrix
		sigma = matrix(c(ests[4,reg]^2,ests[4,reg]*ests[5,reg]*ests[7,reg],ests[4,reg]*ests[6,reg]*ests[8,reg],ests[4,reg]*ests[5,reg]*ests[7,reg],ests[5,reg]^2,ests[5,reg]*ests[6,reg]*ests[9,reg],ests[4,reg]*ests[6,reg]*ests[8,reg],ests[5,reg]*ests[6,reg]*ests[9,reg],ests[6,reg]^2),3)
		
	}
	return(centers)
}

loadReg <-
function(subject, condition, run, experiment = NULL)
{
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#set separator
	sp <- .Platform$file.sep
	cpath <- paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sp,subject,sp,.experiment.conditionDir(experiment),sp,condition,sp,.experiment.dataDir(experiment),sp,.experiment.regDir(experiment),sep='')
	
	regdirs <- list.files(cpath,full.names=T)
	regdirs <- regdirs[which(file.info(list.files(cpath,full.names=T))$isdir)]
	
	reg=NULL
	
	#check run
	if(is.numeric(run)) reg = loadRda(paste(regdirs[run],sp,.experiment.regRda(experiment),sep=''))
	if(is.character(run)) reg = loadRda(paste(cpath,sp,run,sp,.experiment.regRda(experiment),sep=''))
	
	return(reg)
	
}

saveReg <-
function(registration)
{
	save(registration,file=paste(.registration.fullpath(registration),.Platform$file.sep,.registration.filename(registration),sep=''))
	
}
