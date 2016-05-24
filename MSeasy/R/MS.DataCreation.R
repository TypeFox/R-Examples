MS.DataCreation <-
function(DataType="CDF", path="", pathCDF="", mz, N_filt=3, apex=FALSE, quant=FALSE)
{
attach(what=NULL, name="e1")
### MSeasy v 1.3 March 2011 with improved errorhandling and quant option
if (DataType=="Agilent"||DataType=="ASCII")
{
	#a modifier
if (DataType=="Agilent")
{
	st<-strsplit(date(), " ")[[1]]
	stBis<-strsplit(st[4], ":")[[1]]
	Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
	Date<-paste(st[1], st[2], st[3], Hour, sep="_")
	Mypath<-paste("output_MSDataCreation", "_", "result", Date, sep="")
	dir.create(Mypath)

	if ((path!="")==TRUE){
		path<-path
		print(path)
	}
	if ((path!="")==FALSE){
		require("tcltk")
		path=tclvalue(tkchooseDirectory(title="Please, select a directory containing the .D folders"))
	}
	#list rteres.txt files
	list_rteres<-dir(path, pattern="rteres.txt", recursive =TRUE)
	print(list_rteres)

	##if exist delete save_list_temp.rda & initial_DATA.txt
	unlink(paste(Mypath, "/","save_list_temp.rda", sep=""))
	unlink(paste(Mypath, "/","initial_DATA.txt", sep=""))
	
	#All CDF (mzXML, mzML or mzDATA) files must be in the same directory path1
	if ((pathCDF!="")==TRUE){
	path1<-pathCDF
	print(pathCDF)
	}
	if ((pathCDF!="")==FALSE){
	require("tcltk")
	path1=tclvalue(tkchooseDirectory(title="Please, select your CDF directory"))
	}
	filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
						 "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
	filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),
	collapse = "|") 
	#files is a liste of all CDF files found in path1 directory
	files<-list.files(path=path1, pattern=filepattern, recursive = TRUE, full.names=TRUE) #full.names conserve the full path to the file
	list_data<-files
	print(files)
	
	#set mz value to all if not entered
	if (missing(mz)==TRUE || is.null(mz)){
		mz="all"
	}
	
	#specific variables for CDF error handling
	mz_min <- vector()
	mz_max<-vector()
	pbmzmin<-0
	pbmzmax<-0

	#other variables
	rteres<-list()
	Data<-list()
	data_final<-list()
	an_rteres<-vector()
	an_data<-vector()
	an<-vector()
	problem1<-0
	problem2<-0
	pblist<-vector()
	pblistl<-vector()
	pb<-1

	### keep the firsts characters in the file names of the .d files
	# Start change specific for CDF files
	   for (l in 1:length(list_rteres)) an_rteres[l] <- strsplit(list_rteres[l], 
				"/")
			for (l in 1:length(list_data)) an_data[l] <- strsplit(list_data[l], 
				"/")
			an_rteres <- unlist(lapply(an_rteres, "[", 1))
			an_data <- unlist(lapply(an_data, "[", length(an_data[[1]])))
				   

			an_rteres<-strsplit(an_rteres,".D", fixed=TRUE)
			an_data<-strsplit(an_data,filepattern)
			an_rteres <- unlist(lapply(an_rteres, "[", 1))
			an_data <- unlist(lapply(an_data, "[", 1))

		### Check that the list of names in CDF files and rteres files are the same

		test_pb1<-match(an_rteres, an_data)
		test_pb2<-match(an_data, an_rteres)

		   for (p in 1:length(test_pb1))
					 if (is.na(test_pb1[p])==TRUE)
					 {
					 problem1<-1
					 pb1<-p
					 }

		   for (k in 1:length(test_pb2))
					 if (is.na(test_pb2[k])==TRUE)
					 {
					 problem2<-1
					 pb2<-k
					 }

		   if (problem1==1)
				 cat(" There is a problem in the rteres for ", an_rteres[pb1])

		   if (problem2==1)
		 #start change	  
				 #cat("There is a problem in the export3ddata for ", an_data[pb2])
					cat(" There is a problem in the CDF (mzXML, mzML or mzDATA) file for ", 
						an_data[pb2])
		#end change	
		   if (problem1==0)
				 if (problem2==0)
				 {
				 an<-an_rteres
				 print(an)
				 }

		##Search for starting point in rteres.txt files
		skip_value<-vector()
		for (m in 1:length(an))
		{
		ma_test<-as.matrix(readLines(paste(paste(path,"/", sep=""), list_rteres[m], sep=""),n=-1))
		skip_value[m]<-as.numeric(agrep("scan scan scan",ma_test))
		}
		#search for differences in user defined mz vector and reel mz present in CDF files


cdf_error<-function(){

	for (m in 1:length(an)) {
			obj<-xcms::xcmsRaw(files[m]) 
			mz_min[m]<-min(obj@mzrange)
			mz_max[m]<-max(obj@mzrange)
		}
		#case mz="all"
		if (length(mz)==1 && mz=="all"){
					mz<-max(mz_min):min(mz_max)
		}
		#case mz entered by user or estimated after mz="all"
		if (max(mz_min)>mz[1]){
		
					mz<-max(mz_min):mz[length(mz)]
					pbmzmin<-1
		}	
		if (min(mz_max)<mz[length(mz)]){
	
					mz<-mz[1]:min(mz_max)
					pbmzmax<-1
				}
		#error handling
		if (pbmzmin!=0){
					
					print(paste("mz minimum value was set to ", mz[1], cat("\n")))
					pbmzmin<-0
				}
		if (pbmzmax!=0){
					print(paste("mz maximum value was set to ", min(mz_max), cat("\n")))
					pbmzmax<-0
			}
	assign("mz", mz, envir=as.environment("e1"))
	return(mz)
}

mz<-cdf_error()

print(paste("mz:",mz[1],":",mz[length(mz)],cat("\n")))

##Search for errors in rteres.txt

`errorrteres`<-function(){
       Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE)
       Rte<-Rte[1:(dim(Rte)[1]-2),]
# added Y GUITTON 2012-05-17
		Rte[,9]<-as.vector(Rte[,9])
		Rte[,10]<-as.vector(Rte[,10])
		Rte[,11]<-as.vector(Rte[,11])

		for (co in 1:nrow(Rte))
			if (Rte[co, 11]=="")
		{
		Rte[co,11]<-Rte[co,10]
		Rte[co,10]<-Rte[co,9]
		Rte[co,9]<-Rte[co,8]
		Rte[co,8]<-Rte[co,7]
		Rte[co,7]<-1
		}
##end of addition

       Rte<-Rte[,c(1:5, 9, 11)]
       colnames(Rte)<-c("peak", "RT", "first_scan", "max_scan", "last_scan", "quantification1", "quantification2")
## Treat CDF files 3d data
	   ###start change
			obj<-xcms::xcmsRaw(files[m]) 
				Td<-t(obj@env$profile)
				colnames(Td)<-c(min(obj@mzrange):max(obj@mzrange))
				
				if (is.na(match(mz[length(mz)],colnames(Td))) == TRUE){
					Td<-Td[,match(mz[1],colnames(Td)):match(max(obj@mzrange),colnames(Td))]
					
				}else{
					Td<-Td[,match(mz[1],colnames(Td)):match(mz[length(mz)],colnames(Td))]
				}
				colnames(Td)<-c(mz[1]:mz[length(mz)])
				
			End <- match(mz[length(mz)], colnames(Td))
            Start <- match(mz[1], colnames(Td))	
			a_env <- vector()
                av <- 1
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$last_scan))[pic] > 
                  as.numeric(dim(Td)[1])) {# number of the last scan in the file profile is a matrix with col=mz and row =scan number
                  a_env[av] <- pic
                  av <- av + 1
                }
#end of change	
## Treat export3d files removed

}#end test error

for (m in 1:length(an)){result<-try(errorrteres(), silent=TRUE); if(class(result) == "try-error"){pblist[pb]<-list_rteres[m]; pb=pb+1; assign("pblist",pblist ,envir=as.environment("e1"));next;}}
#for (m in 1:length(an)){result<-try(errorrteres(), silent=TRUE); if(class(result) == "try-error"){pblist[pb]<-list_rteres[m]; pb=pb+1; return(pblist);next;}}
pblistl<-length(pblist)


if (pblistl!=0){

       eval(cat("A problem occured in the rteres.txt files listed below \n Check the pkty column for missing value for the first peak \n"), envir=as.environment("e1")); pblist<-data.frame(pblist); eval(print(pblist), envir=as.environment("e1"));
       stop
       }
else{
for (m in 1:length(an)){## rteres files treatment


Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE)
Rte<-Rte[1:(dim(Rte)[1]-2),]
Rte[,9]<-as.vector(Rte[,9])
Rte[,10]<-as.vector(Rte[,10])
Rte[,11]<-as.vector(Rte[,11])

for (co in 1:nrow(Rte))
if (Rte[co, 11]=="")
{
Rte[co,11]<-Rte[co,10]
Rte[co,10]<-Rte[co,9]
Rte[co,9]<-Rte[co,8]
Rte[co,8]<-Rte[co,7]
Rte[co,7]<-1
}

Rte<-Rte[,c(1:5, 9, 11)]
colnames(Rte)<-c("peak", "RT", "first_scan", "max_scan", "last_scan", "quantification1", "quantification2")

## Treat export3d files removed from here

#start change
              
                obj<-xcms::xcmsRaw(files[m]) 
				Td<-t(obj@env$profile)
				colnames(Td)<-c(min(obj@mzrange):max(obj@mzrange))
				if (is.na(match(mz[length(mz)],colnames(Td))) == TRUE){
					Td<-Td[,match(mz[1],colnames(Td)):match(max(obj@mzrange),colnames(Td))]
					
				}else{
					Td<-Td[,match(mz[1],colnames(Td)):match(mz[length(mz)],colnames(Td))]
				}
				colnames(Td)<-c(mz[1]:mz[length(mz)])
				
				
                End <- match(mz[length(mz)], colnames(Td))
                Start <- match(mz[1], colnames(Td))
                Td <- Td[, Start:End] #reduction of Td dimension to the only mz needed


## Remove incomplete peaks that are at the beginning or at the end of the chromatogram, i.e. with missing scans in the Export3d

				a_env <- vector()
                av <- 1
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$last_scan))[pic] > 
                  as.numeric(dim(Td)[1])) {# number of the last scan in the file profile is a matrix with col=mz and row =scan number
                  a_env[av] <- pic
                  av <- av + 1
                }
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$first_scan))[pic] < 
                  1) { #the first scan of the file IS 1 (I hope so)
                  a_env[av] <- pic
                  av <- av + 1
                }
#end CDF specific change 

if (length(a_env)!=0)
Rte<-Rte[-a_env,]

## Each part of the list contains a matrix where each row corresponds to a given peak in the analysis

if (quant ==TRUE)
{
#changes here
#tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td)+2)
tempo <- matrix(nrow = nrow(Rte), ncol = length(mz)+3)
colnames(tempo)<-c("RT", "quantification1", "quantification2", mz)
tempo[,1]<-as.numeric(as.vector(Rte$RT))
tempo[,2]<-as.numeric(as.vector(Rte$quantification1))
tempo[,3]<-as.numeric(as.vector(substr(Rte$quantification2, 1, 4)))

  for (Pi in 1:nrow(Rte))
  {
  print(c(m, Pi))

## If apex=TRUE, look for scan_max

           if (apex==TRUE)
           {
           sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
 # start change
              #a <- as.matrix(Td[Td$scan_number == sc, 2:ncol(Td)])
              a<-as.matrix(Td[sc,])

		   Maxi<-max(a)
           a<-(a/Maxi)*100

           tempo[Pi,4:ncol(tempo)]<-as.vector(a)
           }

## If apex=FALSE, look for scanmin and scanmax for each peak and select 5% scan around scan_max

           if (apex==FALSE)
           {
           Rang<-as.numeric(as.vector(Rte$last_scan))[Pi]-as.numeric(as.vector(Rte$first_scan))[Pi]
           sc5<-Rang*5/100
           sc<-round(as.numeric(as.vector(Rte$max_scan))[Pi]-sc5, digits=0):round(as.numeric(as.vector(Rte$max_scan))[Pi]+sc5, digits=0)
			#a <- Td[Td$scan_number == sc[1], 2:ncol(Td)]
					a<-Td[sc[1],]

                   if (length(sc)>1)
                   {
                       for (La in 2:length(sc))
                       {
                       #a <- rbind(a, Td[Td$scan_number == sc[La], 2:ncol(Td)])
						a <- rbind(a, Td[sc[La],])
                       }

                       a<-as.matrix(a)

                       for (p in 1:nrow(a))
                       {
                        Maxi <- max(a[p, 1:ncol(a)]) #changed 2 in 1:ncol
                       a[p,]<-(a[p,]/Maxi)*100
                       }

                       tempo[Pi,4:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                   }

                   else
                   {
                   sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
					#a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
				   a <- as.matrix(Td[sc, ])
                   Maxi<-max(a)
                   a<-(a/Maxi)*100
                   tempo[Pi,4:ncol(tempo)]<-as.vector(a)
                   }
           }

data_final[[m]]<-tempo
}
}

if (quant ==FALSE)
{
#tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td))
tempo <- matrix(nrow = nrow(Rte), ncol = length(mz)+1)
colnames(tempo)<-c("RT", mz)
tempo[,1]<-as.numeric(as.vector(Rte$RT))



  for (Pi in 1:nrow(Rte))
  {
  print(c(m, Pi))

## If apex=TRUE, look for scan_max

           if (apex==TRUE)
           {
           sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
           #a <- as.matrix(Td[Td$scan_number == sc, 2:ncol(Td)])
            a<-as.matrix(Td[sc,])
           Maxi<-max(a)
           a<-(a/Maxi)*100

           tempo[Pi,2:ncol(tempo)]<-as.vector(a)
           }

## If apex=FALSE, look for scanmin and scanmax for each peak and select 5% scan around scan_max

           if (apex==FALSE)
           {
           Rang<-as.numeric(as.vector(Rte$last_scan))[Pi]-as.numeric(as.vector(Rte$first_scan))[Pi]
           sc5<-Rang*5/100
           sc<-round(as.numeric(as.vector(Rte$max_scan))[Pi]-sc5, digits=0):round(as.numeric(as.vector(Rte$max_scan))[Pi]+sc5, digits=0)
           
#a <- Td[Td$scan_number == sc[1], 2:ncol(Td)]
					 a<-Td[sc[1],]
                   if (length(sc)>1)
                   {
                       for (La in 2:length(sc))
                       {
                        #a <- rbind(a, Td[Td$scan_number == sc[La], 2:ncol(Td)])
						a <- rbind(a, Td[sc[La],])
                       }

                       a<-as.matrix(a)

                       for (p in 1:nrow(a))
                       {
                       Maxi <- max(a[p, 1:ncol(a)]) #changed 2 in 1:ncol
                       a[p,]<-(a[p,]/Maxi)*100
                       }

                       tempo[Pi,2:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                   }

                   else
                   {
                   sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
                   #a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
				   a <- as.matrix(Td[sc, ])
                   Maxi<-max(a)
                   a<-(a/Maxi)*100
                   tempo[Pi,2:ncol(tempo)]<-as.vector(a)
                   }
           }

data_final[[m]]<-tempo
}
}

#save(data_final, file=paste(Mypath, "/","save_list_temp.rda", sep=""))
save(data_final, file=file.path(Mypath,"save_list_temp.rda"))
}#end for

names(data_final)<-an

## Finally all informations for each analysis are grouped in a unique matrix

  for (i in 1:length(data_final))
  {
  temp<-rep(names(data_final)[i], nrow(data_final[[i]]))
  data_final[[i]]<-cbind(temp, data_final[[i]])
  }

data_fin<-data_final[[1]]

 if (length(data_final)>1)
    for (j in 2:length(data_final))
    {
    data_fin<-rbind(data_fin, data_final[[j]])
    }

data_fin<-as.data.frame(data_fin)

 for (i in 3:ncol(data_fin))
 {
 data_fin[,i]<-as.numeric(as.vector(data_fin[,i]))
 }

colnames(data_fin)[1]<-"analysis"


if (quant == TRUE)
{
colnames(data_fin)<-c("analysis", "retention_time", "quantification1", "quantification2", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}
if (quant == FALSE)
{
colnames(data_fin)<-c("analysis", "retention_time", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}



#Rprof(NULL)
#summaryRprof(filename = "Rprof.out")$sampling.time
print(paste("A data file has been generated in the folder:", Mypath, cat("\n")))
return(data_fin)
}#end else


}#end if AGILENT


if (DataType=="ASCII" )

## this is for ASCII files

{

st<-strsplit(date(), " ")[[1]]
	stBis<-strsplit(st[4], ":")[[1]]
	Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
	Date<-paste(st[1], st[2], st[3], Hour, sep="_")
	Mypath<-paste("output_MSDataCreation", "_", "result", Date, sep="")
	dir.create(Mypath)

	if ((path!="")==TRUE){
		path<-path
		print(path)
	}
	if ((path!="")==FALSE){
		require("tcltk")
		path=tclvalue(tkchooseDirectory(title="Please, select your directory containing ASCII files"))
	}

L_result<-list()
L_an<-dir(path)
res_data_temp<-list()
print(L_an)

for (a in 1:length(L_an))
{

#change to file.path instead of paste(path,"/",L_an[a], sep="")
#add check.names=FALSE for mz values
te<-read.table(file.path(path, L_an[a]), header=TRUE, check.names=FALSE)
#check if mz is missing if DataType=ASCII 
if(missing(mz)==TRUE){
	mzall<-colnames(te)[3]:colnames(te)[dim(te)[2]]
	mz<-mzall
	print(paste("mz = ",colnames(te)[3],":",colnames(te)[dim(te)[2]],sep=""))
}
#check if mz values entered by user are correct
if(missing(mz)==FALSE){
	mzall<-colnames(te)[3]:colnames(te)[dim(te)[2]]
	if (min(mz)<min(mzall)){
		mz=min(mzall):max(mz)
		print(paste("mz min changed by: ", min(mzall), sep=""))
	}
	if (max(mzall)<max(mz)){
		mz=mz[1]:max(mzall)
		print(paste("mz max changed by: ", max(mzall), sep=""))
	}
	print(paste("mz = ",mz[1],":",mz[length(mz)],sep=""))
}
result<-vector()

      if (N_filt>3)
      {
      ## If N_filt>3, smoothing of the chromatogram

      filt<-filter(te[,2], c(rep(1, N_filt))/N_filt, method="convolution")
      b<-cbind(te[,1],filt)

      ## detection of peaks in the chromatogram

      for (i in N_filt:(dim(b)[1]-N_filt))
      {
          if (b[i,2]>5000 & (b[i,2]>b[(i-1),2]) & (b[i,2]>b[(i+1),2]) & (b[(i-1),2]>b[(i-2),2]) & (b[(i+1),2]>b[(i+2),2]) & (b[(i-2),2]>b[(i-3),2]) & (b[(i+2),2]>b[(i+3),2]))
          {
          result[i]=b[i,1]
          }
      }
      }

      else
      {
      b<-te[,c(1,2)]

      ## detection of peaks in the chromatogram

      for (i in 4:(dim(b)[1]-4))
      {
            if (b[i,2]>5000 & (b[i,2]>b[(i-1),2]) & (b[i,2]>b[(i+1),2]) & (b[(i-1),2]>b[(i-2),2]) & (b[(i+1),2]>b[(i+2),2]) & (b[(i-2),2]>b[(i-3),2]) & (b[(i+2),2]>b[(i+3),2]))
            {
            result[i]=b[i,1]
            }
      }
      }


result<-result[is.na(result)==FALSE]

## treat the first peak

te_tempo<-match(te[,1], result[1])
names(te_tempo)<-1:length(te_tempo)
tb<-te_tempo[is.na(te_tempo)==FALSE]
li<-names(tb)

   if (apex==FALSE)
   {
   liInf<-min(as.numeric(li))-1
   liSup<-max(as.numeric(li))+1
   li<-c(liInf, li, liSup)
   }

#tempo<-mean(te[li,]) deprecated
tempo<-colMeans(te[li,])

## Do similar analyses for the other peaks

for (j in 2:length(result))
{
te_tempo<-match(te[,1], result[j])
names(te_tempo)<-1:length(te_tempo)
tb<-te_tempo[is.na(te_tempo)==FALSE]
li<-names(tb)

    if (apex==FALSE)
    {
    liInf<-min(as.numeric(li))-1
    liSup<-max(as.numeric(li))+1
    li<-c(liInf, li, liSup)
    }
#change mean colMeans
tempo<-rbind(tempo, colMeans(te[li,]))
}

caj<-rep(L_an[a], dim(tempo)[1])
tempo<-cbind(caj, tempo)
res_data_temp[[a]]<-tempo
print(a)

save(res_data_temp, file=paste(Mypath, "/", "save_list_temp.rda", sep=""))
}

res_data<-res_data_temp[[1]]

    if (length(res_data_temp)>1)
       for (l in 2:length(res_data_temp))
       {
       res_data<-rbind(res_data, res_data_temp[[l]])
       }

res_data<-as.data.frame(res_data)

    for (i in 4:ncol(res_data))
    {
    res_data[,i]<-as.numeric(as.vector(res_data[,i]))
    data_fin<-res_data[,-3]
    }

colnames(data_fin)<-c("analysis", "retention_time", mzall)
rownames(data_fin)<-1:nrow(data_fin)
#just in case mz values entered by user are different from mzall
data_fin<-cbind(data_fin[,1:2],data_fin[,which(colnames(data_fin)==mz[1]):which(colnames(data_fin)==mz[length(mz)])])
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
#Rprof(NULL)
#summaryRprof(filename = "Rprof.out")$sampling.time
print(paste("A data file has been generated in the folder:", Mypath, cat("\n")))
return(data_fin)
}



}#end if agilent or ascii
if (DataType=="CDF")
{
	st<-strsplit(date(), " ")[[1]]
	stBis<-strsplit(st[4], ":")[[1]]
	Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
	Date<-paste(st[1], st[2], st[3], Hour, sep="_")
	Mypath<-paste("output_MSDataCreation", "_", "result", Date, sep="")
	dir.create(Mypath)

	filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
	filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),
	collapse = "|") 

	require("xcms")
	require("tcltk")
	#specific variables for CDF error handling
	mz_min <- vector()
	mz_max<-vector()
	pbmzmin<-0
	pbmzmax<-0


	if(is.null(pathCDF)==TRUE||missing(pathCDF)==TRUE){
		cdffiles<-list.files(path=tclvalue(tkchooseDirectory(title="Select a directory with CDF (mzXML) files AND peaklist files")), pattern=filepattern, recursive = TRUE, full.names = TRUE) 
		if (length(cdffiles)==0){
			cat("Stop No CDF or mzXML files detected \n ")
		}
		
		peakfiles<-list.files(path=dirname(cdffiles), pattern="peaklist.txt", recursive = TRUE, full.names = TRUE)
		#if Agilent rteres.txt files are used instead of peaklist.txt
		if (length(peakfiles)==0  && length(cdffiles)>0){
			peakfiles2<-list.files(path=dirname(cdffiles), pattern="rteres.txt", recursive = TRUE, full.names = TRUE)
			if(length(peakfiles2)>0){
				cat("No peaklist.txt files found, \n There are rteres.txt files from Agilent available in your CDF directory, \n Please start MS.DataCreation with DataType=Agilent \n")
			}
			else
			{
				cat("WARNING No peaklist or rteres files founds ! \n \n")			
			}
		}
	}
	else{
		cdffiles<-list.files(path=pathCDF, pattern=filepattern, recursive = TRUE, full.names = TRUE) 
		if (length(cdffiles)==0){
			cat("Stop No CDF or mzXML files detected \n ")
		}
		peakfiles<-list.files(path=dirname(cdffiles), pattern="peaklist.txt", recursive = TRUE, full.names = TRUE)
		#if Agilent rteres.txt files are used instead of peaklist.txt
		if (length(peakfiles)==0  && length(cdffiles)>0){
			peakfiles2<-list.files(path=dirname(cdffiles), pattern="rteres.txt", recursive = TRUE, full.names = TRUE)
			if(length(peakfiles2)>0){
				cat("No peaklist.txt files found, \n There are rteres.txt files from Agilent available in your CDF directory, \n Please start MS.DataCreationCDF function \n")
				
			}
			else
			{
			cat("\n WARNING No peaklist or rteres files founds ! \n")
			
			}
		}
	}
if(length(cdffiles)>0){
	#check if length cdffiles = peakfiles
	cdf_error<-function(){
	
		for (m in 1:length(cdffiles)) {
				obj<-xcms::xcmsRaw(cdffiles[m]) 
				mz_min[m]<-min(obj@mzrange)
				mz_max[m]<-max(obj@mzrange)
			}
		
		
			if ((mz=="all")==TRUE||is.null(mz)==TRUE){
					mz<-max(mz_min):min(mz_max)
			}
			if (max(mz_min)>mz[1]){
						mz<-max(mz_min):mz[length(mz)]
						pbmzmin<-1
			}	
			if (min(mz_max)<mz[length(mz)]){
						mz<-mz[1]:min(mz_max)
						pbmzmax<-1
					}	
			

			if (pbmzmin!=0){
						
						print(paste("mz minimum value was set to ", mz[1], cat("\n")))
						pbmzmin<-0
					}
			if (pbmzmax!=0){
						print(paste("mz maximum value was set to ", min(mz_max), cat("\n")))
						pbmzmax<-0
				}
		
		assign("mz", mz, envir=as.environment("e1"))
		return(mz)
	}




	if (length(cdffiles)==length(peakfiles)){

		#search for differences in user defined mz vector and reel mz present in CDF files	
		if(missing(mz)==TRUE){
			mz="all"
		}
		mz<-cdf_error()
		print(paste("mz:",mz[1],":",mz[length(mz)],cat("\n")))
		
	#boucle sur les fichiers 
	peaks<-NULL #initialise matrice peaks
	for (fi in 1:length(cdffiles)){
		peaklist<-read.table(file.path(dirname(cdffiles[fi]),"peaklist.txt"), header=TRUE,sep="\t")
		obj<-xcms::xcmsRaw(cdffiles[fi]) #crée un fichier xcms
		scanmin<-matrix()
		scanmax<-matrix()
		scantop<-matrix()
		RT<-matrix()
		quant=FALSE
		
		peakslg=dim(peaklist)[1]
		if (length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))>0){
			quant=TRUE
		}
			
			for (j in 1:peakslg){
			
				#changer ici si un rtmin rtmax est donné et utiliser scantime
				#calculer RT avec scantime vs scanindex
				if(apex==FALSE){
					scanmin[j]=peaklist[j,"firstscan"]
					scanmax[j]=peaklist[j,"lastscan"]
					scantop[j]=peaklist[j,"maxscan"]
					Rang=scanmax[j]-scanmin[j]
					sc5=Rang*5/100
					MS=round(scantop[j]-sc5, digits=0):round(scantop[j]+sc5, digits=0)
				}
				if (apex==TRUE){
					scantop[j]=peaklist[j,"maxscan"]
					MS=scantop[j]
				}
				RT[j]<-peaklist[j,grep("[Rr][Tt]",colnames(peaklist))[1]]
				#tests quant==TRUE
			
				if(quant==TRUE){
				#Max column for area is 2
					quantif<-matrix(nrow=1,  ncol=2)
			
				if (length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))>2){
						
						quantif[1,1]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]], sep="\t")
						quantif[1,2]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[2]], sep="\t")
				}
				else
				{
					    quantif[1,1]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]], sep="\t")
						quantif[1,2]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[2]], sep="\t")
				}
				if (length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))<2){
						
						quantif[1,1]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]], sep="\t")
						#duplication of the area column for further analysis
						quantif[1,2]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]], sep="\t")
				}	
				else
				{
					    quantif[1,1]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]], sep="\t")
						quantif[1,2]<-paste(peaklist[j,grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[2]], sep="\t")
				}
				
				}
				
				
				
				#if peaklist does not contains RT column
				if (is.null(RT[j])==TRUE){
					RT[j]<-obj@scantime[scantop[j]]/60
				}
				
				an<-matrix(ncol=length(mz), nrow=length(MS))
			for (i in 1:length(MS)) {
				for (p in 1:length(mz)){
					an[i, p] <- 0 #set 0 intensity for every m/z
				}
			}
			colnames(an)<-mz
			if(apex==FALSE){
				scan=scanmin[j]
			}
			else
			{
				scan=scantop[j]
			}
			for (i in 1:length(MS)) {
				Ma<-match(round(obj@env$mz[(obj@scanindex[scan]+1):obj@scanindex[scan+1]],0),colnames(an))
				for (p in 1:length(Ma)){
					if (is.na(Ma[p]) == FALSE) {
									 an[i, Ma[p]] <- obj@env$intensity[(obj@scanindex[scan]+1):obj@scanindex[scan+1]][p] #put the right intensity to the m/z
							}
						
				}
			scan=scan+1
			}
			#faire moyenne des intensité/mz puis normalisation sur le max
			#normalisation of m/z intensity 
			#first average intensity for each m/z and search of the m/z with the highest intensity
			max=0
			av<-matrix()
			for (p in 1:length(mz)){
				av[p] <- colMeans(an)[p]
				if (av[p]>max){
					max=av[p]
					mzmax=p
				}
			}

			norm<-matrix(ncol=length(mz))

			peakstemp<-matrix(ncol=length(mz),nrow=1)
			colnames(peakstemp)<-mz
			for (p in 1:length(mz)){
				norm[p]<-(100*av[p])/av[mzmax] #normalisation
				peakstemp[1,p]<-norm[p] #sauver la ligne d'info du pic
			}
			
			#en fonction des colnames de peaklist.txt ajouter RT et si quant=T aire airemax %aire RTmin RTsec
			#ajouter les bons noms de colonnes
			
			if (quant==FALSE){
			peakstemp<-cbind(basename(cdffiles[fi]),RT[j],peakstemp)
			}
			if(quant==TRUE){
		
					peakstemp<-cbind(basename(cdffiles[fi]),RT[j],quantif, peakstemp)
			}
			av<-NULL
			norm<-NULL
			an<-NULL
			gc()
			peaks<-rbind(peaks,peakstemp)
			}#end for peaklg
			
		}#end for cdffiles


	colnames(peaks)[1]<-"analysis" 
	colnames(peaks)[2]<-"retention_time" 
	#ajouter AREA si besoin
	if(quant==TRUE){
	#MS.clust need two area column for the moment
	if (length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))<2){
		print(paste(length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))," One Quantification columns have been found", sep=" "))
		print(colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))])
		colnames(peaks)[3]<-colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))]
		colnames(peaks)[4]<-colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))]
	}
	else
	{
			if (length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))>2){
						print(paste(length(grep("[Qq][Uu][Aa][Nn]",colnames(peaklist)))," Quantification columns have been found, Max is 2 /n Only the 2 first will be used", sep=" "))
						print(paste(colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]],colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[2]], sep=" "))
			}
		colnames(peaks)[3]<-colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[1]]
		colnames(peaks)[4]<-colnames(peaklist)[grep("[Qq][Uu][Aa][Nn]",colnames(peaklist))[2]]
	}
	}
	write.table(peaks, file=file.path(Mypath,"initial_DATA.txt"),quote=c(1,2), row.names=FALSE)
	print(paste("A data file has been generated in the folder:", Mypath, cat("\n")))
	return(data.frame(peaks, check.names=FALSE))
	}#end if meme longueur

	if (length(peakfiles)!=0 & length(cdffiles)!=length(peakfiles)){
	#if there are different number of CDF files vs peaklist.txt files
		if (length(cdffiles)>length(peakfiles)){
			cat("Some peaklist.txt files are missing \n")
		}
		if (length(cdffiles)<length(peakfiles)){
			cat("Some CDF or mzXML files are missing \n")
		}
	}
}#end length sup 0
}#enf if CDF
}#end function

