#' SeqGrapheR - graph based interactive visualization of DNA sequence clusters
#'
#' Function start GUI for interactive graph based visualization of DNA sequence clusters. 
#' The GUI is written using the gWidgets and use RGtk2 graphical toolkit. This function
#' does not take any parameters. Note that SeqGrapheR require external programs
#' blastall and dotter for some functions.  
#' 
#' @export SeqGrapheR
#' @importFrom utils read.table write.table
#' @importFrom grDevices dev.off png rgb
#' @importFrom graphics plot.new points segments plot abline hist par
#' @import Biostrings
#' @importFrom gWidgets svalue gconfirm gfile ginput delete gmessage galert dispose ggroup ggraphics gslider glabel gwindow gtable gbutton size addhandlerdestroy gmenu gframe gimage addSpace glayout gtext gstatusbar visible  
#' @importFrom igraph fastgreedy.community "V<-" membership vcount graph.data.frame simplify V E get.edgelist neighborhood  degree ecount layout.fruchterman.reingold induced.subgraph    
#' @importFrom rggobi edges "edges<-" ids ggobi "glyph_colour<-" "glyph_type<-" "glyph_size<-" glyph_colour ggobi_count ggobi_get glyph_type glyph_size displays pmode ggobi_display_get_tour_projection colorscheme "$<-.GGobi" "[[.GGobi" "ids<-" "[<-.GGobi" "[.GGobi" "[[<-.GGobi" "$<-.GGobi" "[.GGobiData" "[[<-.GGobiData" "[[.GGobiData" "$<-.GGobiData" "$.GGobiData" "variables<-" "shadowed<-" 
#' @import gWidgetsRGtk2
#' @import cairoDevice
## @importMethodsFrom igraph
## @importMethodsFrom rggobi
#' @examples 
#' ## do nut run!
#' ## SeqGrapheR()


SeqGrapheR=function(){    # main function
	options(warn=-1)
	options("guiToolkit"="RGtk2")
	#check for Biostring version
	if (exists("readDNAStringSet")){
			read.DNAStringSet=readDNAStringSet
			read.BStringSet=readBStringSet
			write.XStringSet=writeXStringSet
			
		}
	
	
	########################## INTERNAL  FUNCTIONS ################################################################
	########################install.packages("SeqGrapheR_0.4.8.3.tar.gz",repos=NULL)##############################################################################
	
	
	formatDNA=function(dna,width=60){
		N=nchar(dna)
		s=strsplit(dna,split='')[[1]]
		sout=paste((lapply((split(s,f=floor(seq(0,N-1)/60))),paste,collapse="")),collapse="\n")
	}
	
	writeFasta=function(seqs,file,format=F){
		#this function is to replace write.XStringSet which is too slow (bug in Biostrings??)
		seqnames=names(seqs)
		seqchars=as.character(seqs)
		if (format){
			seqcharsF=sapply(seqchars,formatDNA)
			cat(paste(">",seqnames,"\n",seqcharsF,"\n",sep=''),sep='',file=file)
		}else{
			cat(paste(">",seqnames,"\n",seqchars,"\n",sep=''),sep='',file=file)	
		}
		
	}
	
	fastgreedy=function(GG){
		eb=fastgreedy.community(GG)
		cutoff=which(eb$modularity==max(eb$modularity))-1  # select cutoff
		mb=membership(GG)
		nc=length(unique(mb))   # number of clusters
		N=vcount(GG) 
		new.mb=mb
		new.mb=as.numeric(factor(new.mb,levels=sort(unique(new.mb))[order(rle(sort(new.mb))$length,decreasing=T)] ))-1
		mb.counts=rle(sort(new.mb))$length
		mb.names=rle(sort(new.mb))$values
		CL.names=paste("cluster",mb.names+1,sep="")
		FG=list(mb=new.mb,counts=mb.counts,names=mb.names)
		FG
	}
	
	
	blastall=function(databaseSeq,querySeq,programname='blastn',m=8,params="",protein=FALSE,databaseParams='')
	{
		
		# wrapper for blastall program
		# blastall -p programname -d databasefilename -i queryfilename -o outputfilename
		# m=8 specify tabular output - then type of output is data.frame
		#create temporary files with sequences
		database=tempfile()
		query=tempfile()
		output=tempfile()
		
		#check names
		writeFasta(databaseSeq,file=database)
		writeFasta(querySeq,file=query)
		
		#write.XStringSet(databaseSeq,file=database,format='fasta')   # too slow!!!
		#write.XStringSet(querySeq,file=query,format='fasta')
		
		#create database:
		cmd=paste("formatdb -i",database,"-p",protein,databaseParams)
		system(cmd)
		cmd=paste("blastall -p",programname,"-d",database,"-i",query,"-m",m,params,'-o',output)
		system(cmd)
		if (m==8){
			blastOut=read.table(output,sep="\t",header=FALSE,as.is=TRUE)
			colnames(blastOut)=c('Query','Subject','percIdentity','alignmentLength','mistmatches','gapOpenings','q.start','q.end','s.start','s.end','Evalue','bitScore')
		}else{
			blastOut=scan(output,what=character())
		}
		unlink(output)
		unlink(database)
		unlink(query)
		blastOut
	}
	
	megablast=function(databaseSeq,querySeq,D=3,params="",databaseParams='')
	{	# wrapper for mgblast program
		#megablast -d databasefilename -i queryfilename -D 2 -o outputfile
		# D=3 specify tabular output
		
		if (params=='tgicl'){
			params=' -s 25 -W18 -UF -X40 -JF -F "m D" -v90000000 -b90000000 -H 3000'
		}
		#create temporary files with sequences
		database=tempfile()
		query=tempfile()
		output=tempfile()
		write.XStringSet(databaseSeq,filepath=database,format='fasta')
		write.XStringSet(querySeq,filepath=query,format='fasta')
		#create database:
		cmd=paste("formatdb -i",database,"-p F",databaseParams)
		system(cmd)
		cmd=paste("megablast","-d",database,"-i",query,"-D",D,params,'-o',output)
		system(cmd)
		if (D==3){
			blastOut=read.table(output,sep="\t",header=FALSE,as.is=TRUE)
			colnames(blastOut)=c('Query','Subject','percIdentity','alignmentLength','mistmatches','gapOpenings','q.start','q.end','s.start','s.end','Evalue','bitScore')
		}else{
			blastOut=scan(output,what=character())
		}
		unlink(output)
		unlink(database)
		unlink(query)
		blastOut
	}
	
	
	dotter=function(seq1,seq2=NULL){
		#require external program dotter
		
		if (is.null(seq2)) {
			seq2=seq1	
		}
		
		if (class(seq1)!="DNAStringSet"){
			seq1=DNAStringSet(seq1)
		}
		
		if (class(seq2)!="DNAStringSet"){
			seq2=DNAStringSet(seq2)
		}
		
		
		sf1=tempfile('seq1')
		write.XStringSet(seq1,filepath=sf1)
		sf2=tempfile('seq2')
		write.XStringSet(seq2,filepath=sf2)
		system(paste('dotter',sf1,sf2),wait=FALSE)
		Sys.sleep(2)
		unlink(c(sf1,sf2))
	}
	
	
	blast2graph=function(blastTable,options=list(thr=200,LCOV=55,SCOV=NULL,PID=90,OVL=NULL),method="thr",seqs=NULL)
	{
		# create graph from blast table
		# blast table can be fitred for quality of hit using thr param - minimal bitScore when method is thr
		# or  LCOV,SCOV,OVL,PID parameters when method is cov - then sequences must be provided to get seq length
		N=dim(blastTable)[1]
		ind=1:N
		i=1
		# blastTable .. output of blast -tabular formatdb
		# thr - blast bitScore threshold
		# seqs - sequences - used when thrseq as threshold is used
		# thrseq - threshold as a function of the sequence length
		if (method=='thr'){
			vvw=blastTable[,c(1,2,12)]
			vvw=vvw[vvw[,3]>options$thr,]
			colnames(vvw)=c(1,2,'weight')
			g=graph.data.frame(vvw,directed=FALSE)
			return(g)
		}
		if (method=='cov' & !is.null(seqs)){
			print(i);i=i+1
			
			seqLength=nchar(seqs)
			names(seqLength)=gsub(" .+","",names(seqs))
			print(i);i=i+1
			Qlength=seqLength[blastTable$Query]
			Slength=seqLength[blastTable$Subject]
			print(i);i=i+1
			Qovl=abs(blastTable$q.end-blastTable$q.start)+1
			Sovl=abs(blastTable$s.end-blastTable$s.start)+1
			print(i);i=i+1
			Lind=(Qlength<Slength)+1
			LCOV=(cbind(Qovl,Sovl))[cbind(ind,Lind)]/(cbind(Qlength,Slength))[cbind(ind,Lind)]
			SCOV=cbind(Sovl,Qovl)[cbind(ind,Lind)]/cbind(Slength,Qlength)[cbind(ind,Lind)]
			OVL=base::pmax(Qovl,Sovl)
			PID=blastTable$percIdentity
			
			if (is.null(options$LCOV)) {
				condLcov=TRUE
			}else{
				condLcov=(100*LCOV)>=options$LCOV
			}
			
			if (is.null(options$SCOV)) {
				condScov=TRUE
			}else{
				condScov=(100*SCOV)>=options$SCOV
			}
			if (is.null(options$PID)) {
				condPID=TRUE
			}else{
				condPID=PID>=options$PID
			}
			
			if (is.null(options$OVL)) {
				condOVL=TRUE
			}else{
				condOVL=OVL>=options$OVL
			}
			vvw=blastTable[,c(1,2,12)]
			vvw=vvw[condLcov & condScov & condPID & condOVL,]
			colnames(vvw)=c(1,2,'weight')
			g=simplify(graph.data.frame(vvw,directed=FALSE))
			return(g)
		}	
	}
	
	getACEcontigs=function(aceFile){
		ACE=scan(aceFile,sep="\n",what=character())
		contigNames=grep("^CO ",ACE,value=TRUE)
		pos1=grep("^CO ",ACE)
		pos2=c(pos1[-1],length(ACE))
		contigNames=gsub("^CO ",">",contigNames)
		contigNames=gsub(" .*","",contigNames)
		N=length(contigNames)
		contigReads=list()
		for (i in seq_along(pos1)){
			readNames=grep("^AF ",ACE[pos1[i]:pos2[i]],value=TRUE)
			readNames=gsub("^AF ","",readNames)
			readNames=gsub(" +.+$","",readNames)
			contigReads[[i]]=readNames
		}
		#sort by number of reads
		ord=order(sapply(contigReads,length),decreasing=TRUE)
		contigNames=contigNames[ord]
		contigReads=contigReads[ord]
		contigNames=gsub(" .+$","",contigNames)
		contigNames=gsub(">","",contigNames)
		names(contigReads)=contigNames
		# print contigs summary
		contigReads
	}
	
	
	##### GUI handlers functions ################################################################
	
	setdir=function(File){
		if (!is.na(File)){
			newDir=dirname(File)
			setwd(newDir)
		}
	}
	
	
	IdsEditorCleanUp=function(){ 
# read content of ids editor window, valid separators are ,;space and tab
# remove all unknown ids - not present in graph 
# format Ids in editor window, remove duplicates ids
		if (exists("GL",envir=envGL)){
			# check content of the gui Ids editor window:
			windContent=svalue(gui$IdsEditor)
			#vectorize:
			windContent=unlist(strsplit(windContent,split="[ ,;\t\n]"))
			#check names
			svalue(gui$IdsEditor)=paste(unique(windContent[windContent %in% V(envGL$GL$G)$name]),collapse=", ")
		}
	}
	
	getIdsFromEditor=function(){ 
# read content of ids editor window, valid separators are ,;space and tab
# remove all unknown ids - not present in graph 
		if (exists("GL",envir=envGL)){
			# check content of the gui Ids editor window:
			windContent=svalue(gui$IdsEditor)
			#vectorize:
			windContent=unlist(strsplit(windContent,split="[ ,;\t]"))
			windContent=unique(windContent[windContent %in% V(envGL$GL$G)$name])
			#check names
		}
	}
	
	addIdsToEditor=function(idsNew){   
		idsOriginal=getIdsFromEditor()
		ids=unique(c(idsOriginal,idsNew))
		svalue(gui$IdsEditor)=paste(ids,collapse=", ")
		ids
	}
	
	
	
	removeIdsFromEditor=function(ids){   
		idsOriginal=getIdsFromEditor()
		idsNew=idsOriginal[!(idsOriginal %in% ids)]
		svalue(gui$IdsEditor)=paste(idsNew,collapse=", ")
		ids
	}
	
	
	
# show graph in ggobi 
        # added to show pairs edges in graph:
        get_pairs=function(GL){
                                        # check if last character label pair
            label = substring(V(GL$G)$name,nchar(V(GL$G)$name))
            if (length(unique(label))==2){
                id=gsub(".{1}$","",V(GL$G)$name)
                tbl = table(id, label)
                tbl=tbl[tbl[,1]==1 & tbl[,1]==1,][]
                v1 = paste(rownames(tbl),label[1],sep="")
                v2 = paste(rownames(tbl),label[2],sep="")
                return(cbind(v1,v2))
            }else{
                return(NULL)
            }
            
        }



        showGraph=function(graphLayout,fgcolor=NULL){
            gf=graphLayout$L
            rownames(gf)=V(graphLayout$G)$name
            ggobiObject=ggobi(gf)

            
                                        #display edges:
            if (gconfirm(message='Include edges in the graph?\n\n(can be computationaly / memory demanding)',icon='warning')){
                weight=E(graphLayout$G)$weight
                ggobiObject$whg=data.frame(weight)
                e=get.edgelist(graphLayout$G);rownames=c('src','dest')
                edges(ggobiObject$whg)=e
            }

                                        # pair edges, make it optional - it slows import:
            ep = get_pairs(graphLayout)
            if (!is.null(ep)){
                if (gconfirm(message='Include edges for read pairs?)',icon='warning')){
                    ggobiObject$pairs = data.frame(weight = rep(1,nrow(ep)))
                    edges(ggobiObject$pairs) = ep
                }
            }
            
            if (!is.null(fgcolor)){
                glyph_colour(ggobiObject[1])=fgcolor
            }
                                        #gd=displays(ggobiObject)[[1]]
                                        #pmode(gd) = 'Rotation'
            ggobiObject
	}
	
	defHandlerAdd2Editor=function(h,...){
		# get selected Ids
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)>0){ #check is anything is selected
			selectedIds=unique(unlist(envGL$IdsListIds[selectedIdsIndex]))
			addIdsToEditor(selectedIds)
		}
	}
	
	
	
	
	formatIdsList=function(IdsList,IdsListDf){  # lists and data.frame description ...
		idsFormated=character()
		for (i in 1:length(IdsList)){
			
			
			idsFormated=c(idsFormated,paste(">",IdsListDf[i,1]," ",IdsListDf[i,2]," ",IdsListDf[i,3],sep=''))
			idsFormated=c(idsFormated,paste(paste(IdsList[[i]],sep=" ", collapse=" ")))
		}	
		
		idsFormated
	}
	
	defHandlerSaveMultipleLists=function(h,...){
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)>0){ #check is anything is selected
			selectedLists=envGL$IdsListIds[selectedIdsIndex]
			selectedListsDf=envGL$IdsListDataFrame[selectedIdsIndex,]
			selectedLists=formatIdsList(selectedLists,selectedListsDf)
			lsFile=gfile("Save as ",type='save')
			if (!is.na(lsFile)){
				cat(selectedLists,file=lsFile,sep="\n")
			}
		}
	}
	
	defHandlerLoadMultipleLists=function(h,...){
		lsFile=gfile("Select file",type='open')
		if (!is.na(lsFile)){
			setdir(lsFile)
			importedList=scan(file=lsFile,sep='\n',what=character())
			N=length(importedList)
			df=importedList[seq(1,N,2)]
			df=gsub(">","",df,fixed=TRUE)
			df=gsub("\t"," ",df,fixed=TRUE) # alternative separatorsaddhandlerdestroy
			df=strsplit(df,split=" ")
			Group=sapply(df,FUN=function(x)x[1])
			Subset=sapply(df,FUN=function(x)x[2])
			Number_of_reads=sapply(df,FUN=function(x)as.integer(tail(x,1)))
			listsDf=data.frame(Group=Group,Subset=Subset,Number_of_reads=Number_of_reads,stringsAsFactors=FALSE)
			ids=importedList[seq(1,N,2)+1]
			ids=gsub("\t"," ",ids,fixed=TRUE) # alternative separators
			ids=strsplit(ids,split=' ')
			addList(listsDf,ids)
			
		}
		
	}
	
	defHandlerLoadIdsTable=function(h,...){
		lsFile=gfile("Select file",type='open')
		if (!is.na(lsFile)){
			setdir(lsFile)
			idsTable=read.table(lsFile,sep='',header=FALSE,as.is=FALSE)[,1:2]
			ids=list(c(t(idsTable)))
			Group=ginput(message="Enter Group name",icon="question")
			Subset=ginput(message="Enter Subset name",icon="question")
			listDf=data.frame(Group=Group,Subset=Subset,Number_of_reads=nrow(idsTable),stringsAsFactors=FALSE)
			print(listDf)
			
			
			addList(listDf,ids)
			
		}
	}
	
	defHandlerRemoveFromEditor=function(h,...){
		# get selected Ids
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)>0){ #check is anything is selected
			selectedIds=unique(unlist(envGL$IdsListIds[selectedIdsIndex]))
			removeIdsFromEditor(selectedIds)
		}
	}
	
	
	defHandlerClearEditor=function(h,...){
		svalue(gui$IdsEditor)=""
	}
	
	defHandlerDeleteSelectedList=function(h,...){
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		envGL$IdsListVisible[selectedIdsIndex]=FALSE
		updateIdsListDataWindow()
		
	}
	
	
	defHandlerYellow2editor=function(h,...){
		#get selected:
		if (DataLoaded()) {
			selectedIds=V(envGL$GL$G)$name[glyph_colour(envGL$g[1])==9]
			#click
			svalue(gui$move2editorArrow)=paste(iconDir,"/yellowArrow2.png",sep="")
			Sys.sleep(0.1)
			svalue(gui$move2editorArrow)=paste(iconDir,"/yellowArrow.png",sep='')
			addIdsToEditor(selectedIds)
		}
	}	
	
	defHandlerColors2lists=function(h,...){
		colors=glyph_colour(envGL$g[1])
		colorLists=lapply(1:9,FUN=function(x)V(envGL$GL$G)$name[colors==x])
		csName=ginput("enter Color scheme name")
		colorListsDf=data.frame(Group=rep(csName,9),Subset=paste('col',1:9,sep="_"),Number_of_reads=as.integer(sapply(colorLists,length)),stringsAsFactors=FALSE)
		addList(colorListsDf,colorLists)
	}	
	
	defHandlerMakeNewList=function(h,...){
		ids=getIdsFromEditor()
		Number_of_reads=length(ids)
		if (Number_of_reads>0){
			Group=ginput(message="Enter Group name",icon="question")
			Subset=ginput(message="Enter Subset name",icon="question")
			# put list into gtable
			listsDf=data.frame(Group=Group,Subset=Subset,Number_of_reads=Number_of_reads,stringsAsFactors=FALSE)
			addList(listsDf,list(ids))
		}
		
	}
	
	defHandlerShowClusters=function(h,...){
		if (DataLoaded()) {
			fg=fastgreedy(envGL$GL$G)
			ids=split(V(envGL$GL$G)$name,f=fg$mb)
			listDf=data.frame(Group=fg$names,Subset="fastgreedy",Number_of_reads=fg$counts)
			addList(listDf,ids)
		}
		
	}
	
	defHandlerShowNeighbors=function(h,...){
		#get selected:
		if (DataLoaded()) {
			selids=which(glyph_colour(envGL$g[1])==9)
			selids=unique(unlist(neighborhood(envGL$GL$G,1,selids)))
			glyph_colour(envGL$g[1])[selids]=9
			glyph_colour(envGL$g[1])[!selids]=1
			
		}
	}
	# open GL file and crerate ggobi object,
	defHandlerOpenGL=function(h,...){
		glfile=gfile(text='select graph file', type = "open", filter = list("GL graph" = list(patterns = c("*.GL")), "All files" = list(patterns = c("*"))), handler = function(h,...) return(h$file))
		if (!is.na(glfile)){    # if not canceled)
			setdir(glfile)
			GL=NULL
			try(load(glfile),silent=TRUE) 
			if (!is.null('GL')) { #in loaded object
				# check if g exist in envGL and close and remove everything
				if (exists('g',envir=envGL)){
					while(ggobi_count()>0) close(ggobi_get())  # will close all ggobi instances
					rm_gW('degreeWindow',envir=gui)   # will remove degrreWindow if open
					rm(list=ls(envir=envGL),envir=envGL)
				}
				envGL$GL=GL  
				envGL$g=showGraph(envGL$GL)
				createEmptyIdsLists()
				#initial first row:
				envGL$IdsListDataFrame[1,1:2]=c("ALL","ALL")
				envGL$IdsListDataFrame[1,3]=as.integer(vcount(envGL$GL$G))
				envGL$IdsListDataFrame[1,4]=as.integer(1)
				envGL$IdsListVisible=TRUE
				envGL$IdsListIds[[1]]=V(envGL$GL$G)$name
				envGL$IdsListFreq[[1]]=NA
				updateIdsListDataWindow()
				
				if (length(names(envGL$g))==2){# do not show histograpm when no edges are not imported to ggobi
					degreeHistogram()
					#add degree info to ggobi
					envGL$g$gf$degree=degree(envGL$GL$G)
				}else{
					glyph_type(envGL$g[1])=1
					if (exists("currentGraph",envir=gui)){ #remove it
						delete(gui$graphicsWindow,gui$currentGraph)
						rm('currentGraph',envir=gui)
					}
					
				}
				
			}else{#wrong file
				gmessage("Wrong file format", title="message",icon ="error")
			}
		}
		svalue(gui$StatusBar)=paste(glfile," :  V=",vcount(envGL$GL$G)," E=",ecount(envGL$GL$G),sep='')
	}
	
	## open ncol file, calculate layout, and create ggobi object
	defHandlerOpenNcol=function(h,...){
		nncolFile=gfile(text='select graph file in ncol format', type = "open", filter = list("ncol graph file" = list(patterns = c("*.ncol",'*.NCOL')), "All files" = list(patterns = c("*"))), handler = function(h,...)return(h$file))
		if (!is.na(nncolFile)){
			setdir(nncolFile)
			Gncol=read.table(nncolFile,header=FALSE,sep="",colClasses=c("character","character","numeric"))
			colnames(Gncol)=c(1,2,"weight")
			GL=list()
			GL$G=graph.data.frame(Gncol,directed=FALSE)
			GL$G=simplify(GL$G)
			gal=galert("Graph layout is being calculated...",delay=100000)
			Sys.sleep(1)  # have time to display message
			GL$L=layout.fruchterman.reingold(GL$G,dim=3,verbose=TRUE)
			
			removeAll()
			envGL$GL=GL
			envGL$g=showGraph(envGL$GL)
			createEmptyIdsLists()
			#initial first row:
			envGL$IdsListDataFrame[1,1:2]=c("ALL","ALL")
			envGL$IdsListDataFrame[1,3]=as.integer(vcount(envGL$GL$G))
			envGL$IdsListDataFrame[1,4]=as.integer(1)
			envGL$IdsListVisible=TRUE
			envGL$IdsListIds[[1]]=V(envGL$GL$G)$name
			envGL$IdsListFreq[[1]]=NA
			updateIdsListDataWindow()
			degreeHistogram()
			#add degree info to ggobi
			envGL$g$gf$degree=degree(envGL$GL$G)
			dispose(gal)
		}
	}
	
	
# open sequences, run megablast, create graph, calculate layout, show graph. 
	defHandlerOpenSequences2Graph=function(h,...){
		fastaFile=gfile(text='select fasta file for graph construction', type = "open", filter = list("fasta" = list(patterns = c("*.fas",'*.fna',"fasta")), "All files" = list(patterns = c("*"))), handler = function(h,...)return(h$file))
		if (!is.na(fastaFile)){    # if the action was not canceled
			setdir(fastaFile)
			GL=list()
			seqs=read.DNAStringSet(fastaFile,format="fasta")
			names(seqs)=gsub(" .+","",names(seqs))  # make short names
			print('running blast')
			Sys.sleep(1)  # have time to display message
			megablastTable=megablast(seqs,seqs,params='tgicl')
			GL$G=blast2graph(megablastTable)  # make graph
			GL$L=layout.fruchterman.reingold(GL$G,dim=3,verbose=TRUE)
			# show 
			removeAll()
			
			# make sequence index
			seqIndex=1:length(seqs)
			names(seqIndex)=names(seqs)
			GL$Seq=seqs[seqIndex[V(GL$G)$name]]    # biostring object as part of GL
			V(GL$G)$length=unname(nchar(GL$Seq))   # length as vertex attribute
			
			# put everything into envGL eviroment
			envGL$seqIndex
			envGL$magablastTable=megablastTable
			envGL$GL=GL
			envGL$g=showGraph(envGL$GL)
			createEmptyIdsLists()
			#initial first row:
			envGL$IdsListDataFrame[1,1:2]=c("ALL","ALL")
			envGL$IdsListDataFrame[1,3]=as.integer(vcount(envGL$GL$G))
			envGL$IdsListDataFrame[1,4]=as.integer(1)
			envGL$IdsListVisible=TRUE
			envGL$IdsListIds[[1]]=V(envGL$GL$G)$name
			envGL$IdsListFreq[[1]]=NA
			updateIdsListDataWindow()
			degreeHistogram()
			#add degree info to ggobi
			envGL$g$gf$degree=degree(envGL$GL$G)
			#add Length info to ggobi
			envGL$g$gf$length=V(envGL$GL$G)$length
		}
	}
	
	defHandlerImportSequences=function(h,...){
		#assume that graph is already loades from GL or ncol formatdb
		if (DataLoaded()){
			fastaFile=gfile(text='select fasta file', type = "open", filter = list("fasta" = list(patterns = c("*.fas",'*.fna',"fasta")), "All files" = list(patterns = c("*"))), handler = function(h,...)return(h$file))
			if (!is.na(fastaFile)){    # if the action was not canceled
				setdir(fastaFile)
				seqs=read.DNAStringSet(fastaFile)
				names(seqs)=gsub(" .+","",names(seqs))  # short names only
				# check if contain same names as graph GL$G
				if (sum(!(V(envGL$GL$G)$name %in% names(seqs)))==0){
					envGL$seqIndex=1:length(seqs)
					names(envGL$seqIndex)=names(seqs)
					envGL$GL$Seq=seqs[envGL$seqIndex[V(envGL$GL$G)$name]]    # biostring object as part of GL
					V(envGL$GL$G)$length=unname(nchar(envGL$GL$Seq))   # length as vertex attribute
					#add Length info to ggobi
					envGL$g$gf$length=V(envGL$GL$G)$length
				}else{
					gmessage("wrong sequence file",icon='error')
				}
			}
		}else{
			gmessage("No graph is loaded",icon='error')
		}
	}
	
	
	
	
	
	
	
	
	
	
	defHandlerImportContigs=function(h,...){
		aceFile=gfile(text='select ACE file', type = "open", filter = list("ACE file" = list(patterns = c("*.ace",'*.ACE')), "All files" = list(patterns = c("*"))), handler = function(h,...) return(h$file))
		if (!is.na(aceFile)) {
			setdir(aceFile)
			contigReads=getACEcontigs(aceFile)
			if (length(contigReads)>0){
				Group=ginput("Enter Group Label:",text="Contigs",icon='question')
				if (!is.na(Group)){
					# TO DO -check duplicates ##
					Subset=names(contigReads)
					Number_of_reads=as.integer(sapply(contigReads,length))
					newListDf=data.frame(Group=Group,Subset=Subset,Number_of_reads=Number_of_reads,stringsAsFactors=FALSE)
					addList(newListDf,contigReads)
					
					contigNames=rep(names(contigReads),sapply(contigReads,length))
					names(contigNames)=unlist(contigReads)
					ids_contig=paste(V(envGL$GL$G)$name,contigNames[V(envGL$GL$G)$name],sep="_")   #new name for nodes ids_contig
					rownames(envGL$g[1])=ids_contig
				}
			}
		}
	}
	
	defHanderIdslist_single=function(h,...){
		idsFile=gfile(text='select file with sequence Ids',type='open',handler=function(h,...) return(h$file))
		if (!is.na(idsFile)) {
			setdir(idsFile)
			idsList=unlist(scan(idsFile,what=character()))
			idsList=unique(idsList[idsList %in% V(envGL$GL$G)$name])
			addIdsToEditor(idsList)
			
			
		}
	}
	
	
	
	## gtable with Ids list manipulation:
	####################################################################################
	createEmptyIdsLists=function(){
		envGL$IdsListDataFrame=data.frame(Group=character(),Subset=character(),Number_of_reads=integer(),ID=integer(),stringsAsFactors=FALSE)
		envGL$IdsListIds=list()
		envGL$IdsListVisible=logical()
		envGL$IdsListFreq  # will contain frequence of Ids in list or NA if freq is not available
	}
	
	updateIdsListDataWindow=function(){
		#show actual data in data window
		gui$mainData[]=envGL$IdsListDataFrame[envGL$IdsListVisible,]
	}
	
	
	addList=function(listsDf,IdsList){  # add new list to IdsListDataFrame and also to IdsListIds and create ID 
		# lists validation - include only list with valid IDS
		validList=sapply(IdsList,FUN=function(x)any(x %in% V(envGL$GL$G)$name))
		if (any(validList)){
			listsDf=listsDf[validList,]
			IdsList=IdsList[validList]
			#check for empty lists
			notEmpty=listsDf$Number_of_reads>0
			listsDf=listsDf[notEmpty,]
			IdsList=IdsList[notEmpty]
			# check if list contain Ids frequence
			IdsListFreq=list()
			for (i in seq_along(IdsList)){
				
				if ((listsDf$Number_of_reads[i]*2)==length(IdsList[[i]])){
					IdsListFreq=append(IdsListFreq,list(as.numeric(IdsList[[i]][(1:listsDf$Number_of_reads[i])*2])))
					IdsList[[i]]=IdsList[[i]][(1:listsDf$Number_of_reads[i])*2-1]
					idsFound=which(IdsList[[i]] %in% V(envGL$GL$G)$name)
					IdsList[[i]]=IdsList[[i]][idsFound]
					listsDf$Number_of_reads[i]=length(IdsList[[i]])
				}else{
					IdsListFreq=append(IdsListFreq,NA)   # when no frequence infor is available
					idsFound=which(IdsList[[i]] %in% V(envGL$GL$G)$name)
					IdsList[[i]]=IdsList[[i]][idsFound]
					listsDf$Number_of_reads[i]=length(IdsList[[i]])
				}
				
			}
			#-------------------------	
			N=length(envGL$IdsListIds)
			M=length(IdsList)
			ID=as.integer((N+1):(N+M))
			listsDf=cbind(listsDf,ID)
			envGL$IdsListIds=append(envGL$IdsListIds,IdsList)
			envGL$IdsListFreq=append(envGL$IdsListFreq,IdsListFreq)
			
			envGL$IdsListDataFrame=rbind(envGL$IdsListDataFrame,listsDf)
			envGL$IdsListVisible=c(envGL$IdsListVisible,rep(TRUE,M))
			updateIdsListDataWindow()
		}else{galert('No list with sequences from the graph!')}
	}
	
	
	divideList=function(idsList,prefixLength=1){
		prefix=substr(idsList,1,prefixLength)
		uniquePrefixes=unique(prefix)
		idsSubLists=list()
		for (i in uniquePrefixes){
			idsSubLists[[i]]=idsList[prefix==i]
		}
		idsSubLists
	}
	
	defHandlerCreateSublist=function(h,...){
		prefixLength=as.numeric(ginput("enter prefix length"))
		if (!is.na(prefixLength)){
			selectedIdsIndex=svalue(gui$mainData,index=FALSE)
			if (length(selectedIdsIndex)>0){ #check is anything is selected
				for (j in 1:length(selectedIdsIndex)){
					selectedListsDataFrame=envGL$IdsListDataFrame[selectedIdsIndex,][j,]  
					selectedIds=envGL$IdsListIds[selectedIdsIndex][[j]]
					idsSubLists=divideList(selectedIds,prefixLength)
					for (i in 1:length(idsSubLists)){
						listsDf=selectedListsDataFrame
						listsDf[,1:2]=paste(listsDf[,1:2],"_",names(idsSubLists)[i],sep='')
						listsDf[,3]=as.integer(length(idsSubLists[[i]]))
						addList(listsDf[,-4],idsSubLists[i])
					}
				}
			}
		}
	}
	
	
	defHandlerDoubleClickHighlight=function(h,...){
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)>0){ #check is anything is selected
			selectedIds=unique(unlist(envGL$IdsListIds[selectedIdsIndex]))
			selids=V(envGL$GL$G)$name %in% selectedIds
			glyph_colour(envGL$g[1])[selids]=9
			glyph_colour(envGL$g[1])[!selids]=1
			
		}
		
	}
	
	## set color_schape and size
	setColor=function(h,...){
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)>0){ #check is anything is selected
			selectedIds=unique(unlist(envGL$IdsListIds[selectedIdsIndex]))
			if (h$action[[2]]=='color'){
				glyph_colour(envGL$g[1])[V(envGL$GL$G)$name %in% selectedIds]=h$action[[1]]	
			}
			if (h$action[[2]]=='size'){
				glyph_size(envGL$g[1])[V(envGL$GL$G)$name %in% selectedIds]=h$action[[1]]	
			}
			if (h$action[[2]]=='shape'){
				glyph_type(envGL$g[1])[V(envGL$GL$G)$name %in% selectedIds]=h$action[[1]]	
			}
			#if (length(selectedIdsIndex)==1){   # this take a long time 
			#	# add Ids group_subset label 
			#	groupName=paste(envGL$IdsListDataFrame[selectedIdsIndex,1:2],sep="_",collapse='_')
			#	ids_groupName=paste(V(envGL$GL$G)$name[V(envGL$GL$G)$name %in% selectedIds],groupName,sep="_")
			#	rownames(envGL$g[1])[V(envGL$GL$G)$name %in% selectedIds]=ids_groupName
			#}
			
		}
	}
	
	
	
	DataLoaded=function(){
		exists("GL",envir=envGL) & exists('g',envir=envGL)
	}
	
	
	ShowGGobi=function(){
		#check if ggobbi window is open 
		#if not show it
		
	}
	
	
	
	## Graphics window:
	defHandlerShowDegreeHistogram=function(h,...){
		#check if graph is loaded:
		if (DataLoaded()){
			degreeHistogram()
		}
	}
	
	defHandlerShowLengthHistogram=function(h,...){
		#check if graph is loaded:
		if (DataLoaded()){
			lengthHistogram()
		}
	}
	
	defHandlerDotter=function(h,...){
		#check if graph is loaded:
		if (DataLoaded()){
			showDotter()
		}
	}
	showDotter=function(){
		if (!is.null(envGL$GL$Seq)){
			selectedIds=V(envGL$GL$G)$name[glyph_colour(envGL$g[1])==9]
			if (length(selectedIds)!=0){
				dotter(envGL$GL$Seq[names(envGL$GL$Seq) %in% selectedIds])
				
			}else{
				gmessage('No sequences selected!')
			}
			
			
		}else{
			gmessage('No sequences loaded!')
		}
		
		
	}
	
	cleanFileName=function(filename,suffix){
		if (length(grep(paste("[.]",suffix,"$",sep=''),filename))==0){
			filename=paste(filename,".",suffix,sep='')
		}else{
			#filename is correct
		}
		filename
	}
	
	defHandlerShowDegreeLengthScatter=function(h,...){
		#check if graph is loaded:
		if (DataLoaded()){
			degreeLengthScatter()
		}
	}
	
	degreeHistogram=function(){
		#rm_gW('degreeWindow',envir=gui)
		#everything is kept inside function only pointer is put to gui enviroment
		updateGraph=function(h,...){
			plot(histogram,xlab="read degree",main="select threshold")
			abline(v=svalue(h$obj),col='red',lwd=2)
			glyph_colour(envGL$g[1])=(graph.degree>=svalue(h$obj))*8+1
		}
		
		if (exists("currentGraph",envir=gui)){ #remove it
			delete(gui$graphicsWindow,gui$currentGraph)
			rm('currentGraph',envir=gui)
		}
		graph.degree=degree(envGL$GL$G)
		gui$currentGraph=ggroup(container=gui$graphicsWindow,horizontal=FALSE,expand=TRUE)
		degreeHist=ggraphics(container=gui$currentGraph,dpi=10,ps=8)
		brks=seq(0,max(graph.degree)+1,length.out=100)
		histogram=hist(graph.degree,plot=FALSE,breaks=brks)
		graph.summary=paste("\nTotal number of vertices : ",vcount(envGL$GL$G),"\nTotal number of edges :",ecount(envGL$GL$G))
		gSl=gslider(from=min(brks)-1,to=max(brks)+1,by=1,handler=updateGraph,container=gui$currentGraph)
		infoText=glabel(text=graph.summary,container=gui$currentGraph)
		Sys.sleep(1) # time to resize window
		plot(histogram,xlab="Degree",main="select degree threshold")
		abline(v=min(brks),col='red',lwd=2)
	}
	
	
	defHandlerPlotBlastResults=function(h,...){
		#check if blast resultsexists
		if (exists('lastBlast',envir=envGL)){
			blastHistogram(envGL$lastBlast)
		}else{
			gmessage('Nothing to show',icon='error')
		}
	}
	
	defHandlerPlotFreqHistogram=function(h,...){
# check what list is selected - one only
		selectedIdsIndex=svalue(gui$mainData,index=FALSE)
		if (length(selectedIdsIndex)==1){ #check one list is selected
			freqTable=envGL$IdsListFreq[[selectedIdsIndex]]
			if (!any(is.na(freqTable))){
				names(freqTable)=envGL$IdsListIds[[selectedIdsIndex]]
				zeroFreqIds= V(envGL$GL$G)$name[!(V(envGL$GL$G)$name %in% envGL$IdsListIds[[selectedIdsIndex]])]
				zeroFreq=rep(0,length(zeroFreqIds))
				names(zeroFreq)=zeroFreqIds
				freqTable=c(freqTable,zeroFreq)
				freqTable=freqTable[V(envGL$GL$G)$name]
				plotIdsFreqHistogram(freqTable)
			}else{
				gmessage("No frequency information for selected Ids list",icon='error')
			}
		}else{
			gmessage("only one list can be selected for Frequency histogram",icon='error')
		}
	}
	
	plotIdsFreqHistogram=function(IdsFreqTable){
		
		updateGraph=function(h,...){
			plot(histogram,xlab="Frequency in Ids list",main="select Ids Frequency threshold")
			abline(v=svalue(h$obj),col='red',lwd=2)
			glyph_colour(envGL$g[1])=(graph.freq>=svalue(h$obj))*8+1
		}
		
		if (exists("currentGraph",envir=gui)){ #remove it
			delete(gui$graphicsWindow,gui$currentGraph)
			rm('currentGraph',envir=gui)
		}
		graph.freq=IdsFreqTable
		gui$currentGraph=ggroup(container=gui$graphicsWindow,horizontal=FALSE,expand=TRUE)
		freqHist=ggraphics(container=gui$currentGraph,dpi=10,ps=8)
		histogram=hist(graph.freq,plot=FALSE,breaks=50)
		graph.summary=paste("\nTotal number of vertices : ",vcount(envGL$GL$G),"\nTotal number of edges :",ecount(envGL$GL$G))
		steps=(max(graph.freq)-min(graph.freq))/100
		gSl=gslider(from=min(graph.freq),to=max(graph.freq),by=steps,handler=updateGraph,container=gui$currentGraph)
		infoText=glabel(text=graph.summary,container=gui$currentGraph)
		Sys.sleep(1) # time to resize window
		plot(histogram,xlab="Frequency in Ids list",main="select Ids Frequency threshold")
		abline(v=min(graph.freq),col='red',lwd=2)
		
		
		
	}
	
	
	
	
	blastHistogram=function(blastTable){
		updateGraph=function(h,...){
			plot(histogram,xlab="bit score",main="Blast results - select threshold")
			abline(v=svalue(h$obj),col='red',lwd=2)
			glyph_colour(envGL$g[1])=(nodesBitScore>=svalue(h$obj))*8+1
		}
		showBlastTable=function(h,...){
			#show hits for selected nodes
			selNodeNames=V(envGL$GL$G)$name[glyph_colour(envGL$g[1])==9]
			Blast2Show=blastTable[blastTable$Query %in% selNodeNames,]
			Blast2Show=Blast2Show[order(Blast2Show$bitScore,decreasing=TRUE),]
			Blast2Show=Blast2Show[!duplicated(Blast2Show$Query),]
			# show results in separate window
			if (exists("blastTableWindow",envir=gui)) {
				try(dispose(gui$blastTableWindow))
			}
			gui$blastTableWindow=ggroup(container=gwindow("Table with blast results"),horizontal=FALSE,expand=TRUE)
			gui$blastViewer=gtable(Blast2Show,container=gui$blastTableWindow,multiple=TRUE,expand=TRUE)
			gui$blast2editorButton=gbutton("move selected to Ids editor",handler=function(h,...)addIdsToEditor(svalue(gui$blastViewer)) ,container=gui$blastTableWindow)
			
			
		}
		createSublist=function(h,...){
			selNodeNames=V(envGL$GL$G)$name[glyph_colour(envGL$g[1])==9]
			Blast2Show=blastTable[blastTable$Query %in% selNodeNames,]
			Blast2Show=Blast2Show[order(Blast2Show$bitScore,decreasing=TRUE),]
			Blast2Show=Blast2Show[!duplicated(Blast2Show$Query),]
			groupID=gsub("[__].+$","",Blast2Show$Subject)
			groupReads=split(Blast2Show$Query,groupID)
			Group="Blast"  # ask fo it
			Subset=names(groupReads)
			Number_of_reads=sapply(groupReads,length)
			newListDf=data.frame(Group=Group,Subset=Subset,Number_of_reads=Number_of_reads,stringsAsFactors=FALSE)
			addList(newListDf,groupReads)
			
		}
		
		if (exists("currentGraph",envir=gui)){ #remove it
			delete(gui$graphicsWindow,gui$currentGraph)
			rm('currentGraph',envir=gui)
		}
		#remove dublicated hits keep higher score
		blastTable=blastTable[order(blastTable$bitScore,decreasing=TRUE),]
		blastTable=blastTable[!duplicated(blastTable$Query),]
		
		bitScore=blastTable$bitScore
		names(bitScore)=blastTable$Query
		#add names for no hits
		noHitNames=V(envGL$GL$G)$name[!(V(envGL$GL$G)$name %in% names(bitScore))]
		zeroBitScore=rep(0,length(noHitNames))
		names(zeroBitScore)=noHitNames
		bitScore=append(zeroBitScore,bitScore)
		
		#Score in correct order
		nodesBitScore=bitScore[V(envGL$GL$G)$name]
		
		#creat graph
		gui$currentGraph=ggroup(container=gui$graphicsWindow,horizontal=FALSE,expand=TRUE)
		blastHist=ggraphics(container=gui$currentGraph,dpi=8,ps=6)
		histogram=hist(nodesBitScore,plot=FALSE,breaks=50)
		graph.summary=paste("\nTotal number of vertices : ",vcount(envGL$GL$G),"\nTotal number of edges :",ecount(envGL$GL$G))
		gSl=gslider(from=min(nodesBitScore),to=max(nodesBitScore),handler=updateGraph,container=gui$currentGraph)
		lowerPanel=ggroup(container=gui$currentGraph,horizontal=TRUE)
		#infoText=glabel(text=graph.summary,container=lowerPanel)
		showInfoButton=gbutton("Show hits for selected\n nodes in table",container=lowerPanel,handler=showBlastTable)
		showInfoButton2=gbutton("Create Ids sublists \nfrom selected nodes",container=lowerPanel,handler=createSublist)
		Sys.sleep(1) # time to resize window
		plot(histogram,xlab="bit score",main="Blast results - select threshold")
		abline(v=min(nodesBitScore),col='red',lwd=2)
	}
	
	
	
	
	
	
	degreeLengthScatter=function(){
		if(!is.null(V(envGL$GL$G)$length)){
			if (exists("currentGraph",envir=gui)){ #remove it
				delete(gui$graphicsWindow,gui$currentGraph)
				rm('currentGraph',envir=gui)
			}
			graph.degree=degree(envGL$GL$G)
			readLength=V(envGL$GL$G)$length
			gui$currentGraph=ggroup(container=gui$graphicsWindow,horizontal=FALSE,expand=TRUE)
			degreeHist=ggraphics(container=gui$currentGraph,dpi=10,ps=8)
			Sys.sleep(1) # time to resize window
			plot(graph.degree,readLength,xlab="degree",ylab="length [bp]",col="#00000030",pch=18)
		}else{
			gmessage("No sequence length information,\n load the sequences first!",icon='error')
		}
	}
	
	
	lengthHistogram=function(){
		if(!is.null(V(envGL$GL$G)$length)){
			updateGraph=function(h,...){
				plot(histogram,xlab="read length",main="select threshold")
				abline(v=svalue(h$obj),col='red',lwd=2)
				glyph_colour(envGL$g[1])=(readLength>svalue(h$obj))*8+1
			}
			#remove current graph
			if (exists("currentGraph",envir=gui)){ #remove it
				delete(gui$graphicsWindow,gui$currentGraph)
				rm('currentGraph',envir=gui)
			}
			readLength=V(envGL$GL$G)$length
			gui$currentGraph=ggroup(container=gui$graphicsWindow,horizontal=FALSE,expand=TRUE)
			lengthHist=ggraphics(container=gui$currentGraph,dpi=10,ps=8)
			histogram=hist(readLength,plot=FALSE,breaks=50)
			graph.summary=paste("\nTotal number of vertices : ",vcount(envGL$GL$G),"\nTotal number of edges :",ecount(envGL$GL$G))
			gSl=gslider(from=min(readLength),to=max(readLength),handler=updateGraph,container=gui$currentGraph)
			infoText=glabel(text=graph.summary,container=gui$currentGraph)
			Sys.sleep(1) # time to resize window
			plot(histogram,xlab="Degree",main="select length threshold")
			abline(v=min(readLength),col='red',lwd=2)
		}else{
			gmessage("No sequence length information,\n load the sequences first!",icon='error')
		}
	}
	
	
	defHandlerBlastx=function(h,...){
		#check if graph seqs are loaded
		if(!is.null(envGL$GL$Seq)){
			oridir=getwd()
			if (!is.null(envGL$gui$wdirs$blastx)){
				setwd(envGL$gui$wdirs$blastx)
			}
			protfile=gfile(text='select fasta file with protein sequences ', type = "open",  handler = function(h,...) return(h$file))
			if (!is.na(protfile)){    # if not canceled)
				envGL$gui$wdirs$blastx=dirname(protfile)
				protseq=read.BStringSet(protfile,format='fasta')
				cat("running blastx\n")
				blastxresults=blastall(protseq,envGL$GL$Seq,programname='blastx',m=8,params="",protein=TRUE)
				# params=' -m 8 -p blastx -a 7 -b 1 -v 1 -K 1'
				# put bitScore results into histogram
				blastHistogram(blastxresults)
				envGL$lastBlast=blastxresults
			}
			setwd(oridir)
		}else{
			gmessage("No sequence information,\n load the sequences first!",icon='error')
		}
	}
	
	defHandlerBlastn=function(h,...){
		#check if graph seqs are loaded
		if(!is.null(envGL$GL$Seq)){
			oridir=getwd()
			if (!is.null(envGL$gui$wdirs$blastn)){
				setwd(envGL$gui$wdirs$blastn)
			}
			DNAfile=gfile(text='select fasta file with DNA sequences ', type = "open",  handler = function(h,...) return(h$file))
			if (!is.na(DNAfile)){    # if not canceled)
				envGL$gui$wdirs$blastn=dirname(DNAfile)
				DNAseq=read.DNAStringSet(DNAfile,format='fasta')
				cat("running blastn\n")
				blastnresults=blastall(DNAseq,envGL$GL$Seq,programname='blastn',m=8,params="",protein=FALSE)
				# params=' -m 8 -p blastx -a 7 -b 1 -v 1 -K 1'
				# put bitScore results into histogram
				blastHistogram(blastnresults)
				envGL$lastBlast=blastnresults
			}
			setwd(oridir)
		}else{
			gmessage("No sequence information,\n load the sequences first!",icon='error')
		}
	}
	
	defHandlerImportBlastTable=function(h,...){
		blastFile=gfile('Select file with blast results in tabular format',type='open')
		if (!is.null(blastFile)){   # if not canceled
			setdir(blastFile)
			blast=read.table(blastFile,sep="\t",header=FALSE,as.is=TRUE)
			#check blast validity:
			cond1=dim(blast)[2]==12   # corect dimension
			cond2=blast[,1] %in% V(envGL$GL$G)$name  # to check ids
			if (cond1 & sum(cond2)>0){
				blast=blast[cond2,]
				colnames(blast)=c('Query','Subject','percIdentity','alignmentLength','mistmatches','gapOpenings','q.start','q.end','s.start','s.end','Evalue','bitScore')		
				blastHistogram(blast)
				envGL$lastBlast=blast	
			}
		}
	}
	
	
	defHandlerExportBlastTable=function(h,...){
		if (!is.null(envGL$lastBlast)){
			blastFile=gfile("Save blast results as..",type='save')
			if (!is.null(blastFile)){   # if not canceled
				setdir(blastFile)
				write.table(envGL$lastBlast,file=blastFile,sep='\t',quote=F,col.names=F,row.names=F)
			}
		}else{
			galert('no blast result!',icon='error')
		}
	}
	
	subGL=function(GL,ids){
		#make subgraph from GL based on Ids list
		allIds=V(GL$G)$name
		include=which(allIds %in% ids)
		subGL=list()
		subGL$G=induced.subgraph(GL$G,include)
		subGL$L=GL$L[include,]
		subGL
	}
	
	defHandlerSaveSubgraph=function(h,...){
		if (DataLoaded()) {
			selectedIds=V(envGL$GL$G)$name[glyph_colour(envGL$g[1])==9]
			GL=subGL(envGL$GL,selectedIds)
			GLfile=gfile("save selected Subgraph as..", type='save',handler = function(h,...) return(h$file))
			setdir(GLfile)
			GLfile=paste(gsub("[.]GL$","",GLfile),'.GL',sep='')
			save(GL,file=GLfile)
		}
	}
	
	
	defHandlerNewLayout=function(h,...){
		newL=layout.fruchterman.reingold(envGL$GL$G,dim=3,verbose=TRUE)
		N=dim(envGL$g$gf)[2]
		envGL$g$gf[[paste('V',N+1,sep='')]]=newL[,1]
		envGL$g$gf[[paste('V',N+2,sep='')]]=newL[,2]
		envGL$g$gf[[paste('V',N+3,sep='')]]=newL[,3]
	}
	
	plotGL=function(GL,edge.color=NULL,vertex.color=NULL,wlim=NULL,cex=.8,...){
		GG=GL$G
		LL=GL$L[,1:3]
		e=get.edgelist(GG,names=FALSE)
		w=E(GG)$weight
		if (!is.null(wlim)) {e=e[w>wlim,]; w=w[w>wlim]}
		X0=LL[e[,1],1]
		Y0=LL[e[,1],2]
		X1=LL[e[,2],1]
		Y1=LL[e[,2],2]
		if (is.null(edge.color))
		{edge.color="#AAAAAA"}
		if (is.null(vertex.color))
		{vertex.color="#FF8888"}
		op=par(mar=c(0,0,0,0))
		plot(range(LL[,1]),range(LL[,2]),xlab="",ylab="",axes=FALSE,type="n",...)
		segments(X0,Y0,X1,Y1,lwd=.5,col=edge.color)
		LL1=data.frame(LL,vertex.color,stringsAsFactors=FALSE)
		LL2=data.frame(LL,'#000000',stringsAsFactors=FALSE)
		colnames(LL1)=colnames(LL2)=c('v1','v2','v3','color')
		LL=rbind(LL2,LL1)
		LL=LL[order(LL[,3]),]
		points(LL[,1:2],pch=c(1,19),col=LL$color,cex=cex)
		par(op)
	}
	
	ggobi_current_layout=function(g,L){
		gd=displays(g)[[1]]
		#check if the mode is correct:
		if (pmode(gd) %in% c('Rotation','2D Tour','2x1D Tour')){
			pr=ggobi_display_get_tour_projection(gd)
			if (dim(pr)[1]>3){
				pr=pr[1:3,]		#remove other dimensions
			}
			x1=sum(pr[,1])
			x2=sum(pr[,2])
			x3=x1-pr[1,1]+pr[1,2]
			mat=rbind(t(pr),c(x1,x2,x2)-pr[,1]-pr[,2])
			L2=L
			L2[,1]=mat[1,1]*L[,1]+mat[1,2]*L[,2]+mat[1,3]*L[,3]
			L2[,2]=mat[2,1]*L[,1]+mat[2,2]*L[,2]+mat[2,3]*L[,3]
			L2[,3]=mat[3,1]*L[,1]+mat[3,2]*L[,2]+mat[3,3]*L[,3]
			
			return(L2)
		}else{
			return(L)
		}
	}
	
	getColors=function(g){
		clrs=colorscheme(g)$colors
		col=sapply(clrs,FUN=function(x)rgb(x[1],x[2],x[3],1))
		colors.rbg=col[glyph_colour(g[1])]
		colors.rbg
	}
	
	defHandlerSaveImage=function(h,...){
		imageFile=gfile("Save as",type="save")
		Lnow=ggobi_current_layout(envGL$g,envGL$GL$L)
		GLnow=envGL$GL
		GLnow$L=Lnow
		png(paste(imageFile,".png",sep=""),width=1500,height=1500)
		plotGL(GLnow,vertex.color=getColors(envGL$g))
		dev.off()
		
	}
	
	defHandlerSaveCurrentProjection=function(h,...){
		Lnow=ggobi_current_layout(envGL$g,envGL$GL$L)
		GL=envGL$GL
		GL$L=Lnow
		glFile=gfile("Save as",type="save")
		save(GL,file=paste(glFile,".GL",sep=''))
	}
	
	
	
	
	
	# to remove gWidget and its corresponding var stored in envir
	rm_gW=function(objname,envir){
		obj=try(get(objname,envir=envir,inherits=FALSE),silent=TRUE)
		try(dispose(obj),silent=TRUE)
		if (exists(objname)){
			try(rm(list=objname,envir=envir),silent=TRUE)
		}
	}
	
	
# to remove ggobi windows and other individual windows and clean up envGL enviroment with data
	removeAll=function(){
		if (exists("blastTableWindow",envir=gui)) {
			try(dispose(gui$blastTableWindow))
		}
		
		if (exists('g',envir=envGL)){
			while(ggobi_count()>0) close(ggobi_get())  # will close all ggobi instances
			rm_gW('degreeWindow',envir=gui)   # will remove degree Window if open
			rm(list=ls(envir=envGL),envir=envGL)
		}
		
	}
	
	defHandleReopenGgobi=function(h,...){
		if (DataLoaded()){
			if (gconfirm("This will close all existing ggobi windows!")){
				if (exists('g',envir=envGL)){
					while(ggobi_count()>0) close(ggobi_get())  # will close all ggobi instances
				}
				envGL$g=showGraph(envGL$GL)
				#add degree info to ggobi
				envGL$g$gf$degree=degree(envGL$GL$G)
				if (!is.null(V(envGL$GL$G)$length)){
					envGL$g$gf$length=V(envGL$GL$G)$length
				}
			}
		}else{
			gmessage('No data loaded!',icon='error')
		}
	}
	
	
	
	
	
	
	defHandlerSaveProject=function(h,...){
		prFile=gfile("Save project as...",type='save')
		setdir(prFile)
		prFile=cleanFileName(prFile,suffix="SGR")
		save(list=ls(envir=envGL),envir=envGL,file=prFile)
	}
	
	defHandlerOpenProject=function(h,...){
		prFile=gfile("Select project file",type='open')
		if (!is.na(prFile)){
			setdir(prFile)
			#clean project enviroment
			removeAll()
			load(envir=envGL,file=prFile)
			#check validity of file
			if (exists("GL",envir=envGL)){
				if (class(envGL$GL$G)=='igraph' & is.matrix(envGL$GL$L) & is.list(envGL$IdsListIds) & is.data.frame(envGL$IdsListDataFrame) & !is.null(envGL$IdsListVisible)){
					#update widgets:
					envGL$g=showGraph(envGL$GL)
					if (length(names(envGL$g))==2){# do not show histograpm when no edges are not imported to ggobi
						degreeHistogram()
						#add degree info to ggobi
						envGL$g$gf$degree=degree(envGL$GL$G)
					}else{
						plot.new()
						glyph_type(envGL$g[1])=1
					}
					
					
					updateIdsListDataWindow()
					svalue(gui$StatusBar)=paste(prFile," :  V=",vcount(envGL$GL$G)," E=",ecount(envGL$GL$G),sep='')
				}else{
					gmessage("Not a valid SeqGrapheR project",icon='error')
				}
			}else{
				gmessage("Not a valid SeqGrapheR project",icon='error')
			}
		}
	}
	
	
	
	newIdsList=function(h,...){
		#click
		
	}
	
	
	defHandler=function(h,...){
		print(h)
	}
	
	defHandlerQuit=function(h,...){
		removeAll()
		try(dispose(gui$winMain),silent=TRUE)
		#rm(list=c('gui','envGL'))
	}
	
	#########################################################################
	################################ GUI MAIN  ##############################
	#########################################################################
	
	
	iconDir=paste(system.file('images',package='SeqGrapheR'),'/',sep='')   # for package
	
# for testing as source:
	#iconDir="/home/petr/Dropbox/Rscripts/icons"  # home
	#iconDir="/mnt/raid/petr/icons"	#work						# work
	#iconDir="/mnt/raid/spolecny/petr/icons"   
	
	
	
	## enviroments for DATA
	envGL <- new.env()   # enviromenrt for graph object and related ids lists:
	## enviroment for GUI
	gui <- new.env()
	
	createEmptyIdsLists()    # creates envGL$IdsListDataFrame and	envGL$IdsListIds
	
	gui$winMain = gwindow("SeqGrapheR",visible=FALSE)
	size(gui$winMain)=c(900,700,exapand=TRUE)
	
	
	
	gui$groupMenu = ggroup(container = gui$winMain, horizontal = FALSE)  # main group
	
	## menu list:
	mbl = list()
	# open menu - import graph from various format
	mbl$File$"Open graph"$"GL format"$handler=defHandlerOpenGL   #done
	mbl$File$"Open graph"$"ncol format"$handler=defHandlerOpenNcol  #done
	mbl$File$"Open graph"$"from DNA sequences"$handler=defHandlerOpenSequences2Graph  #done
	
	mbl$File$Project$Open$handler=defHandlerOpenProject
	mbl$File$Project$Save$handler=defHandlerSaveProject
	
	mbl$File$Import$Sequences$handler=defHandlerImportSequences  #done
	mbl$File$Import$"Single Ids list to editor"$handler=defHanderIdslist_single   #done
	mbl$File$Import$"Single Ids (tabular format)"$handler=defHandlerLoadIdsTable #done
	mbl$File$Import$"Multiple Ids lists"$handler=defHandlerLoadMultipleLists #done
	
	mbl$File$Import$Contigs$handler=defHandlerImportContigs			#done
	mbl$File$Import$"Blast results (tabular format)"$handler=defHandlerImportBlastTable
	#mbl$File$Import$"Color scheme"$handler=defHandler
	
	mbl$File$Export$"Selected Ids lists"$handler=defHandlerSaveMultipleLists
	mbl$File$Export$"Current graph projection"$handler=defHandlerSaveCurrentProjection#
	mbl$File$Export$"Graph image"$handler=defHandlerSaveImage   #done
	mbl$File$Export$"Selected nodes as subgraph"$handler=defHandlerSaveSubgraph
	mbl$File$Export$"Blast results"$handler=defHandlerExportBlastTable
	
	
	mbl$Tools$Plot$"Read degree"$handler=defHandlerShowDegreeHistogram #done
	mbl$Tools$Plot$"Read length"$handler=defHandlerShowLengthHistogram  #done
	mbl$Tools$Plot$"Lengths vs. degree"$handler=defHandlerShowDegreeLengthScatter #done
	mbl$Tools$Plot$"Frequency of Ids in selected list"$handler=defHandlerPlotFreqHistogram
	mbl$Tools$Plot$"Blast results"$handler=defHandlerPlotBlastResults   #done
	
	
	mbl$Tools$Similarity_search$'Protein blastx'$handler=defHandlerBlastx
	mbl$Tools$Similarity_search$'DNA blastn'$handler=defHandlerBlastn
	mbl$Tools$'Graph analysis'$Layout$New$handler=defHandlerNewLayout
	mbl$Tools$'Show dotplot for selected sequences'$handler=defHandlerDotter
	
	mbl$Tools$'Graph analysis'$'Find clusters'$handler=defHandlerShowClusters
	mbl$Tools$'Graph analysis'$'Show neighbors'$handler=defHandlerShowNeighbors
	
	mbl$GGobi$Reopen$handler=defHandleReopenGgobi
	
	mbl$Quit$Quit$handler=defHandlerQuit
	addhandlerdestroy(gui$winMain,handler=defHandlerQuit)
	gui$mb = gmenu(mbl, container=gui$groupMenu)
	gui$desktop= gframe(container=gui$groupMenu,horizontal=TRUE,expand=TRUE)
	
	
	## group for data window 
	#gui$gdl = glayout(cont = gui$groupMenu,homogeneous=TRUE,spacing=1)
	
	
	
	
	
	## Table for Ids selection
	gui$gfb=gframe(text=" Ids lists selector:",container=gui$desktop,horizontal=FALSE,expand=TRUE)
	gui$mainData=gtable(envGL$IdsListDataFrame,multiple=TRUE,container=gui$gfb,expand=TRUE,chosencol=4,handler=defHandlerDoubleClickHighlight)
	
	## bottons for Ids list manipulation
	gui$gdc=ggroup(container=gui$desktop,expand=FALSE,horizontal=FALSE)
	glabel(" Ids/lists manipulation",container=gui$gdc)
	gui$rightArrow=gbutton("list(s) ---> editor",handler=defHandlerAdd2Editor,container=gui$gdc,size='button')
	gui$leftArrow=gbutton("Remove selected lists \n<--from editor",handler=defHandlerRemoveFromEditor,container=gui$gdc,size='button')
	gui$newList=gbutton("New list\nfrom editor",handler=defHandlerMakeNewList,container=gui$gdc,size='button')
	gui$deleteList=gbutton("Delete selected list!",handler=defHandlerDeleteSelectedList,container=gui$gdc,size='button')
	
	# yellow arrow - brush selection
	gui$move2editorGroup=gframe(container=gui$gdc,expand=FALSE,horizontal=FALSE)
	gui$move2editorLabel=glabel("Move selected yellow\nnodes into editor",handler=defHandlerYellow2editor,container=gui$move2editorGroup)
	gui$move2editorArrow=gimage(filename="yellowArrow.png",dirname=iconDir,handler=defHandlerYellow2editor,container=gui$move2editorGroup)
	
	gui$makeListSubsets=gbutton("Make sublist from\n selected Ids list\n using Ids prefix",handler=defHandlerCreateSublist,container=gui$gdc,size='button')
	gui$clearEditor=gbutton("Clear editor",handler=defHandlerClearEditor,container=gui$gdc,size='button')
	addSpace(gui$gdc,10,horizontal=FALSE)
	
	## frame  for node Style selections
	gui$gfa=gframe(text=" node style\n options:",container=gui$gdc,horizontal=FALSE,expand=FALSE)
	
	## icons for color selection:
	addSpace(gui$gfa,5,horizontal=FALSE)
	gui$clrLab=glabel("Set color :   ",container=gui$gfa,expand=FALSE)
	gui$colors=glayout(container=gui$gfa,border=TRUE,homogeneous=FALSE,expand=FALSE,spacing=1)
	iconsColor=c("color1.png", "color2.png", "color3.png", "color4.png", "color5.png","color6.png", "color7.png", "color8.png", "color9.png")
	iconsPositions=matrix(c(1,2,3,1,2,3,1,2,3,1,1,1,2,2,2,3,3,3),byrow=FALSE,ncol=2)
	iconsPositions=matrix(c(1,2,3,4,5,6,7,8,9,1,1,1,1,1,1,1,1,1),byrow=FALSE,ncol=2)
	for (i in 1:9){
		gui$colors[iconsPositions[i,2],iconsPositions[i,1]]=assign(iconsColor[i],value=gimage(filename=iconsColor[i],dirname=iconDir,handler=setColor,container=gui$colors,action=list(i,'color')),envir=gui)
	}
	gui$color22list2=gbutton("<---colors to Ids lists",handler=defHandlerColors2lists,container=gui$gfa,size='button')
	
	addSpace(gui$gfa,20,horizontal=FALSE)
	## icons for shape selection:
	gui$shpLab=glabel("Set shape :   ",container=gui$gfa,expand=FALSE)
	gui$shapes=glayout(container=gui$gfa,border=TRUE,homogeneous=FALSE,expand=FALSE,spacing=1)
	iconsShape=c("shape1.png", "shape2.png", "shape3.png", "shape4.png", "shape5.png","shape6.png", "shape7.png")
	for (i in 1:7){
		gui$shapes[iconsPositions[i,2],iconsPositions[i,1]]=assign(iconsShape[i],value=gimage(filename=iconsShape[i],dirname=iconDir,handler=setColor,container=gui$colors,action=list(i,'shape')),envir=gui)
	}
	
	addSpace(gui$gfa,5,horizontal=FALSE)
	## icons for size selection:
	gui$sizLab=glabel("Set size :   ",container=gui$gfa,expand=FALSE)
	gui$sizes=glayout(container=gui$gfa,border=TRUE,homogeneous=FALSE,expand=FALSE,spacing=1)
	iconsSize=c("size1.png", "size2.png", "size3.png", "size4.png", "size5.png","size6.png", "size7.png",'size8.png')
	for (i in 1:8){
		gui$sizes[iconsPositions[i,2],iconsPositions[i,1]]=assign(iconsSize[i],value=gimage(filename=iconsSize[i],dirname=iconDir,handler=setColor,container=gui$colors,action=list(i,'size')),envir=gui)
	}
	
	
	## Ids list editor:
	gui$gfd=gframe(container=gui$desktop,horizontal=FALSE,expand=TRUE)
	gui$gfdEditor=gframe(text="Ids editor",container=gui$gfd,horizontal=FALSE,expand=TRUE)
	gui$IdsEditor=gtext(container=gui$gfdEditor,expand=TRUE,horizontal=FALSE)
	gui$graphicsWindow=gframe(container=gui$gfdEditor,horizontal=FALSE,expand=TRUE)  # window to put various plots in
	gui$StatusBar=gstatusbar(text = "", container = gui$winMain,expand=TRUE)
	
	size(gui$gfb)[1]=0
	size(gui$graphicsWindow)[2]=100
	size(gui$gfdEditor)[1]=0
	visible(gui$winMain)=TRUE
	output=list(envGL=envGL,gui=gui)
}
