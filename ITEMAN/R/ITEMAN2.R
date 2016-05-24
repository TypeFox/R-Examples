#######################################################################################
#
#     An R Routine for Item Analysis
#
#     Cengiz Zopluoglu
#     Assistant Professor
#     Research, Measurement, and Evaluation
#     University of Miami
#
#     Please report any programming bug: c.zopluoglu1@miami.edu
# 
#     09/29/2015
#
#     Depends on: "ggplot"
#     Depends on: "polycor"
##########################################################################################


##########################################################################################
#
#
#   PRESS CTRL+A to highlight everything and hit ENTER
#  
#   Then, close this file, and open "SampleItemAnalysis.r" for an illustration
#
##########################################################################################


ITEMAN2 <- function (data, options,ngroup=ncol(data)+1,correction=TRUE) 
{

      #########################################################################################
	# data, a data frame with N rows and n columns, where N denotes the number of subjects 
      #        and n denotes the number of items. All items should be scored using nominal/ordinal 
      #        response categories. All variables (columns) must be "numeric". 

  # options, numbers representing the response categories (e.g.,0,1,2,3)
	#          make sure each item is consistent, and includes the same response options
  # Recommend that the numerical codes are recoded such that the minimum score is 0

      # ngroup, number of score groups
	#
	# correction, TRUE or FALSE, if TRUE item and distractor discrimination is corrected for
	# 		  spuriousnes by removing the item score from the total score
      #########################################################################################


    total.score <- rowMeans(data,na.rm=TRUE)*ncol(data)

    pbis <- c()
    pbis.corrected <- c()
    bis  <- c()
    bis.corrected <- c()

	for(k in 1:ncol(data)) { 
		pbis[k]=cor(data[,k],total.score,use="pairwise.complete.obs")
		pbis.corrected[k]=cor(data[,k],
		                      rowMeans(data[,-k],na.rm=TRUE)*(ncol(data)-1),
                          use="pairwise.complete.obs")
		bis[k]=polyserial(total.score,data[,k])
    bis.corrected[k]=polyserial(rowMeans(data[,-k],na.rm=TRUE)*(ncol(data)-1),data[,k])
	}


    item.stat <- matrix(nrow=ncol(data),ncol=4)
	  colnames(item.stat) <- c("Mean Score","Item Difficulty","Point-Biserial","Polyserial")

    rnames <- ("Item 1")
    for(i in 2:ncol(data)){ rnames <- c(rnames,paste("Item ",i,sep=""))}
    rownames(item.stat) <- rnames	
	    item.stat[,1]=colMeans(data,na.rm=TRUE)
    item.stat[,2]=colMeans(data,na.rm=TRUE)/max(options)
    if(correction==TRUE){ item.stat[,3]=pbis.corrected } else { item.stat[,3]=pbis }
	  if(correction==TRUE){ item.stat[,4]=bis.corrected } else { item.stat[,4]=bis }


    sgroups <- cut(total.score,breaks=ngroup)
    slevels <- levels(sgroups)

    sgnum <- rowMeans(cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", slevels) ),
                   upper = as.numeric( sub("[^,]*,([^]]*)\\]","\\1",slevels))))

	

    SG <- vector("list",ngroup)
    
    for(j in 1:ngroup){
	SG[[j]]=which(sgroups==slevels[j])
    }
 
    prop <- vector("list",ncol(data))
    names(prop) <- rnames
   
    for(i in 1:ncol(data)) {

	dist <- matrix(nrow=length(options),ncol=ngroup)
	colnames(dist) <- slevels
	rownames(dist) <- options 

	for(g in 1:ngroup){
	  for(o in 1:length(options)){
		dist[o,g]=length(which(data[SG[[g]],i]==options[o]))/length(SG[[g]])
	  }
	}

	prop[[i]]=dist

    }

    dist.sel <- matrix(nrow=ncol(data),ncol=length(options))  
    dist.disc <- matrix(nrow=ncol(data),ncol=length(options))
    dist.disc2 <- matrix(nrow=ncol(data),ncol=length(options))
    colnames(dist.disc) <- options
    rownames(dist.disc) <- rnames
    colnames(dist.disc2) <- options
    rownames(dist.disc2) <- rnames
    colnames(dist.sel) <- options
    rownames(dist.sel) <- rnames

	for(i in 1:ncol(data)){
	  for(o in 1:length(options)) {      
		  temp <- ifelse(data[,i]==options[o],1,0)
		  temp[is.na(temp)]=0
		  dist.sel[i,o]=mean(temp,na.rm=TRUE)
		  if(correction==FALSE){
		    dist.disc[i,o]=cor(temp,total.score,use="pairwise.complete.obs")
		    dist.disc2[i,o]=polyserial(total.score,temp)
		  } else {
		    dist.disc[i,o]=cor(temp,rowMeans(data[,-i],na.rm=TRUE)*(ncol(data)-1),use="pairwise.complete.obs")
		    dist.disc2[i,o]=polyserial(rowMeans(data[,-i],na.rm=TRUE)*(ncol(data)-1),temp)
		  }
	  }
    
	 }
	

	plots <- vector("list",ncol(data))

	for(i in 1:ncol(data)) {

		options.d <- c()
		for(u in 1:length(options)){ 
			if(correction==TRUE){
				options.d[u] <- paste(options[u],"( ",round(dist.disc2[i,u],2)," )",sep="")
			} else { options.d[u] <- paste(options[u],"( ",round(dist.disc[i,u],2)," )",sep="") }
		}
		
		d <- as.data.frame(cbind(sg=sgnum,p=prop[[i]][1,]))
		for(u in 2:length(options)){ d <- rbind(d,cbind(sg=sgnum,p=prop[[i]][u,]))}
		optt <- c()
		for(u in 1:length(options)){ optt <- c(optt,rep(options.d[u],ngroup))}
		d$opt <- optt
		
		
		pp <- ggplot(data=d,aes_string(x="sg",y="p",group="opt",shape="opt"))+
                 geom_line()+
                 geom_point(size=3)+
		             ggtitle(paste("Item ",i,sep=""))+
                 theme(panel.background = element_blank(),legend.title=element_blank(),legend.key = element_blank())+
		             scale_x_continuous(limits = c(0,ncol(data)*max(options)),breaks=seq(0,ncol(data)*max(options),ceiling(ncol(data)/10)))+
		             scale_y_continuous(limits = c(0,1))+xlab("Score Groups")+ylab("Proporion of Being Selected")
       #  theme(legend.justification=c(0,1),legend.position=c(0,1),legend.text=element_text(size=12,face="bold"))
		
		plots[[i]] <- pp
	}

	###############################################################

	cat("************************************************************************","\n")
	cat("ITEMAN: An R routine for Classical Item Analysis","\n")
	cat("","\n")
	cat("Cengiz Zopluoglu","\n")
	cat("","\n")
	cat("University of Miami","\n")
	cat("Department of Educational and Psychological Studies","\n")
	cat("Research, Measurement, and Evaluation Program","\n")
	cat("","\n")
	cat("c.zopluoglu@miami.edu","\n")
	cat("","\n")
	cat("Please report any programming bug or problem you experience to improve the code.","\n")
	cat("*************************************************************************","\n")

	cat("Processing Date: ",date(),"\n")

	cat(sprintf("%50s","ITEM STATISTICS"),"\n")
	cat("","\n")
	print(round(item.stat,3))
	cat("","\n")
	cat("    * Item difficulty is the ratio of mean score to possible maximum score","\n")
  cat("      and assumes the minimum score is 0","\n")
  cat("","\n")
	cat("","\n")

	cat(sprintf("%50s","RESPONSE CATEGORY SELECTION PROPORTIONS"),"\n")
	cat("","\n")
	print(round(dist.sel,3))
	cat("","\n")
	cat("","\n")
	cat("","\n")

  cat(sprintf("%50s","RESPONSE CATEGORY Point-Biserial Correlation"),"\n")
	cat("","\n")
	print(round(dist.disc,3))
	cat("","\n")
	cat("","\n")
  
	cat(sprintf("%50s","RESPONSE CATEGORY Biserial Correlation"),"\n")
	cat("","\n")
	print(round(dist.disc2,3))
	cat("","\n")
	cat("","\n")
	cat("","\n")
	return(list(plots=plots))
}







