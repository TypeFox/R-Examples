#' Index of paired faces with directions
#'
#' This does some heavy lifting to pull together faces which are paired together. This is needed for many later functions for compiling OPC
#' @param plyFile a stanford PLY file
#'
#' 
#' index_paired_directed_faces()

index_paired_directed_faces <- function(plyFile) {

Faces <- t(plyFile$it) ### Extract the face list from the plyFile
rownames(Faces) <- seq(1:length(Faces[,1])) ## Add names to the face list to keep them strait
Dbins <- plyFile$Directional_Bins ### Extracting the Binned Directions

Faces <- cbind(Faces, Dbins) ### Make a faces matrix to manipulate

binlist <- list() ### Make a list for holding the separete directional bins

### Sort into binlist on basis of face direction
for (i in unique(Faces[,4])) {
	
	binlist[[as.character(i)]] <- subset(Faces[,1:3], Faces[,4]==as.numeric(i))	
}

### Order the vertex reference numbers from smallest to largest###
for (i in names(binlist)) {
	binlist[[as.character(i)]] <- t(apply(binlist[[as.character(i)]],1,sort))
	
}

direct_face_lists <- list() ### Create Edge List for holding edge names per face

### into the rabbit hole ###
for (i in names(binlist)) {
	
	Direction <- binlist[[as.character(i)]]
	paired_face_list <- list()
	for (j in rownames(Direction)) {
		Caterpiller1 <- Direction[as.character(j),][1]
		Caterpiller2 <- Direction[as.character(j),][2]
		Caterpiller3 <- Direction[as.character(j),][3]
		
		QofH1 <- paste(Caterpiller1, Caterpiller2, sep='_')
		QofH2 <- paste(Caterpiller1, Caterpiller3, sep='_')
		QofH3 <- paste(Caterpiller2, Caterpiller3, sep='_')
		paired_face_list[[j]] <- paste(QofH1, QofH2, QofH3, sep=",")
	}
	
	for (j in names(paired_face_list)) {
		paired_face_list[[j]] <- unlist(strsplit(paired_face_list[[j]], ','))
	}
	
	whiterabbit <- matrix(0, length(names(paired_face_list)), 3)
	rownames(whiterabbit) <- names(paired_face_list)
	for (j in names(paired_face_list)) {
		cheshirecat <- paired_face_list[[as.character(j)]]
		whiterabbit[j,1] <- cheshirecat[1]
		whiterabbit[j,2] <- cheshirecat[2]
		whiterabbit[j,3] <- cheshirecat[3]
	}
	alice <- c(whiterabbit[,1], whiterabbit[,2], whiterabbit[,3])
	
	madhatter <- data.frame(Faces=names(alice), edges=alice)
	madhatter <- aggregate(madhatter, list(madhatter$edges), FUN=paste)
	madhatter <- as.matrix(madhatter) 
	madhatter <- madhatter[,-3]
	colnames(madhatter) <- c('Edge', 'Faces')
	
	direct_face_lists[[as.character(i)]] <- madhatter
}
### out of the rabbit hole ###

index_paired_directed_faces <- direct_face_lists
return(index_paired_directed_faces)
}

