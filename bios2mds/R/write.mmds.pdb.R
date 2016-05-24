write.mmds.pdb <- function (x, project = NULL, axis = c(1, 2, 3), file.pdb = "R.pdb", file.pml = NULL) {

  if (!inherits(x, "mmds"))
    stop("object of class 'mmds' expected")
  if (length(axis) != 3)
    stop("axis have to be of length 3")
  if (any(axis > length(x$eigen.perc)))
    stop("wrong axis")
  pdb.coord <- function(i,nbgroup, x, y, z) {
    format <- "%s%5s%3s%6s%3s%-7s%8.3f%8.3f%8.3f%6s%6s%12s"
    sprintf(format, "HETATM", i, "O", "HOH", paste(groupId[nbgroup]," ",sep="") , i, x, y, z, "1.00", "0.00", "O")
  } 

  pdb.axis <- function(i, x, y, z, atom) {
    format <- "%s%5s%3s%9s%-7s%8.3f%8.3f%8.3f%6s%6s%12s"
    sprintf(format, "HETATM", i, atom, "Z ", i, x, y, z, "1.00", "0.00", atom)
  }

  pdb.conect <- function(i, j) {
    format <- "%s%5s%5s"
    sprintf(format, "CONECT", i, j)
  } 
  groupName<-"NoGroup"
  groupId<-c("Z","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P",
	"Q","R","S","T","U","V","W","X","Y")
    
  x$coord <- x$coord * 50
  active.nb <- nrow(x$coord)
  ncol <- nrow(x$col)
  if (!is.null(project) && inherits(project,"project")) {
    project$coord <- project$coord * 50
    sup.nb <- nrow(project$coord)
  }
 resInGroup<-list()
  resInGroup$NoGroup<-vector()

  if(!is.null(file.pml)){
   pml.lines <- NULL
  pml.header <- "bg_color white"
  active.names<-rownames(x$coord)
  if(x$col[1]==""){
	  active.col<-matrix(rep(as.vector(col2rgb(c("blue"))/255),length(x$coord[,1])),ncol=3,byrow=TRUE)
	}
  else	{
	active.col<-matrix(as.vector(col2rgb(x$col[,3])/255),ncol=3,byrow=TRUE)
  }
  if (!is.null(project) && inherits(project,"project")) {
  sup.names<-rownames(project$coord)
  if(project$col[1]==""){
	  sup.col<-matrix(rep(as.vector(col2rgb(c("red"))/255),length(project$coord[,1])),ncol=3,byrow=TRUE)
	}
  else	{
	sup.col<-matrix(as.vector(col2rgb(project$col[,3])/255),ncol=3,byrow=TRUE)
  }
  }
}

  pdb.lines <- NULL

  #pdb for active
  j <- 1
    for (i  in 1:nrow(x$group)) {
	if(!length(which(x$group[i]==groupName))==1){
		groupName<-c(groupName,x$group[i])
		resInGroup[[x$group[i]]]<-vector()
        }
    }
  for (i in 1:active.nb) {
    if(ncol!=1){
      nbgroup=which(x$col[i,2]==groupName)
    }
    else {
      nbgroup=which(x$col[1,2]==groupName)
    }
    if(nbgroup > 26){
	nbgroup<-1
    }
    pdb.lines <- rbind(pdb.lines, pdb.coord(j,nbgroup,x$coord[i, axis[1]], x$coord[i, axis[2]],
       x$coord[i, axis[3]]))
 if(!is.null(file.pml)){
    k<-j
   
    pml.lines<-rbind(pml.lines, paste("set_color ref_", groupId[nbgroup],k, ", [", active.col[i,1],
    ", ", active.col[i,2], ", ", active.col[i,3], "]", sep = ""))
	pml.lines <- rbind(pml.lines, paste("label (///",groupId[nbgroup],"/",k, "), ", "\"", active.names[i], "\"", sep = ""))
    pml.lines <- rbind(pml.lines, paste("color ref_", groupId[nbgroup],k, ", ///",groupId[nbgroup],"/",k, sep = ""))
    if(ncol!=1){    
    pml.lines <- rbind(pml.lines, paste("select ref_",x$col[i,2],"_",active.names[i],",(///",groupId[nbgroup],"/",k,")",sep=""))
    resInGroup[[x$col[i,2]]]<-c(resInGroup[[x$col[i,2]]],j)
	}
else {
    pml.lines <- rbind(pml.lines, paste("select ref_",x$col[1,2],"_",active.names[i],",(///",groupId[nbgroup],"/",k,")",sep=""))
    resInGroup[[x$col[1,2]]]<-c(resInGroup[[x$col[1,2]]],j)
}
}
      j <- j + 1
  }
if(!is.null(file.pml)){
  pml.header <- rbind(pml.header,paste("select ref, resi 1-",j-1,sep=""))
}     
  #pdb for sup if given
  if (!is.null(project) && inherits(project,"project")) {
    ncol <- nrow(project$col)
   for (i  in 1:nrow(project$group)) {
	if(!length(which(project$group[i]==groupName))==1){
		groupName<-c(groupName,project$group[i])
		resInGroup[[project$group[i]]]<-vector()
        }
    }
    j <-5001
    for (i in 1:sup.nb) {
      if(ncol!=1){
	nbgroup=which(project$col[i,2]==groupName)
      }
      else {
     
	nbgroup=which(project$col[1,2]==groupName)
      }
      if(nbgroup > 26){
  	nbgroup<-1
      }
      pdb.lines <- rbind(pdb.lines, pdb.coord(j,nbgroup, project$coord[i, axis[1]], project$coord[i, axis[2]],
         project$coord[i, axis[3]]))
 if(!is.null(file.pml)){
    pml.lines<-rbind(pml.lines, paste("set_color sup_", groupId[nbgroup],j, ", [", sup.col[i,1],
    ", ", sup.col[i,2], ", ", sup.col[i,3], "]", sep = ""))
    pml.lines <- rbind(pml.lines, paste("label (///",groupId[nbgroup],"/", j, "), ", "\"", sup.names[i], "\"", sep = ""))
    pml.lines <- rbind(pml.lines, paste("color sup_", groupId[nbgroup],j, ", ///",groupId[nbgroup],"/", j, sep = ""))
if(ncol!=1){
    pml.lines <- rbind(pml.lines, paste("select sup_",project$col[i,2],"_",sup.names[i],",(///",groupId[nbgroup],"/",j,")",sep=""))
    resInGroup[[project$col[i,2]]]<-c(resInGroup[[project$col[i,2]]],j)
}
else {
    pml.lines <- rbind(pml.lines, paste("select sup_",project$col[1,2],"_",sup.names[i],",(///",groupId[nbgroup],"/",j,")",sep=""))
    resInGroup[[project$col[1,2]]]<-c(resInGroup[[project$col[1,2]]],j)
}
}
        j <- j + 1
    }
if(!is.null(file.pml)){
  pml.header <- rbind(pml.header,paste("select sup, resi 5001-",j-1,sep=""))
}
  }

  m <- (max(x$coord[, axis[1]], x$coord[, axis[2]], x$coord[, axis[3]]) * 2)
  if (!is.null(project) && inherits(project,"project")) {
    m <- (max(x$coord[, axis[1]], x$coord[, axis[2]], x$coord[, axis[3]], 
         project$coord[, axis[1]], project$coord[, axis[2]], project$coord[, axis[3]]) * 2)
  }

  all.nb <- active.nb
  if (!is.null(project) && inherits(project,"project"))
    all.nb <- (active.nb + sup.nb)

  pdb.lines <- rbind(pdb.lines, pdb.axis(all.nb + 1, m, 0, 0, "O"))
  pdb.lines <- rbind(pdb.lines, pdb.axis(all.nb + 2, 0, m, 0, "C"))
  pdb.lines <- rbind(pdb.lines, pdb.axis(all.nb + 3, 0, 0, m, "N"))
  pdb.lines <- rbind(pdb.lines, pdb.axis(all.nb + 4, 0, 0, 0, "Po"))

  pdb.lines <- rbind(pdb.lines, pdb.conect(all.nb + 1, all.nb + 4))
  pdb.lines <- rbind(pdb.lines, pdb.conect(all.nb + 2, all.nb + 4))
  pdb.lines <- rbind(pdb.lines, pdb.conect(all.nb + 3, all.nb + 4))
  if(!is.null(file.pml)){
   nameTmp<-sort(names(resInGroup))
   for(i in 1:length(nameTmp)){
	if(length(resInGroup[[nameTmp[i]]])>0){
		vec<-list()
		vec[[1]]<-resInGroup[[nameTmp[i]]][1]
		if(length(resInGroup[[nameTmp[i]]])>1){
			for(j in 2:length(resInGroup[[nameTmp[i]]])){
				if(resInGroup[[nameTmp[i]]][j] != (1+resInGroup[[nameTmp[i]]][j-1])){
					vec[[length(vec)]]<-paste(vec[[length(vec)]],resInGroup[[nameTmp[i]]][j-1],sep="-")
					vec[[length(vec) + 1]] <-resInGroup[[nameTmp[i]]][j]
				}
			}
			vec[[length(vec)]]<-paste(vec[[length(vec)]],resInGroup[[nameTmp[i]]][length(resInGroup[[nameTmp[i]]])],sep="-")
		}
		vec<-unlist(vec)
		pml.header <- rbind(pml.header, paste("select ",nameTmp[i],", (resi ",paste(vec,collapse=","),")",sep=""))
	}
}
pml.lines<- rbind(pml.header,pml.lines)
 pml.lines <- rbind(pml.lines, "hide labels")
  pml.lines <- rbind(pml.lines, "deselect")
  pml.lines <- rbind(pml.lines, "as spheres, ref")
if (!is.null(project) && inherits(project,"project")) {
  pml.lines <- rbind(pml.lines, "as nonbonded, sup")
}
  pml.lines <- rbind(pml.lines, "set sphere_scale, 5")
  pml.lines <- rbind(pml.lines, "set nonbonded_size, 10")
  cat(pml.lines, file = file.pml, sep = "\n")
}
  cat(pdb.lines, file = file.pdb, sep = "\n")
}
