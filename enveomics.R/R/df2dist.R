
enve.df2dist <- function(
	### Transform a dataframe (or coercible object, like a table) into a `dist` object.
	x,
	### A table (or coercible object) with at least three columns: (1) ID of the object 1,
	### (2) ID of the object 2, and (3) distance between the two objects.
	obj1.index=1,
	### Index of the column containing the ID of the object 1.
	obj2.index=2,
	### Index of the column containing the ID of the object 2.
	dist.index=3,
	### Index of the column containing the distance.
	default.d=NA,
	### Default value (for missing values)
	max.sim=0
	### If not-zero, assumes that the values are similarity (not distance)
	### and this is the maximum similarity (corresponding to distance 0).
	### Applies transformation: distance = (max.sim - values)/max.sim.
	){
   x <- as.data.frame(x);
   a <- as.character(x[, obj1.index]);
   b <- as.character(x[, obj2.index]);
   d <- as.double(x[, dist.index]);
   if(max.sim!=0) d <- (max.sim - d)/max.sim
   ids <- unique(c(a,b));
   m <- matrix(default.d, nrow=length(ids), ncol=length(ids), dimnames=list(ids, ids));
   diag(m) <- 0.0
   for(i in 1:nrow(x)){
      m[a[i], b[i]] <- d[i];
      m[b[i], a[i]] <- d[i];
   }
   return(as.dist(m));
   ### Returns a `dist` object.
}



enve.df2dist.group <- function(
	### Transform a dataframe (or coercible object, like a table) into a `dist` object, where
	### there are 1 or more distances between each pair of objects.
	x,
	### A dataframe (or coercible object) with at least three columns: (1) ID of the object 1,
	### (2) ID of the object 2, and (3) distance between the two objects.
	obj1.index=1,
	### Index of the column containing the ID of the object 1.
	obj2.index=2,
	### Index of the column containing the ID of the object 2.
	dist.index=3,
	### Index of the column containing the distance.
	summary=median,
	### Function summarizing the different distances between the two objects.
	empty.rm=TRUE
	### Remove rows with empty or NA groups
	){
   x <- as.data.frame(x);
   if(empty.rm) x <- x[ !(is.na(x[,obj1.index]) | is.na(x[,obj2.index]) | x[,obj1.index]=='' | x[,obj2.index]==''), ]
   a <- as.character(x[, obj1.index]);
   b <- as.character(x[, obj2.index]);
   d <- as.double(x[, dist.index]);
   ids <- unique(c(a,b));
   if(length(ids)<2) return(NA);
   m <- matrix(NA, nrow=length(ids), ncol=length(ids), dimnames=list(ids, ids));
   diag(m) <- 0
   for(i in 2:length(ids)){
      id.i <- ids[i];
      for(j in 1:(i-1)){
	 id.j <- ids[j];
	 d.ij <- summary(c( d[ a==id.i & b==id.j], d[ b==id.i & a==id.j] ));
	 m[id.i, id.j] <- d.ij;
	 m[id.j, id.i] <- d.ij;
      }
   }
   return(as.dist(m));
   ### Returns a `dist` object.
}

enve.df2dist.list <- function(
	### Transform a dataframe (or coercible object, like a table) into a `dist` object.
	x,
	### A dataframe (or coercible object) with at least three columns: (1) ID of the object 1,
	### (2) ID of the object 2, and (3) distance between the two objects.
	groups,
	### Named array where the IDs correspond to the object IDs, and the values correspond to
	### the group.
	obj1.index=1,
	### Index of the column containing the ID of the object 1.
	obj2.index=2,
	### Index of the column containing the ID of the object 2.
	dist.index=3,
	### Index of the column containing the distance.
	empty.rm=TRUE,
	### Remove incomplete matrices
	...
	### Any other parameters supported by `enve.df2dist.group`.
	){
   x <- as.data.frame(x);
   a <- as.character(x[, obj1.index]);
   b <- as.character(x[, obj2.index]);
   d <- as.numeric(x[, dist.index]);
   ids.all <- unique(c(a,b));
   l <- list();
   same_group <- groups[a]==groups[b];
   same_group <- ifelse(is.na(same_group), FALSE, TRUE);
   for(group in unique(groups)){
      ids <- ids.all[ groups[ids.all]==group ];
      if(length(ids)>1 & group!=""){
	 x.sub <- x[ same_group & (groups[a]==group) & (groups[b]==group), ]
	 if(nrow(x.sub)>0){
	    d.g <- enve.df2dist(x.sub, obj1.index, obj2.index, dist.index, ...);
	    if(!empty.rm | !any(is.na(d.g))) l[[ group ]] <- d.g;
	 }
      }
   }
   return(l);
   ### Returns a `list` of `dist` object.
}

