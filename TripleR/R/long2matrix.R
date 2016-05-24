# function takes data in the "long" format and returns a list with quadratic
# round robin matrices for each group
#' @export
long2matrix <- function(formule, data, minData=1, verbose=TRUE, reduce=TRUE, skip3=TRUE, g.id=NULL, exclude.ids=NULL, ...) {
	
	extra <- list(...)

	# parse formula
	#remove spaces from formula
	lhs <- strsplit(gsub(" ","",as.character(formule)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(formule)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]
	
	var.id <- lhs
	actor.id <- rhs[1]
	partner.id <- rhs[2]
	group.id <- NULL
	if (length(rhs)>=3) {group.id <- rhs[3]}
	
	# What are the group ids?
	# if only one group, add group variable
	if (!is.null(group.id)) {if (is.na(group.id)) {group.id <- NULL}}
	if (!is.null(group.id)) {
			# Check if IDs are unique across groups
			t1 <- table(data[, group.id], data[, actor.id])
			unique <- apply(t1, 2, function(col) {
				sum(col != 0)
			})
			if (any(unique > 1)) stop("Some of your actor or partner IDs are not unique (i.e., different persons in different groups have the same ID). Please make them unique. Calculations are aborted.")
			
			group.ids <- names(table(data[,group.id]))
	} else {
		if (!is.null(g.id)) {
			group.id <- g.id
			group.ids <- levels(factor(data[,group.id]))
		} else {
			group.ids <- "1"
			group.id <- "group.id"
			data$group.id <- "1"
		}
	}
	
	# reduce data.frame to relevant variables
	data <- data[, colnames(data) %in% c(actor.id, partner.id, var.id, group.id)]
	
	res <- list()
	del.IDs <- c()	# the complete list of deleted participants
	del.groups <- c()	# the complete list of deleted groups
	
	for (g in group.ids) {
		
		#print(paste("Processing group",g))
		block <- data[data[,group.id]==g,]
		# block <- block[!is.na(block[,var.id]),]
		
		# reduce ids to factor levels which are actually present
		block[,actor.id] <- factor(as.character(block[,actor.id]), ordered=TRUE)
		block[,partner.id] <- factor(as.character(block[,partner.id]), ordered=TRUE)
		block[,group.id] <- factor(as.character(block[,group.id]), ordered=TRUE)
						
		# get names of participants which served both as actors and observers
		if (reduce==TRUE) {
			
			# Schnittmenge aus actors und partners herausfinden
			p <- intersect(levels(as.factor(block[,actor.id])), levels(as.factor(block[,partner.id])))
			block.clean <- block[(block[,actor.id] %in% p) & (block[,partner.id] %in% p), which(colnames(block) != group.id)]
			
			# reduce levels again
			block.clean[,actor.id] <- factor(as.character(block.clean[,actor.id]), ordered=TRUE)
			block.clean[,partner.id] <- factor(as.character(block.clean[,partner.id]), ordered=TRUE)
			
			delComplete.IDs <- setdiff(union(levels(block[,actor.id]), levels(block[,partner.id])), p)
			p2 <- p
			
		} else {
			delComplete.IDs <- c()
			p2 <- union(levels(as.factor(block[,actor.id])), levels(as.factor(block[,partner.id])))
		}
			
			
		if (length(p2) <= 1) {
			if (verbose==TRUE){
				warning(paste("Warning: The provided data of group", g, "are not in round robin format!"), call.=FALSE);
			}
			next();
		}
			
		block.clean <- block[block[,actor.id] %in% p2 & block[,partner.id] %in% p2, which(colnames(block) != group.id)]

		
		# finally: construct the quadratic matrix
		f1 <- as.formula(paste(actor.id,"~",partner.id))
		
		# check for double entries
		if (any(table(block.clean[, c(actor.id, partner.id)]) > 1)) {
			stop(paste0("There are double entries in group ", g, ". Calculations aborted, please remove these double entries."))
		}

		
		#box <- as.matrix(cast(block.clean, f1, value=var.id, fill=NA))		
		box <- as.matrix(reshape2::acast(block.clean, f1, value.var=var.id, fill=NA))
		
		
		if (reduce==FALSE) {
			colmiss <- setdiff(levels(as.factor(block[,actor.id])), levels(as.factor(block[,partner.id])))
			rowmiss <- setdiff(levels(as.factor(block[,partner.id])), levels(as.factor(block[,actor.id])))
			if (length(colmiss)>0) box <- cbind(box, matrix(NA, nrow=nrow(box), ncol=length(colmiss), dimnames=list(NULL, colmiss)))
			if (length(rowmiss)>0) box <- rbind(box, matrix(NA, ncol=ncol(box), nrow=length(rowmiss), dimnames=list(rowmiss, NULL)))
			
			box <- box[order(rownames(box)), order(colnames(box))]
		}
		
		# extract self ratings if present
		self <- diag(box)
		diag(box) <- NA
		
		# Voyeure entfernen (nur actors oder nur partners) - iterativ
		delVoyeur.IDs <- c()
		if (reduce) {
			repeat {
				del.row <- apply(box, 1, function(x) return((length(x) - sum(is.na(x))) < minData))
				del.col <- apply(box, 2, function(x) (length(x) - sum(is.na(x))) < minData)
				dels <- del.row | del.col
				if (sum(dels==TRUE) == 0) {break();}
				delVoyeur.IDs <- c(delVoyeur.IDs, colnames(box)[dels])
				box <- box[!dels,!dels]
				self <- self[!dels]
				if (is.null(box) | length(box) <= 1) break();
			}
		}
		
		if (length(box) <= 1) box <- NULL
		
		# Exclude the exclude.ids
		ex2 <- intersect(rownames(box), exclude.ids)
		self <- self[!rownames(box) %in% ex2]
		if (!is.null(box)) box <- box[!rownames(box) %in% ex2,!colnames(box) %in% ex2]
		
		## Warnung ausgeben
		del.IDs <- c(del.IDs, union(delComplete.IDs, union(delVoyeur.IDs, ex2)))
		del2 <- union(delComplete.IDs, union(delVoyeur.IDs, ex2))
		if (verbose==TRUE & length(del2) > 0) {
			
			if (is.null(extra[["bistyle"]])) {
				warning(paste(var.id,": ",length(del2)," participant(s) have been excluded from group",g,"due to exceedingly missing data; id(s) =",paste(del2, collapse=", "),"."), call.=FALSE)
			} else {
				warning(paste(length(del2)," participant(s) have been excluded from group",g,"due to exceedingly missing data in at least one variable; id(s) =",paste(del2, collapse=", "),"."), call.=FALSE)
			}
			
		}
		
		# | length(box)==1
		if (is.null(box)) {
			del.groups <- c(del.groups, g)
			if (verbose==TRUE) {
				if (is.null(extra[["bistyle"]])) {
					warning(paste(var.id,": group",g,"has 3 or fewer subjects - the group is excluded from the analyses."), call.=FALSE)
				} else {
					warning(paste("Group",g,"has 3 or fewer subjects - the group is excluded from the analyses."), call.=FALSE)
				}
			}
		} else {
			if (!skip3 | nrow(box)>3) {
				res[[g]] <- box
				attr(res[[g]], "group.id") <- g
				attr(res[[g]], "varname") <- var.id
				if (any(!is.na(self))) {
					attr(res[[g]], "self.ratings") <- self
				}
			}
			if (skip3 & nrow(box)<=3) {
				del.groups <- c(del.groups, g)
				if (verbose==TRUE) {
					if (is.null(extra[["bistyle"]])) {
						warning(paste(var.id,": group",g,"has 3 or fewer subjects - the group is excluded from the analyses."), call.=FALSE)
					} else {
						warning(paste("Group",g,"has 3 or fewer subjects - the group is excluded from the analyses."), call.=FALSE)
					}
				}
			}
		}
	
	}
	
	# After processing all groups: attach excluded participant IDs
	if (length(del.IDs)>0) {attr(res, "excluded.participants") <- sort(del.IDs)} else {attr(res, "excluded.participants") <- NULL}
	if (length(del.groups)>0) {attr(res, "excluded.groups") <- sort(del.groups)} else {attr(res, "excluded.groups") <- NULL}
	
	
	if (length(res)==0) {return()} 
	else {return(res)}
}



# function takes data in matrix format and turns them into long format
#' @export
matrix2long <- function(M, new.ids=TRUE, var.id="value") {
	M <- as.matrix(M)
	if (new.ids) {
		rownames(M) <- colnames(M) <- seq(1, nrow(M))
	}
	m1 <- reshape2::melt(as.matrix(M))
	colnames(m1) <- c("actor.id", "partner.id", var.id)
	
	return(m1)
}


