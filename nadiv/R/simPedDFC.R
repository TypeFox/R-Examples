simPedDFC <- function(F, gpn = 4, fsn = 4, s = 2, prefix = NULL)
{

if(gpn < 2) stop("Number of founding grand-parents ('gpn') must be greater than or equal to 2")
if(fsn < 4) stop("Full-sib family size ('fsn') must be greater than or equal to 4")
if(s < 2) stop("Number of sires per full-sib family in P generation ('s') must be greater than or equal to 2")
if(floor(fsn/2) != (fsn/2)) stop("Full-sib family size ('fsn') must be an even number")

unitFun <- function(Fx){
   design <- matrix(NA, nrow = fsn, ncol = gpn)
   rc <- cbind(c(1:(floor(fsn/s)*s)), rep(1:floor(fsn/s), each = s))
   sires <- paste(rep(paste("u", Fx, "_s", seq.int(fsn/s), sep = ""), each = s), letters[1:(floor(fsn/2)*s)], sep = "")

   for(x in 1:dim(rc)[1]){
      design[x, rc[x,2]] <- sires[x]
      damNumb <- which(is.na(design[x, ]))
      design[x, damNumb] <- paste("u", Fx, "_d", damNumb, letters[x], sep = "")
   }
   sexDesign <- grepl(paste("^u", Fx, "_s", sep = ""), design)

   if(is.null(prefix)){
     unitPed <- data.frame(id = as.character(c(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), paste("u", Fx, "_gd", seq.int(gpn), sep = ""), c(design), paste("u", Fx, "_m", rep(seq.int((fsn*(gpn-1))), each = fsn), rep(c("m", "f"), each = (fsn/2)), rep(seq.int(fsn/2), (fsn*(gpn-1))), sep = ""))),
	dam = as.character(c(rep(NA, (2*gpn)), rep(paste("u", Fx, "_gd", seq.int(gpn), sep = ""), each = fsn), rep(unlist(sapply(t(design), FUN = function(x) {x[grepl(paste("^u", Fx, "_d", sep = ""), x)]})), each = fsn))),
	sire = as.character(c(rep(NA, (2*gpn)), rep(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), each = fsn), rep(sires, each = ((gpn-1)*fsn)))),
	sex = c(rep("M", gpn), rep("F", gpn), sapply(sexDesign, FUN = function(x) {if(x) "M" else "F"}), rep(rep(c("M", "F"), each = (fsn/2)), ((gpn-1)*fsn))))
   } else{
       p <- as.character(prefix)   
       unitPed <- data.frame(id = as.character(paste0(p, c(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), paste("u", Fx, "_gd", seq.int(gpn), sep = ""), c(design), paste("u", Fx, "_m", rep(seq.int((fsn*(gpn-1))), each = fsn), rep(c("m", "f"), each = (fsn/2)), rep(seq.int(fsn/2), (fsn*(gpn-1))), sep = "")))),
	dam = as.character(c(rep(NA, (2*gpn)), paste0(p, c(rep(paste("u", Fx, "_gd", seq.int(gpn), sep = ""), each = fsn), rep(unlist(sapply(t(design), FUN = function(x) {x[grepl(paste("^u", Fx, "_d", sep = ""), x)]})), each = fsn))))),
	sire = as.character(c(rep(NA, (2*gpn)), paste0(p, c(rep(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), each = fsn), rep(sires, each = ((gpn-1)*fsn)))))),
	sex = c(rep("M", gpn), rep("F", gpn), sapply(sexDesign, FUN = function(x) {if(x) "M" else "F"}), rep(rep(c("M", "F"), each = (fsn/2)), ((gpn-1)*fsn))))
     }


   unitPed
   }


 ped_out <- do.call(rbind, lapply(seq.int(F), FUN = unitFun))


return(ped_out)
}

