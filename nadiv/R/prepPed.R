prepPed <- function(pedigree, gender = NULL, check = TRUE){

 if(check){      
   if(length(which(pedigree[, 2] == 0)) > 0){
     pedigree[which(pedigree[, 2] == 0), 2] <- NA
     warning("Zero in the dam column interpreted as a missing parent")
   }
   if(length(which(pedigree[, 3] == 0)) > 0){
     pedigree[which(pedigree[, 3] == 0), 3] <- NA
     warning("Zero in the sire column interpreted as a missing parent")
   }
   if(length(which(pedigree[,2] == "*")) > 0) pedigree[which(pedigree[, 2] == "*"), 2] <- NA
   if(length(which(pedigree[,3] == "*")) > 0) pedigree[which(pedigree[, 3] == "*"), 3] <- NA
   if(all(is.na(pedigree[, 2])) & all(is.na(pedigree[, 3]))){
     stop("All dams and sires are missing")
   }
   if(any(pedigree[, 1] == 0 | pedigree[, 1] == "0" | pedigree[, 1] == "*" | is.na(pedigree[, 1]))){
     warning("Missing value in the ID column - row discarded")
     warning("Check to ensure first three columns of the pedigree object are ID, Dam, and Sire")
     pedigree <- pedigree[-which(pedigree[, 1] == 0 | pedigree[, 1] == "0" | pedigree[, 1] == "*" | is.na(pedigree[, 1])), ]
   }
 }

  
 facflag <- is.factor(pedigree[, 1])
 udam <- if(facflag) as.character(unique(pedigree[, 2])) else unique(pedigree[, 2])
 udam <- udam[!is.na(udam)]
 missdam <- udam[which(is.na(match(udam, pedigree[, 1])))]
 usire <- if(facflag) as.character(unique(pedigree[, 3])) else unique(pedigree[, 3])
 usire <- usire[!is.na(usire)]
 misssire <- usire[which(is.na(match(usire, pedigree[, 1])))]

 if(length(missdam) == 0 & length(misssire) == 0){
   ped_fixed <- pedigree
 } else{
     topPed <- data.frame(c(missdam, misssire), rep(NA, length(missdam) + length(misssire)), rep(NA, length(missdam) + length(misssire)), matrix(NA, nrow = (length(missdam) + length(misssire)), ncol = ncol(pedigree) - 3))
  names(topPed) <- names(pedigree)
     if(!is.null(gender)){
        if(is.factor(pedigree[, gender])){
          damgender <- as.character(pedigree[which(pedigree[, 1] == udam[which(!udam %in% missdam)][1]), gender]) 
          siregender <- as.character(pedigree[which(pedigree[, 1] == usire[which(!usire %in% misssire)][1]), gender]) 
        } else{
             damgender <- pedigree[which(pedigree[, 1] == udam[which(!udam %in% missdam)][1]), gender] 
             siregender <- pedigree[which(pedigree[, 1] == usire[which(!usire %in% misssire)][1]), gender]
          } 
        topPed[, gender] <- c(rep(damgender, length(missdam)), rep(siregender, length(misssire)))
     }
     ped_fixed <- rbind(topPed, pedigree)
   }
 npf <- nrow(ped_fixed)

 if(sum((na.omit(ped_fixed[, 2]) %in% ped_fixed[, 1]) == FALSE) > 0 & any(is.na(ped_fixed[, 2]) == FALSE)){
   stop("Something wicked happened: individuals appearing as dams but not added to pedigree. Report as possible bug to <matthewwolak@gmail.com>")
 }
 if(sum((na.omit(ped_fixed[, 3]) %in% ped_fixed[, 1]) == FALSE) > 0 & any(is.na(ped_fixed[, 3]) == FALSE)){
   stop("Something wicked happened: individuals appearing as sires but not added to pedigree. Report as possible bug to <matthewwolak@gmail.com>")
 }
 if(sum(duplicated(ped_fixed[, 1])) > 0){
   stop("some individuals appear more than once in the pedigree")
 }

 nPed_fixed <- numPed(ped_fixed[, 1:3], check = FALSE)
 Cout <- .C("gaUnsort",
	as.integer(nPed_fixed[, 2] - 1),
	as.integer(nPed_fixed[, 3] - 1),
        as.integer(rep(0, npf)),
	as.integer(rep(0, npf)),
	as.integer(npf))

 ped_fixed_ord <- ped_fixed[order(Cout[[3]], Cout[[4]]), ]
 itwork <- try(expr = numPed(ped_fixed_ord[, 1:3]), silent = TRUE)
 if(class(itwork) == "try-error"){
   G <- Matrix(FALSE, npf, npf, sparse = TRUE)
   G[cBind(c(nPed_fixed[which(nPed_fixed[, 2] != -998), 2], nPed_fixed[which(nPed_fixed[, 3] != -998), 3]), c(nPed_fixed[which(nPed_fixed[, 2] != -998), 1], nPed_fixed[which(nPed_fixed[, 3] != -998), 1]))] <- TRUE
   Gtmp <- G
   gconv <- Matrix(TRUE, nrow = 1, ncol = npf, sparse = TRUE)
   gendepth <- rep(0, npf) + as((gconv %*% Gtmp), "ngCMatrix") 
   while(nnzero(Gtmp) > 0){
     Gtmp <- Gtmp %*% G
     gendepth <- gendepth + as((gconv %*% Gtmp), "ngCMatrix") 
   }
   ped_fixed_ord <- ped_fixed[order(as(gendepth, "matrix")), ]
 }

 return(ped_fixed_ord) 
}

