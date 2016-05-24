 classify.qtl<-function (cross, peak, etrait.coord, data.gmap) 
{
    require(qtl)
    if (length(class(cross)) < 2 || class(cross)[2] != "cross") 
        stop("Input should have class \"cross\".")
    if (!all(attr(peak, "class", exact = TRUE) %in% c("peak", "list"))) 
        stop("Input should have class \"peak\".")
    scanone <- get(attr(peak, "scanone", exact = TRUE))
    if (!all(levels(scanone$chr) == names(cross$geno))) 
        stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which was performed on cross object")
    if (!all(scanone$pos == pseudo.map(cross))) 
        stop("Arguments peak and/or cross misspecified: attributes(peak)$scanone should describe a scanone object which was performed on cross object")
    if (class(data.gmap)[1] != "data.frame" || any(!(names(data.gmap) %in% c("Marker", "chr", "PP")))) 
        stop("data.gmap should have class \"data.frame\" with columns names: 'Marker','chr','PP'")
    if (missing(data.gmap)) 
        stop("Argument data.gmap unspecified")
    if (!all(as.vector(data.gmap$Marker) %in% mnames.map(cross))) 
        stop("Arguments cross and/or data.gmap misspecified: data.gmap should have the same marker names as cross")
    if (class(etrait.coord)[1] != "data.frame" || !all(names(etrait.coord) %in% c("etrait.name", "chr", "start", "stop"))) 
        stop("etrait.coord should have class \"data.frame\" and columns names: 'etrait.name','chr','start','stop' \n")
    if ("type" %in% attr(peak, "features", exact = TRUE)) {
        cat("\nWARNING: Feature 'type' is already exits. It will be replace by the new one.\n\n")
        peak <- drop.peakfeat(peak, "type")
    }
    res <- list(un = NA)
    for (i in 1:length(peak)) {
        trait <- names(peak[i])
        resbytrait <- list(un = NA)
        bool <- toupper(etrait.coord$etrait.name) %in% toupper(trait)
        if (any(bool)) {
            start <- etrait.coord$start[bool]
            stop <- etrait.coord$stop[bool]
        }
        else {
            cat("trait:", trait, "\tSTATUT: not found in etrait.coord\n")
            start <- NA
            stop <- NA
        }
        for (y in 1:length(peak[[i]])) {
            if (!is.na(peak[[i]][y][1])) {
                chr <- as.numeric(names(peak[[i]][y]))
                ty <- NA
                for (z in seq(peak[[i]][[y]]$mname.peak)) {
                  if (is.na(start)) {
                    type <- data.frame(type = NA)
                    type <- cbind(peak[[i]][[y]], type)
                    type <- list(type)
                  }
                  else {
                    if (names(peak[[i]][y]) != paste(etrait.coord$chr[bool])) 
                      ty <- c(ty, paste("trans"))
                    else {
                      inf <- peak[[i]][[y]]$inf.cM[z]
                      sup <- peak[[i]][[y]]$sup.cM[z]
                      ic <- c(inf, sup)

                      for (w in 1:2) {
                        mflank <- find.flanking(cross, chr, ic[w])
                        posg_left <- find.markerpos(cross, mflank$left)$pos
                        posg_right <- find.markerpos(cross, mflank$right)$pos
                        bool2 <- data.gmap$Marker %in% paste(mflank$left) 
                        posp_left <- data.gmap$PP[bool2]
                        bool2 <- data.gmap$Marker %in% paste(mflank$right)
                        posp_right <- data.gmap$PP[bool2]
                        if (posg_left == posg_right) {
                          if( posg_left == 0 ){
                            ratio <- 1
                          } else {
                            ratio <- posg_left / posg_right
                          }
                        } else {
                          ratio <- (posg_left - posg_right)/(posp_left - posp_right)
                        }
                        #if(ic[w] == 0) ic[w] <- ic[w]+1
                        ic[w] <- (ic[w] - posg_left) / ratio + posp_left                 
                      }
                   
                      if ((ic[1] <= start && start <= ic[2]) || (ic[1] <= stop && stop <= ic[2])) {
                        ty <- c(ty, paste("cis"))
                      }
                      else ty <- c(ty, paste("trans"))
                    }
                  }
                }
                if (!is.na(start)) {
                  type <- data.frame(type = ty[-1])
                  type <- cbind(peak[[i]][[y]], type)
                  type <- list(type)
                }
            }
            else {
                type <- NA
            }
            attributes(type)$names <- paste(names(peak[[i]][y]))
            resbytrait <- c(resbytrait, type)
        }
        resbytrait <- list(resbytrait[-1])
        attributes(resbytrait)$names <- trait
        res <- c(res, resbytrait)
    }
    res <- res[-1]
    attributes(res)$class <- c("peak", "list")
    attributes(res)$features <- c(attributes(peak)$features, 
        "type")
    attributes(res)$scanone <- attributes(peak)$scanone
    attributes(res)$lod.th <- attributes(peak)$lod.th
    attributes(res)$si <- attributes(peak)$si
    attributes(res)$window <- attributes(peak)$window
    return(res)
}
