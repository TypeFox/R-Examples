sf2nms<-function (event.key, sform){
    sform.vec <- unlist(strsplit(sform, ""))
    sform.vec.val<-sform.vec[which(sapply(sform.vec,function(x) !grepl("[[:punct:]]|[[:digit:]]",x)))]
    str.rep <- event.key[sapply(sform.vec.val, function(x) which(x==event.key[, 1])), 2]
    sfm.rep <- which(sform.vec %in% event.key[, 1])
    sform.vec[sfm.rep] <- str.rep
    paste(sform.vec, collapse = "")
}
