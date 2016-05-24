readratdatafixed <- function(filename, pnames = c(0, 0), pschoolid = c(1, 
    2), pgroupid = c(3, 4), prespid = c(5, 6), pratings = c(8, 17), rowsEQassessors = T) {
  
    readgroupdata <- function(group_data, groupnr = 1) {
        nrresp = length(group_data)
        if (!all(pnames == 0)) 
            {
              resplabel = substr(group_data, pnames[1], pnames[2])
              resplabel = sub("\\s+$", "", resplabel)
            }  # delete trailing whitespace
        if (!all(pschoolid == 0)) 
            schoolid = as.numeric(substr(group_data, pschoolid[1], pschoolid[2]))
        else schoolid=1
        if (!all(pgroupid == 0)) 
            groupid = as.numeric(substr(group_data, pgroupid[1], pgroupid[2]))
        if (!all(prespid == 0)) 
            respid = as.numeric(substr(group_data, prespid[1], prespid[2]))
        if (!all(pratings == 0)) {
            lr = min(pratings[2], nchar(group_data[1]))
            if (lr <= 2) stop(paste("Unexpected empty line in ", groupnr, "th group.", sep=''))
            ratings = substring(group_data, pratings[1], lr)
            nrratings = lr - pratings[1] + 1
            ratingnames=sprintf("r%02.0f", 1:nrratings)
            r = unlist(strsplit(ratings, split = NULL))
            r[r == "0"] = NA
            if (rowsEQassessors) {
                rmat = t(suppressWarnings(matrix(as.numeric(r), nrow = nrratings, 
                  ncol = nrresp)))
            } else {
                rmat = suppressWarnings(matrix(as.numeric(r), nrow = nrratings, 
                  ncol = nrresp))
            }
        }
        
        DF = as.data.frame(cbind(schoolid, groupid, respid, rmat))
        colnames(DF) <- c("schoolid", "groupid", "respid", ratingnames)
        if (all(pnames > 0)) {
            DF = cbind(resplabel, DF)
        }
        
        return(DF)
    }
    
    
    all_data = readLines(filename)
    nlines = length(all_data)
    
    emptylines = which(sapply(all_data, nchar) == 0)
    nrc = length(emptylines)
    
    start = 1
    nrgroup = 0
    if (nrc == 0) {
        group_data = all_data
        DF = readgroupdata(group_data, 1)
    } else {
        # i=1
        for (i in 1:nrc) {
            group_data = all_data[start:(emptylines[i] - 1)]
            nrgroup = nrgroup + 1
            if (nrgroup == 1) {
                DF = readgroupdata(group_data, nrgroup)
            } else {
                temp = readgroupdata(group_data, nrgroup)
                DF = rbind.fill(DF, temp)
            }
            start = emptylines[i] + 1
        }
        nrgroup = nrgroup + 1
        group_data = all_data[start:nlines]
        temp = readgroupdata(group_data, nrgroup)
        DF = rbind.fill(DF, temp)
    }
    return(DF)
}






