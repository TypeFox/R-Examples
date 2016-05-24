ccda.main <-
function (dataset, names_vector, nr, nameslist, prior = "proportions", 
    return.RCDP = FALSE) 
{
    length_dataset = length(dataset[, 1])
    dataset = cbind(names_vector, dataset)
    st = 2
    dim_dataset = length(dataset[1, ])
    ngroups = sum(table(dataset[, 1]) != 0)
    if (length(names_vector) != length_dataset) {
        stop("length of names vector does not equal to the length of dataset")
    }
    if (length(nameslist) != ngroups) {
        stop("number of names does not equal the number of various names in names vector")
    }
    if (ngroups <= 1) {
        stop("all data have the same origin")
    }
	if(prior!="proportions"&prior!="equal"){
		stop("prior has to be set to proportions or equal")
	}
	if((R.Version()$major>=3) & (R.Version()$minor>=1.0)){
		mth="ward.D"
	}
	else{
		mth="ward"
	}
	
    mean_mtx = matrix(rep(0, ngroups * (-st + dim_dataset + 1)), 
        nrow = ngroups)
    for (i in 1:ngroups) {
        mean_mtx[i, ] = colMeans(dataset[(dataset[, 1] == nameslist[i]), 
            st:dim_dataset])
    }
    cluster = hclust(dist(scale(mean_mtx)) * dist(scale(mean_mtx)), 
        method = mth)
    clusterings = cutree(cluster, k = c(1:ngroups))
    grouping_mtx = matrix(rep(1, times = ngroups * length_dataset), 
        ncol = ngroups)
    for (x in 1:ngroups) {
        for (z in 1:ngroups) {
            for (i in 1:length_dataset) {
                if (dataset[i, 1] == nameslist[z]) {
                  grouping_mtx[i, x] = clusterings[z, x]
                }
            }
        }
    }
    grouping_percentage = c(rep(1, times = ngroups))
    for (x in 2:ngroups) {
        grouping_percentage[x] = percentage(dataset[, st:dim_dataset], 
            grouping_mtx[, x], prior)
    }
    random_grouping_percentage = matrix(rep(1, times = ngroups * 
        nr), ncol = nr)
    for (y in 1:nr) {
        random_grouping_mtx = grouping_mtx[sample(1:length_dataset), 
            ]
        for (x in 2:ngroups) {
            random_grouping_percentage[x, y] = percentage(dataset[, 
                st:dim_dataset], random_grouping_mtx[, x], prior)
        }
    }
    q95 = rep(1, ngroups)
    for (i in 1:ngroups) {
        q95[i] = quantile(random_grouping_percentage[i, ], prob = 0.95)
    }
    ratio = grouping_percentage
    D = ratio - q95
    k = min(which(D == max(D)))
    groups = matrix(paste("sub-group", clusterings[, k]), ncol = 1, 
        nrow = ngroups)
    rownames(groups) = nameslist
    cat(paste("\n", "Number of optimal groups: ", k, "\n", sep = ""))
    cat(paste("\n", "Maximal difference between q95 and ratio ", 
        round(max(D), digits = 3), "\n", "\n", sep = ""))
    cat(paste("further investigation of the following ", k, "sub-groups  recommended:"))
    print(groups)
    groups = matrix(paste("sub-group", clusterings[, k]), ncol = 1, 
        nrow = ngroups, dimnames = list(nameslist, paste("further investigation of the following ", 
            k, "sub-groups  recommended:")))
    if (return.RCDP == TRUE) {
        invisible(list(nameslist = nameslist, q95 = q95, ratio = ratio, 
            difference = D, sub_groups = groups, cluster = cluster, 
            RCDP = random_grouping_percentage))
    }
    else {
        invisible(list(nameslist = nameslist, q95 = q95, ratio = ratio, 
            difference = D, sub_groups = groups, cluster = cluster))
    }
}
