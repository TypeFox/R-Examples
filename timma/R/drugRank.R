#' Generate the list of ranked drug combinations
#' 
#' A function to provide a list of drug combinations ranked by their synergy scores
#' 
#' @param profile_select the selected drug-target interaction data
#' @param predicted_matrix the predicted efficacy matrix
#' @param sens the drug sensitivity vector.
#' 
#' @return a matrix contains the information about the list of drug combinations ranked by their synergy scores. 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' float<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[, 1], max_k = 8)
#' k_select<-float$k_sel
#' x<-data.frame(tyner_interaction_binary)
#' kinase_names <- dimnames(x)[[2]]
#' select_kinase_names <- findSameSet(x, k_select, kinase_names)
#' gc_timma <- graycode3(length(k_select))
#' gc_names <- graycodeNames(length(k_select), select_kinase_names, gc_timma$gc_row, gc_timma$gc_col)
#' nr <- gc_names$nr
#' nc <- t(gc_names$nc)
#' timma_row <- nrow(nr) + nrow(nc)
#' timma_col <- ncol(nr) + ncol(nc)
#' timma <- array("", dim = c(timma_row, timma_col))
#' timma[(nrow(nc) + 1):timma_row, 1:ncol(nr)] <- nr
#' timma[1:nrow(nc), (ncol(nr) + 1):timma_col] <- nc
#' timma[(nrow(nc) + 1):timma_row, (ncol(nr) + 1):timma_col] <- float$timma$dummy
#' profile_select<-data.frame(tyner_interaction_binary)[, k_select]
#' drug_combo_rank<-drugRank(profile_select, timma, tyner_sensitivity[, 1])
#' }

drugRank <- function(profile_select, predicted_matrix, sens) {
    profile_select <- cbind(dimnames(profile_select)[[1]], profile_select)
    ndrugs = dim(profile_select)[1]
    ntargets = dim(profile_select)[2] - 1
    
    timma <- predicted_matrix
    len1=(ntargets-floor(ntargets/2))
    len2=2^(floor(ntargets/2))
    len3=floor(ntargets/2)
    len4=2^(ntargets-floor(ntargets/2))
    
    if(len3==1) row_target = matrix(timma[c((1 + len1):(len2 + len1)),1],ncol=1) else row_target = apply(timma[c((1 + len1):(len2 + len1)), c(1:len3)], 2, as.character)
    if(len1==1) col_target = matrix(timma[1, c((1 + len3):(len4 + len3))],nrow=1) else col_target = apply(timma[c(1:len1), c((1 + len3):(len4 + len3))], 2, as.character)
    
    efficacy_mat = apply(timma[c((1 + len1):(len2 + len1)), c((1 + len3):(len4 + len3))], 2, as.numeric)
    # efficacy_mat = efficacy_mat/max(efficacy_mat)  # normalize the efficacy into range [0 1]
    efficacy_vec = matrix(efficacy_mat, length(c(efficacy_mat)), 1)
    
    # identical targets separated by '/'.  Replace ';' with '/' in timma1.csv Replace ';' with '/' in
    # selectedKinasse1.csv empty entries denoted by NA. Replace '-' with NA in timma1.csv
    row_target[which(row_target == "-")] = NA
    col_target[which(col_target == "-")] = NA
    row_target = gsub(";", "/", row_target)
    col_target = gsub(";", "/", col_target)
    
    names_mat = matrix(0,dim(row_target)[1], dim(col_target)[2])  # target combinations
    number_mat = matrix(0,dim(row_target)[1], dim(col_target)[2])  # number of target nodes in the combination
    
    for (i in 1:dim(row_target)[1]) {
        for (j in 1:dim(col_target)[2]) {
            str1 = paste(paste(row_target[i, ], collapse = " "), paste(col_target[, j], collapse = " "), collapse = "", 
                sep = " ")
            # str2 = gsub('NA','',str1)
            str2 = gsub("\\<NA\\>", "", str1)  # exact match
            str3 = gsub("^\\s+|\\s+$", "", str2)  # remove leading and trailing space
            str4 = strsplit(str3, "\\s+")
            names_mat[i, j] = paste(unlist(str4), collapse = ";")
            number_mat[i, j] = length(str4[[1]])
        }
    }
    names_vec = matrix(names_mat, length(c(names_mat)), 1)
    number.vec = matrix(number_mat, length(c(number_mat)), 1)
    data_target = data.frame(cbind(names_vec, efficacy_vec, number.vec))
    colnames(data_target) = c("Gene", "timma", "Number")
    data_target$timma = as.numeric(as.character(data_target$timma))
    
    # deal with only single and pairwise targets
    data_single = data_target[which(data_target$Number == 1), ]
    data_pairwise = data_target[which(data_target$Number == 2), ]
    
    
    single_inhibition = matrix(0,dim(data_pairwise)[1], 5)
    colnames(single_inhibition) = c("timma1", "timma2", "timma.add", "timma.multi", "timma.highest")
    for (i in 1:dim(data_pairwise)[1]) {
        pair = data_pairwise$Gene[i]
        pair_name = unlist(strsplit(as.character(pair), ";"))
        # index = lapply(pair_name, function(x) grep(x, data_single$Gene, fixed = TRUE))
        index = lapply(pair_name, function(x) which(data_single$Gene==x))
        timma1 = data_single$timma[unlist(index[1])]
        timma2 = data_single$timma[unlist(index[2])]
        timma.add = data_pairwise$timma[i] - timma1 - timma2
        timma.multi = data_pairwise$timma[i] - timma1 * timma2
        # timma.multi2 = timma.add+timma1*timma2
        timma.highest = data_pairwise$timma[i] - max(timma1, timma2)
        single_inhibition[i, ] = c(timma1, timma2, timma.add, timma.multi, timma.highest)
    }
    data_pairwise = cbind(data_pairwise, single_inhibition)
    
    # find the drugs that are corresponding to the targets
    names(profile_select)[1] = "Drugs"
    colnames(profile_select) = gsub(";", "/", colnames(profile_select), fixed = TRUE)
    table_drug = list(0)
    for (i in 1:ntargets) {
        target = colnames(profile_select)[i + 1]
        index = which(profile_select[, i + 1] == 1)
        drugs = profile_select$Drugs[index]
        # table_drug[[i]] = expand.grid(target,drugs)
        table_drug[[i]] = expand.grid(target, index)
    }
    table_drug = do.call("rbind", table_drug)
    
    
    index = table_drug[, 2]
    table_drug = cbind(table_drug, profile_select$Drugs[index], sens[index])
    colnames(table_drug) = c("Gene", "Index", "Drug", "Sens")
    
    score_drug = list(0)
    for (i in 1:dim(data_pairwise)[1]) {
        pair = data_pairwise$Gene[i]
        pair_str = strsplit(as.character(pair), ";")[[1]]
        gene1 = pair_str[1]
        gene2 = pair_str[2]
        # drug1_index = grep(gene1, table_drug$Gene, fixed = TRUE)  # indices of drugs that are targeting gene1
        drug1_index = which(table_drug$Gene==gene1)
        # drug2_index = grep(gene2, table_drug$Gene, fixed = TRUE)  # indices of drugs that are targeting gene2
        drug2_index = which(table_drug$Gene==gene2)
        # drugs that are targeting gene1
        drug1 = table_drug$Drug[drug1_index]
        # drugs that are targeting gene2
        drug2 = table_drug$Drug[drug2_index]
        n = expand.grid(drug1, drug2)
        n = t(apply(n, 1, sort))  # sort drug names alphabetically
        new_order = t(apply(expand.grid(drug1, drug2), 1, order))
        
        tmp = mat.or.vec(dim(n)[1], 8)
        for (j in 1:dim(n)[1]) {
            tmp[j, 1] = as.character(n[j, 1])
            tmp[j, 2] = as.character(n[j, 2])
            tmp[j, 3] = table_drug$Sens[which(table_drug$Drug == n[j, 1])[1]]
            tmp[j, 4] = table_drug$Sens[which(table_drug$Drug == n[j, 2])[1]]
            tmp[j, 5] = data_pairwise[i, 3 + new_order[j, 1]]
            tmp[j, 6] = data_pairwise[i, 3 + new_order[j, 2]]
            tmp[j, 7] = pair_str[new_order[j, 1]]
            tmp[j, 8] = pair_str[new_order[j, 2]]
        }
        
        tmp = cbind(tmp, data_pairwise[rep(i, each = dim(n)[1]), ])
        
        score_drug[[i]] = tmp
    }
    
    score_drug_mat = do.call("rbind", score_drug)
    colnames(score_drug_mat)[1:8] = c("Drug1", "Drug2", "Drug1.sens", "Drug2.sens", "timma1.s", "timma2.s", 
        "Target1", "Target2")
    score_drug_mat$Drug1.sens = as.numeric(as.character(score_drug_mat$Drug1.sens))
    score_drug_mat$Drug2.sens = as.numeric(as.character(score_drug_mat$Drug2.sens))
    score_drug_mat$Drug1 = as.character(score_drug_mat$Drug1)
    score_drug_mat$Drug2 = as.character(score_drug_mat$Drug2)
    score_drug_mat$timma1.s = as.numeric(as.character(score_drug_mat$timma1.s))  # sorted timma1
    score_drug_mat$timma2.s = as.numeric(as.character(score_drug_mat$timma2.s))  # sorted timma2
    
    drug_combinations = c(0)
    for (i in 1:dim(score_drug_mat)[1]) {
        drug_combinations[i] = paste(sort(c(score_drug_mat$Drug1[i], score_drug_mat$Drug2[i])), collapse = " ")
    }
    score_drug_mat$drug_combinations = as.factor(drug_combinations)
    comp1 = aggregate(score_drug_mat[, c(1:2)], list(drug_combinations), unique)
    comp2 = aggregate(score_drug_mat[, c(3:6, 10, 12:16)], list(drug_combinations), mean)
    comp3 = aggregate(score_drug_mat[, 7], list(drug_combinations), function(i) paste(unique(i), collapse = ","))
    comp4 = aggregate(score_drug_mat[, 8], list(drug_combinations), function(i) paste(unique(i), collapse = ","))
    colnames(comp3)[2] = "Target1"
    colnames(comp4)[2] = "Target2"
    score_drug_comb = cbind(comp1, comp2, comp3, comp4)
    score_drug_comb = score_drug_comb[,c(2,3,9,12,13,14,16,18)]
    colnames(score_drug_comb) = c("Drug1", "Drug2", "Sensitivity", "Synergy.add", "Synergy.multi", "Synergy.highest", 
    "Target1", "Target2")
    return(score_drug_comb)
}
 
