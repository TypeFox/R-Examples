toLatexTable = function(ANOVA, EF, fixed.names) {
    
    matchRowNames = rownames(ANOVA)
    random.ColNames = colnames(ANOVA)
    fixed.ColNames = colnames(EF)

    if(length(fixed.names) < 2){
      if (is.na(fixed.names))
          fixed.names = c("\\gamma", "\\tau", "\\rho", "\\phi")
    }
    
    fixed.names = fixed.names[1:(length(fixed.ColNames)/2)]

     #check for interaction effects
    if(any(grepl(":",  fixed.ColNames))){
      fixed.names[grep(":",  fixed.ColNames[1:(length(fixed.ColNames)/2)])]=
      sapply(strsplit(fixed.ColNames[grep(":",  fixed.ColNames[1:(length(fixed.ColNames)/2)])], ":"),
        function(x)  paste(fixed.names[match(x,fixed.ColNames)], collapse = ""))
    }
    
    if(length(grep("MS", random.ColNames)) == 1){ 
      ColNamesLen = length(random.ColNames) - 1
    } else {
      ColNamesLen = length(random.ColNames)
    }
    
    if(length(ColNamesLen)> 2){
    random.names = c("\\sigma^2", paste("\\sigma_{", 
        substr(sapply(strsplit(random.ColNames[3:ColNamesLen], "\\:"), 
                  function(x) x[length(x)]), 1, 1), "}^2", sep = ""))
    } else {
      random.names = c("\\sigma^2")    
    }
    
    ANOVA = ifelse(ANOVA == "0", "", ANOVA)
    finalTable = cbind(ANOVA, EF[match(matchRowNames, rownames(EF)), ])

     #avoid repeat in the fixed componenets
    for(i in 1:nrow(finalTable)){
      finalTable[i, (ncol(ANOVA) +1):ncol(finalTable)] = EF[grep(matchRowNames[i], rownames(EF))[1],]
      rownames(EF)[grep(matchRowNames[i], rownames(EF))[1]] <- ""
    }
    
    finalTable = ifelse(is.na(finalTable), "", finalTable)

    output = "\\begin{table}[ht]\n\\centering\n \\caption{Theoretical ANOVA table}\n"
    
     if(length(grep("MS", random.ColNames)) == 1){ 
      output = c(output, paste("\\begin{tabular}[t]{lrll", paste( rep("l", length(fixed.names)), collapse = ""),"} \n", sep = ""))
    }else{ 
      output = c(output, paste("\\begin{tabular}[t]{lrl", paste( rep("l", length(fixed.names)), collapse = ""),"} \n", sep = ""))
    }
      
    output = c(output, "\\toprule \n")

    firstRow = paste("\\multicolumn{1}{l}{\\textbf{Source of Variation}} & \\multicolumn{1}{l}{\\textbf{DF}} & \\multicolumn{1}{l}{\\textbf{EMS}}&", sep = "")
      
    if(length(grep("MS", random.ColNames)) == 1) 
      firstRow = rbind(firstRow, "\\multicolumn{1}{l}{\\textbf{MS}} &")
    
    firstRow = rbind(firstRow, paste(paste("\\multicolumn{1}{l}{$\\bm{E_{", fixed.names,
        sep = "", collapse = "}}$}&"), "}}$}\\\\ \n", sep = ""))

    output = c(output, firstRow)
    output = c(output, "\\midrule \n")
    for (i in 1:length(matchRowNames)) {
        
        # row names
        SV = matchRowNames[i]
        SV = gsub("   ", "\\\\quad ", SV)
        
        # DF
        DF = paste("$", finalTable[i, 1], "$", sep = "")
        DF = ifelse(DF == "$$", "", DF)

        # Random VC
        coef.VC = finalTable[i, 2:(1 + length(random.names))]
        random.VC = ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "")
                next
            if (coef.VC[j] == 1)
                coef.VC[j] = ""

            random.VC = paste(random.VC, coef.VC[j], random.names[j], sep = "")

            random.VC = paste(random.VC, "+", sep = "")
        }

        # Fixed VC
         if(length(grep("MS", random.ColNames)) == 1){ 
        coef.VC = finalTable[i, (3 + length(random.names)):(2 + length(random.names) + length(fixed.names))]
        } else{
        coef.VC = finalTable[i, (2 + length(random.names)):(1 + length(random.names) + length(fixed.names))]
        }
        
        fixed.VC = ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "")
                next
                
            if (coef.VC[j] == 1)
                coef.VC[j] = ""

            fixed.VC = paste(fixed.VC, coef.VC[j], "\\theta_{", fixed.names[j], "}", sep = "")

            fixed.VC = paste(fixed.VC, "+", sep = "")
        }

        total.VC = ifelse(random.VC == "", "", paste("$", random.VC, "+", fixed.VC, "$", sep = ""))
        total.VC = ifelse(fixed.VC == "", paste("$", random.VC, "$", sep = ""), total.VC)
        total.VC = gsub("\\+\\+", "\\+", total.VC )
        total.VC = gsub("\\+\\$", "\\$", total.VC )
        total.VC = gsub("\\$\\$", "", total.VC )
         
        eff = paste("$", finalTable[i, (2 + length(random.names) + length(fixed.names)):ncol(finalTable)], "$", sep = "", collapse = " & ")

        eff = gsub("\\$\\$", "", eff)


         finalTable[i, 2:(1 + length(random.names))]
        
        if(length(grep("MS", random.ColNames)) == 1){ 
          currentRow = paste(SV, " & ", DF, " & ", total.VC," & ", finalTable[i, 2+length(random.names)] , eff, "\\\\ \n", sep = "")
        } else {        
          currentRow = paste(SV, " & ", DF, " & ", total.VC, " &",eff, "\\\\ \n", sep = "")
        }
        
        output = c(output, currentRow)
    }

    output = c(output, "\\bottomrule \n \\end{tabular} \n \\label{tab:} \n\\end{table} \n")
    output = c(output, "")

    return(cat(output))
}

