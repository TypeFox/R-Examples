##' A function to find the domain, element and item code for a
##' specific FAOSTAT query.
##'
##' @export

FAOsearch = function(){
    with(FAOmetaTable, {
        ## Find the Group code
        gc = NA
        ## while loop iterates until a valid value is supplied
        while(length(gc)==0 || is.na(gc)){
            cat(paste(paste("(", 1:length(groupTable$groupName), ") ",
                            groupTable$groupName, sep = ""), collapse = "\n"))
            gcn = readline("\nWhich Group are you looking for: ")
            gc = groupTable[as.numeric(gcn), "groupCode"]
        }

        ## Find the Domain code
        subdomainTable = subset(domainTable, groupCode == gc)
        dc = NA
        ## while loop iterates until a valid value is supplied
        while(length(dc)==0 || is.na(dc)){
            cat(paste(paste("(", 1:length(subdomainTable$domainName), ") ",
                            subdomainTable$domainName, sep = ""),
                      collapse = "\n"))
            dcn = readline("\nWhich Domain are you looking for: ")
            dc = subdomainTable[as.numeric(dcn), "domainCode"]
        }

        ## Individual or aggregated item
        useAgg = NA
        while(is.na(useAgg) || !useAgg %in% c("0", "1")){
            cat("(0) Individual item (e.g. Apples, Wheat)\n")
            cat("(1) Aggregated item (e.g. Total cereals, Total meat\n")
            useAgg = readline(paste("Are you looking for individual item or",
                              "aggregated item:"))
        }

        if(as.numeric(useAgg)){
            ## Find the Item Aggregated code
            subitemTable = subset(itemAggTable, domainCode == dc)
            ic = NA
            ## while loop iterates until a valid value is supplied
            while(length(ic)==0 || is.na(ic)){
                cat(paste(paste("(", 1:length(subitemTable$itemName), ") ",
                                subitemTable$itemName, sep = ""),
                          collapse = "\n"))
                icn = readline(paste("\nWhich Item are you looking for?",
                                     "('All' for everything):"))
                if(icn == "All")
                    icn = 1:length(subitemTable$itemName)
                ic = subitemTable[as.numeric(icn), "itemCode"]
            }
        } else {

            ## Find the Item code
            subitemTable = subset(itemTable, domainCode == dc)
            ic = NA
            ## while loop iterates until a valid value is supplied
            while(length(ic)==0 || is.na(ic)){
                cat(paste(paste("(", 1:length(subitemTable$itemName), ") ",
                                subitemTable$itemName, sep = ""),
                          collapse = "\n"))
                icn = readline(paste("\nWhich Item are you looking for?",
                                     "('All' for everything): "))
                if(icn == "All")
                    icn = 1:length(subitemTable$itemName)
                ic = subitemTable[as.numeric(icn), "itemCode"]
            }
        }

        ## Find the Element code
        subelementTable = subset(elementTable, domainCode == dc)
        ec = NA
        ## while loop iterates until a valid value is supplied
        while(length(ec)==0 || is.na(ec)){
            cat(paste(paste("(", 1:length(subelementTable$elementName), ") ",
                            subelementTable$elementName, sep = ""),
                      collapse = "\n"))
            ecn = readline(paste("\nWhich Element are you looking for?",
                                 "('All' for everything):"))
            if(ecn == "All")
                ecn = 1:length(subelementTable$elementName)
            ec = subelementTable[as.numeric(ecn), "elementCode"]
        }

        tmp = expand.grid(dc, ic, ec, stringsAsFactors = FALSE)
        colnames(tmp) = c("domainCode", "itemCode", "elementCode")
        tmp = merge(tmp, domainTable[, c("domainCode", "domainName")],
            all.x = TRUE)
        tmp = merge(tmp, subitemTable[, c("itemCode", "itemName")],
            all.x = TRUE)
        final.df = merge(tmp, subelementTable[, c("elementCode", "elementName")],
            all.x = TRUE)
        final.df$name =
            with(final.df, paste(domainName, itemName, elementName, sep = "_"))
        final.df$domainName = NULL
        final.df$itemName = NULL
        final.df$elementName = NULL
        .LastSearch <<- final.df
        cat("\n** Search result saved as .LastSearch**\n")

    }
    )
}

utils::globalVariables(names = c("FAOmetaTable"))