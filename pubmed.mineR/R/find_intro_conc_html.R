find_intro_conc_html = function (y, themes, all) 
{
    s = readline(prompt = "Enter the path: ")
    s = as.character(s)
    target = HTMLInitFile(s, filename = "sample", BackGroundColor = "#FFFFFF")
    for (i in 1:length(themes)) {
        checks1 = searchabsL(y, include = themes[i])
        checks2 = input_for_find_intro_conc_html(checks1, all)

        if (all == FALSE){check1 = lapply(checks2, function(x) {SentenceToken(x[1])})
        check1a = lapply(checks2, function(x) {return(x[2])})} else
         if (all == TRUE ) {check1 = lapply(checks2, function(x) {SentenceToken(x[3])})
          check1a = lapply(checks2, function(x) {return(x[4])}); check1b = lapply(checks2, function(x) {return(x[1])})  ;check1c = lapply(checks2, function(x) {return(x[2])})} 

        check2 = lapply(check1, function(x) {check3 = regexpr("conclusion", x, ignore.case = T)
            check4 = which(check3 != -1)
            if (length(check4) != 0) {check5 = x[check4:length(x)]
                return(c(x[1:5], "=========\n", "=========\n", 
                  check5))
            }
            else {
                if (length(x) > 10) {
                  check5 = x[(length(x) - 3):length(x)]
                  return(c(x[1:5], "=========\n", "=========\n", 
                    check5))
                }
                else return(x)
            }
        })
       
            if (all == FALSE){check6 = lapply(1:length(check1a), function(a) {
            return(c(check2[[a]], paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=", 
                check1a[[a]], sep = "")))})}

              else if (all == TRUE){check6 = lapply(1:length(check1a), function(a) {
            return(c(c(check1b[[a]],check1c[[a]],check2[[a]]), paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=", 
                check1a[[a]], sep = "")))})}

        HTML.title(themes[i], file = target)
        HTML(check6, file = target)
    }
}
