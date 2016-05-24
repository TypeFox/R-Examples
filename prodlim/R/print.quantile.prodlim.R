##' @export 
print.quantile.prodlim <- function(x,digits=2,na.val="--",...){
    printx <- function(u){
        ifelse(is.na(u),na.val,round(u,digits))
    }
    lapply(1:length(x),function(i){
        tab <- x[[i]]
        cat("\n")
        if(length(names(x)[i])){cat(names(x)[i],"\n\n")}
        if(0.5 %in% tab$q ||all(c(0.25,0.75) %in% tab$q)){
            if(0.5 %in% tab$q){
                cat("Median time (CI-95%): ",
                    printx(tab[tab$q==0.5,"quantile"]),
                    " (",
                    printx(tab[tab$q==0.5,"lower"]),
                    ";",
                    printx(tab[tab$q==0.5,"upper"]),
                    ")",
                    "\n",sep="")
            }
            if(all(c(0.25,0.75) %in% tab$q)){
                if (attr(x,"model")=="survival")
                    cat("IQR (time):",
                        " (",
                        printx(tab[tab$q==0.75,"quantile"]),
                        ";",
                        printx(tab[tab$q==0.25,"quantile"]),
                        ")",
                        "\n",sep="")
                else
                    cat("IQR (time):",
                        " (",
                        printx(tab[tab$q==0.25,"quantile"]),
                        ";",
                        printx(tab[tab$q==0.75,"quantile"]),
                        ")",
                        "\n",sep="")
                    
            }
        }
        else{
            print(tab,...)
        }
        ## cat("Range (time):",
        ## " (",
        ## printx(tab[tab$q==1,"quantile"]),
        ## ";",
        ## printx(tab[tab$q==0,"quantile"]),
        ## ")",
        ## "\n")
    })
    invisible(x)
}
