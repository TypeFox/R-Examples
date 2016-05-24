ordtaxa <- function (taxa,site)
{
    print(taxa)
    repeat {
        plots <- readline(' enter the plots    : ')
        if (plots == "") {
            break
        } else {
            pnt <- as.numeric(readline(' in front of        : '))
        }

        for (i in strsplit(plots,",")[[1]]){
            ord <- 1:nrow(taxa)
            x <- match(i,row.names(taxa))
            if (!is.na(x)) {
                ord <- ord[-x]
                y <- match(pnt,row.names(taxa[ord,]))
                if (!is.na(y)) {
                        if (y==1) {
                            ord <- c(x,ord)
                        } else {
                            first <- ord[1:(y-1)]
                            last <- ord[y:length(ord)]
                            ord <- c(first,x,last)
                        }
                        taxa <- taxa[ord,]
                        site <- site[ord,]
                        print(taxa)
                    } else {
                        print(paste('plot',pnt,'does not exist'))
                    }
                } else {
                    print(paste('plot',i,'does not exist'))
                }
            }
            repeat {
                species <- readline(' enter the species  : ')
                if (species == "") {
                    break
                } else {
                    pnt <- readline(' in front of        : ')
                }
                for (i in strsplit(species,",")[[1]]){
                    ord <- 1:ncol(taxa)
                    x <- match(i,names(taxa))
                    if (!is.na(x)) {
                        ord <- ord[-x]
                        y <- match(pnt,names(taxa[,ord]))
                        if (!is.na(y)) {
                            if (y==1) {
                                ord <- c(x,ord)
                            } else {
                                first <- ord[1:(y-1)]
                                last <- ord[y:length(ord)]
                                print(first)
                                print(last)
                                ord <- c(first,x,last)
                            }
                            taxa <- taxa[,ord]
                            print(taxa)
                        } else {
                            print(paste('species',pnt,'does not exist'))
                        }
                    } else {
                        print(paste('species',i,'does not exist'))
                    }
                }
            }
        }
    out <- list(taxa=taxa,site=site)
    invisible(out)
}
