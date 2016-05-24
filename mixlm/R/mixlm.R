# Startup
# Set contrasts and report
.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  options(contrasts=c('contr.sum','contr.poly'))
  if (interactive()) {
    packageStartupMessage('mixlm ',as.character(utils::packageVersion("mixlm")))
    packageStartupMessage('Setting default contrast to sum-to-zero.')
  }
}

# Effect labels
effect.labels <- function(t,data){
	csum <- ifelse(options("contrasts")[[1]][1] == "contr.sum",	TRUE, FALSE)
	effects   <- attr(t, "term.labels")
	factors   <- attr(t, "factors")
	intercept <- attr(t, "intercept")
	n.eff     <- length(effects)
	if(n.eff==0){
		return(NULL)
	}
	split.effects <- strsplit(effects,":")
	names     <- "(Intercept)"

	for(i in 1:n.eff){ 
		cur <- split.effects[[i]]
		if(i == 1 && intercept == 0){
			levs <- levels(data[[cur]])
			names <- paste(cur,"(",levs,")",sep="")
		} else {
			inter <- list()
			for(j in 1:length(cur)){
				if("factor"%in%class(data[[cur[j]]])){ # Handle factor main effect
					levs <- levels(data[[cur[j]]])
					if(csum){
						n.lev <- length(levs)
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							# inter[[j]] <- paste(cur[j],"(",levs[-n.lev],"-",levs[n.lev],")",sep="")
							inter[[j]] <- paste(cur[j],"(",levs[-n.lev],")",sep="")
						}
					} else {
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							inter[[j]] <- paste(cur[j],"(",levs[-1],")",sep="")
						}
					}			
				} else {
					inter[[j]] <- cur[j]
				}
			}
			names <- c(names, apply(expand.grid(inter),1,paste,sep="",collapse=":"))
		}
	}
	names
}


effect.source <- function(t,data){
  csum <- ifelse(options("contrasts")[[1]][1] == "contr.sum",	TRUE, FALSE)
  effects   <- attr(t, "term.labels")
  factors   <- attr(t, "factors")
  intercept <- attr(t, "intercept")
  n.eff     <- length(effects)
  if(n.eff==0){
    return(NULL)
  }
  split.effects <- strsplit(effects,":")
  names     <- "(Intercept)"
  
  for(i in 1:n.eff){ 
    cur <- split.effects[[i]]
    if(i == 1 && intercept == 0){
      levs <- levels(data[[cur]])
      names <- paste(cur,"(",levs,")",sep="")
    } else {
      inter <- list()
      for(j in 1:length(cur)){
        if("factor"%in%class(data[[cur[j]]])){ # Handle factor main effect
          levs <- levels(data[[cur[j]]])
          if(csum){
            n.lev <- length(levs)
            if(factors[cur[j],i]==2){
              inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
            } else {
              # inter[[j]] <- paste(cur[j],"(",levs[-n.lev],"-",levs[n.lev],")",sep="")
              inter[[j]] <- paste(cur[j],"(",levs[-n.lev],")",sep="")
            }
          } else {
            if(factors[cur[j],i]==2){
              inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
            } else {
              inter[[j]] <- paste(cur[j],"(",levs[-1],")",sep="")
            }
          }			
        } else {
          inter[[j]] <- cur[j]
        }
      }
      tmp <- apply(expand.grid(inter),1,paste,sep="",collapse=":")
      names <- c(names, rep(paste(cur,sep="",collapse=":"),length(tmp)))
    }
  }
  names
}

