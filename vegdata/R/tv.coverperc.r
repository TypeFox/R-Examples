tv.coverperc <- function (db, obs, RelScale, tv_home, tvscale, quiet=FALSE, ...) 
{
  if(missing(tv_home)) {
    tv_home <- tv.home()
  }
  
  if(missing(tvscale)) {
#    print(list.files(file.path(tv_home, "Popup", tv.dict(db))))
    tvscale <- read.dbf(file.path(tv_home, "Popup", tv.dict(db), "tvscale.dbf") )
  }
  tvscale <- tvscale[!is.na(tvscale$SCALE_NR),]
  rownames(tvscale) <- tvscale[, 1]
  if (missing(RelScale)) {
      ow <- options('warn')
      options(warn = -1)
      RelScale <- tv.site(db=db, tv_home=tv_home, verbose = quiet)[, c("RELEVE_NR", "COVERSCALE")]
      options(ow)
      }
  if (missing(obs))
      obs <- tv.obs(db, tv_home, as.is=TRUE)
  obs$COVERSCALE <- RelScale$COVERSCALE[match(obs$RELEVE_NR, RelScale$RELEVE_NR)]
#  obs$COVER_CODE[is.na(obs$COVERSCALE) | obs$COVERSCALE == '9x']
  g <- obs$COVERSCALE
  if(any(is.na(g)))  {
    print(unique(obs[is.na(g),'COVER_CODE']))
    stop('These releves miss a cover scale value in the header data.')
    }
  #### Split ###
  obs <- split(obs, g, drop = FALSE)
  for (i in names(obs)) {
    if (i == "00") {
    	obs[[i]]$COVER_CODE <- replace(as.character(obs[[i]]$COVER_CODE), obs[[i]]$COVER_CODE == '9X', '100')
    	if(any(is.na(as.numeric(obs[[i]]$COVER_CODE)))) stop('Not all percentage cover values in your databse are numeric, please check in Turboveg.')
      obs[[i]] <- data.frame(obs[[i]], COVER_PERC = as.numeric(as.character(obs[[i]][, "COVER_CODE"])))
    }
    else {
      p <- which(is.na(tvscale[i,]))[1]
      if(is.na(p)) p <- ncol(tvscale)
      scala <- tvscale[i,]
      if(is.na(scala[1])) stop('Can not find cover scale "', i, '" in ', file.path('Turbowin','Popup', tv.dict(db),'tvscale.dbf'))
      code <- iconv(t(scala[seq(4,(p-1),2)]), from="CP437", to='UTF-8')
      perc <- scala[seq(5,p,2)][1,]
      d.f <- data.frame(code=code[,1], perc = as.numeric(perc))
      if(!quiet) {
        cat('Cover code used: ',i , as.character(tvscale[i, 2]))
      #  write.table(t(d.f), col.names = FALSE, sep = "\t", quote = FALSE)
      #  print(table(t(d.f), col.names = FALSE, sep = "\t", quote = FALSE))
        print(as.table(t(d.f)), col.names = FALSE, sep = "\t", quote = FALSE)
      }
      obs[[i]]["COVER_PERC"] <- d.f$perc[match(obs[[i]][,"COVER_CODE"], d.f$code)]
  }
  }
  obs <- unsplit(obs, g)
  if(any(is.na(obs$COVER_PERC))) {
      print(obs[is.na(obs$COVER_PERC),'COVER_CODE'])
      stop("Invalid cover codes, please check tvabund.dbf and tvscale.dbf!")
  }
  obs
}

