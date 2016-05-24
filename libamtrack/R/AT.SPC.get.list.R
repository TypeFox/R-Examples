AT.SPC.get.list <- function(path.to.files, endian = "little"){
  pattern <- tail(strsplit(path.to.files, split = c("/", "\\"), fixed = TRUE)[[1]], 1)
  path    <- gsub(pattern     = pattern, 
                  replacement = "",
                  x           = path.to.files,
                  fixed       = TRUE)
  files <- list.files(path = path, pattern = pattern)
  df    <- expand.grid( file.name           = paste(path, files, sep = ""),
                        n.depth.steps       = 0,
                        projectile          = "",
                        target.material     = "",
                        energy.MeV.u        = 0,
                        peak.position.g.cm2 = 0,
                        endian              = endian,
                        stringsAsFactors    = FALSE)

  for (i in 1:nrow(df)){
      # i <- 2
      tmp    <- AT.SPC.read(file.name   =  df$file.name[i],
                            flavour     = "C",
                            endian      = endian,
                            header.only = TRUE)
      df$n.depth.steps[i]       <- tmp$n.depth.steps
      df$projectile[i]          <- tmp$projectile
      df$target.material[i]     <- tmp$target.material
      df$energy.MeV.u[i]        <- tmp$energy.MeV.u
      df$peak.position.g.cm2[i] <- tmp$peak.position.g.cm2
      cat("Read spc header", i, "of", nrow(df), ".\n")
    }
  return(df)
}