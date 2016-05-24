write.pdb <- function(x, file="Rpdb.pdb")
{
  if(is.pdb(x)) x <- list(x)
  if(!all(unlist(lapply(x, is.pdb))))
      stop("'x' must be an object of class 'pdb'")

  lines <- NULL

  title <- unique(unlist(lapply(x, function(y) return(y$title))))
  if(!is.null(title))
  {
    title[ substr(title, 1, 6) != "TITLE " ] <- paste0("TITLE ",title[ substr(title, 1, 6) != "TITLE " ])
    lines <- c(lines,title)
  }

  remark <- unique(unlist(lapply(x, function(y) return(y$remark))))
  if(!is.null(remark))
  {
    remark[ substr(remark, 1, 6) != "REMARK" ] <- paste0("REMARK",remark[ substr(remark, 1, 6) != "REMARK" ])
    lines <- c(lines,remark)
  }

  if(!is.null(x[[1]]$cryst1))
  {
    abc <- paste(format(x[[1]]$cryst1$abc,width=9,nsmall=3,justify="right"),collapse="")
    abg <- paste(format(x[[1]]$cryst1$abg,width=7,nsmall=2,justify="right"),collapse="")
    lines <- c(lines,paste("CRYST1",abc,abg,sep=""))
  }

  for(model in seq_along(x)) {
    if(length(x) > 1) lines <- c(lines, paste0("MODEL ",format(model, width = 4)))

    if(basis(x[[model]]) != "xyz") x[[model]] <- abc2xyz.pdb(x[[model]])
  
    X <- round(x[[model]]$atoms$x1, digits = 3)
    Y <- round(x[[model]]$atoms$x2, digits = 3)
    Z <- round(x[[model]]$atoms$x3, digits = 3)
    occ  <- round(x[[model]]$atoms$occ , digits = 2)
    temp <- round(x[[model]]$atoms$temp, digits = 2)
  
    recname <- format(x[[model]]$atoms$recname, justify = "left" , width = 6)
    eleid   <- format(x[[model]]$atoms$eleid  , justify = "right", width = 5)
    elename <- format(x[[model]]$atoms$elename, justify = "left" , width = 4)
    alt     <- format(x[[model]]$atoms$alt    , justify = "left" , width = 1)
    resname <- format(x[[model]]$atoms$resname, justify = "left" , width = 3)
    chainid <- format(x[[model]]$atoms$chainid, justify = "left" , width = 1)
    resid   <- format(x[[model]]$atoms$resid  , justify = "right", width = 4)
    insert  <- format(x[[model]]$atoms$insert , justify = "left" , width = 1)
    X       <- format(X                       , justify = "right", width = 8)
    Y       <- format(Y                       , justify = "right", width = 8)
    Z       <- format(Z                       , justify = "right", width = 8)
    occ     <- format(occ                     , justify = "right", width = 6)
    temp    <- format(temp                    , justify = "right", width = 6)
    segid   <- format(x[[model]]$atoms$segid  , justify = "left" , width = 3)
  
    lines <- c(lines,paste(recname,eleid," ",elename,alt,resname," ",chainid,resid,insert,"   ",X,Y,Z,occ,temp,"      ",segid,sep=""))
    if(length(x) > 1) lines <- c(lines,"ENDMDL")
  }
  
  if(!is.null(x[[1]]$conect))
  {
    eleid.1 <- format(x[[1]]$conect$eleid.1,width=5,justify="right")
    eleid.2 <- format(x[[1]]$conect$eleid.2,width=5,justify="right")
    conect <- paste("CONECT",eleid.1,eleid.2,sep="")
    lines <- c(lines,conect)    
  }
  writeLines(lines,file)
}