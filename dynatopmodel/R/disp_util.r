# active the dev.num(th) device of the required type. if dev.id is supplied, attempt to
# activate the window with that id
activate.device <- function(dev.num,
                            dev.id=NULL,
                            dev=getOption("device"),
                            get.title.fun=function(id){return("")})
{
  old.dev <- getOption("device")
  options("device"=dev)
  on.exit(options(device=old.dev))

  if(!is.null(dev.id))
  {
    # null device is always #1 so we need to exclude that as a allowable id
    if(dev.id==1){dev.id <- 2}
    idev <- which(dev.list()==dev.id)
    if(length(idev)==0)
    {
      # dev.lis numbers the open devices starting with 2
      max.win <- max(dev.list(), 1)
      if(dev.id > max.win)
      {
        for(id in (max.win+1):(dev.id))
        {
          dev.new(title=get.title.fun(id))
        }
      }
      #open.devices()
      #		dev.new(title=get.title.fun(id))
      #		dev.set(which=dev.id)
      #	open.devices(n=idev, dev=dev)
    } else
    {
    }
    dev.set(dev.id)
  }
  else
  {
    #	dev.nm.2 <-
    nms <- names(dev.list())
    # replace synonyms
    nms <- gsub("windows","X11", nms)
    devs <- dev.list()[grep(dev, nms)]
    if(dev.num > length(devs))
    {
      open.devices(n=dev.num, dev=dev)
    }

    idev <- dev.num
    try(res <- dev.set(devs[idev]),silent=T)
  }
  return(dev.cur())
}

# open enough windows to display storages and spatial output
open.devices <- function(n=1, title="", width=8, height=8,
                         dev=getOption("device"))
{
  # determine how many of devices of the given type are open
  dev.nm.2 <- gsub("X11", "windows" , dev)
  devs <- which(names(dev.list())==dev.nm.2)

  # determine number of new windows required by comparing list of open devices with default device opened by call to dev.new
  n.dev <- length(devs)
  n.new <- max(0, n-n.dev)

  i<-0
  old.dev <- getOption("device")
  options("device"=dev)
  on.exit(options(device=old.dev))
  while(i<n.new)
  {
    # open required number of devices of  type
    try(dev.new(title=title, width=width,height=height), silent=F)
    i <- i+1
  }
  id <- NULL
  # return id of nth window of default type
  try(id<-dev.list()[devs][n])
  return(id)

}
