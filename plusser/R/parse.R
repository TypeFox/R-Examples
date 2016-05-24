##' Parsing a Google+ post
##' 
##' This function turns a Google+ post into a (1 row) data frame extracting or
##' computing a number of fields. See \code{Details}.
##' 
##' This function extracts or computes the following fields:
##' \describe{
##'   \item{\code{ti}}{Date and time the post was published.}
##'   \item{\code{age}}{The age of the post as difference between now and
##'                     \code{ti} in (floating point) days.}
##'   \item{\code{id}}{The post's unique Google+ post ID.}
##'   \item{\code{au}}{The post's author's Google+ user ID.}
##'   \item{\code{ve}}{The action describing the post.}
##'   \item{\code{nC}}{The number of comments the post has attracted so far.}
##'   \item{\code{nP}}{The number of +1s the post has attracted so far.}
##'   \item{\code{nR}}{The number of times the post has been reshared so far.}
##'   \item{\code{at}}{The type of attachment (article, photo, ...) as reported by the API}
##'   \item{\code{msg}}{The post's content.}   
##' }
##'   
##' @param p A raw post as returned from e.g. the \code{\link{harvestPage}}
##'   function.
##' @return A 1 row data frame filled with the information from the post parsed.
##' @export
##' @examples
##' \dontrun{
##' myPosts <- harvestPage("115046504166916768425", ret="list")
##' myPosts.df <- ldply(myPosts, parsePost)
##' }
parsePost <- function(p) {
  ti <- ymd_hms(p$published)
  tod <- hour(ti)
  if (tod >= 6 & tod < 9) {
    todf <- "earlyMorning"
  } else if (tod >= 9 & tod < 12) {
    todf <- "lateMorning"
  } else if (tod >= 12 & tod < 17) {
    todf <- "afternoon"
  } else if (tod >= 17 & tod < 20) {
    todf <- "evening"
  } else if (tod >= 20 | tod < 1) {
    todf <- "night"
  } else {
    todf <- "sleep"
  }
  id <- p$id
  au <- p$actor$id
  ve <- p$verb
  msg <- as.character(p$object$content)
  msg <- gsub("<.*?>", "", msg) ## removing all HTML tags from a G+ msg
  nC <- p$object$replies$totalItems
  nP <- p$object$plusoners$totalItems
  nR <- p$object$resharers$totalItems
  at <- unlist(ifelse(test = length(p$object$attachments)!=0, yes = p$object$attachments[[1]]["objectType"], no = NA))
  df <- data.frame(ti=ti,
                   age=as.numeric((now() - ti), units="days"),
                   id=id,
                   au=au,
                   ve=ve,
                   nC=nC,
                   nP=nP,
                   nR=nR,
                   at=at,
                   msg=msg)
  return(df)
}

##' Parse profile
##' 
##' This function takes a raw list of profiles (as retrieved by
##' \code{\link{harvestProfile}}) and parses the contained information into a
##' one row data frame.
##' 
##' The following fields will be filled with data (if available) or \code{NA}
##' otherwise:
##' \describe{
##'   \item{\code{id}}{The Google+ user ID.}
##'   \item{\code{sex}}{The user's gender: \code{male}, \code{female}, or 
##'                     \code{other}.}
##'   \item{\code{ln}}{The user's last name.}
##'   \item{\code{fn}}{The user's first name.}
##'   \item{\code{verified}}{Logical. \code{TRUE} if it is a verified Google+ 
##'                          profile.}
##'   \item{\code{website}}{A URL listed in the profile.}
##'   \item{\code{ageMin, ageMax}}{Google+ provides only age ranges for some
##'                                profiles. This will contain the lower and
##'                               upper bound of the age range of the user.}
##'   \item{\code{bday}}{The birthday of the user (YYYY-MM-DD).}
##'   \item{\code{nCircled}}{The number of Persons the user circled by.}
##'   \item{\code{currentLoc}}{The user's current location.}
##'   \item{\code{lang}}{The primary language the user reported.}
##'   \item{\code{p1count}}{The number of +1s the user awarded.}
##'   \item{\code{relationship}}{The user's relationship status.}
##'   \item{\code{bio}}{The `About Me' short autobiography.}
##'   \item{\code{tagline}}{The tagline of a profile.}
##'   \item{\code{type}}{The type of a profile: \code{person} or \code{page}.}
##'   \item{\code{brag}}{The `bragging rights' section of the profile.}
##'   \item{\code{occ}}{The person's occupation.}
##'   \item{\code{skills}}{The person's skills.}
##'   }
##' 
##' @param p a raw profile as retrieved e.g. by \code{\link{harvestProfile}}.
##' @return a one row data frame with a number of fields. See Details.
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/people/get}
##' @export
##' @examples
##' \dontrun{
##' gProfile <- parseProfile(harvestProfile("+google", parseFun=NULL))
##' }
parseProfile <- function(p) {
  urls <- p$urls
  if (is.null(urls)) {
    ws <- NA
  } else {
    ut <- sapply(urls, function(x) x["type"]) == "website"
    if (any(ut)) {
      ws <- urls[ut][[1]]["value"]
    } else {
      ws <- NA
    }
  }
  this.ext <- list(id=p$id,
                   sex=p$gender,
                   ln=p$name[1],
                   fn=p$name[2],
                   verified=p$verified,
                   website=ws,
                   ageMin=p$ageRange[2],
                   ageMax=p$ageRange[1],
                   bday=p$birthday,
                   nCircled=p$circledByCount,
                   currentLoc=p$currentLocation,
                   lang=p$language,
                   p1count=p$plusOneCount,
                   relationship=p$relationshipStatus,
                   bio=p$aboutMe,
                   tagline=p$tagline,
                   type=p$objectType,
                   brag=p$braggingRights,
                   occ=p$occupation,
                   skills=p$skills)
  this.ext[sapply(this.ext, is.null)] <- NA
  this.ext <- as.data.frame(this.ext, stringsAsFactors=FALSE)
  rownames(this.ext) <- NULL
  return(this.ext)
}

##' Parsing of post attachments
##' 
##' This function takes a raw list of posts (as retrieved by \code{\link{harvestPage}}) and extracts any attachments it might find.
##' It uses private (i.e. not exported) parsing functions for some known attachment types.
##' At present, these are articles, albums, photos and videos.
##' Other attachment types will just be cast generically to data.frames.
##' The rownames of all these data frames are the \code{id}s of the posts that attachment belongs to.
##' The columns of the returned data frames should be pretty much self-explanatory.
##' If in doubt, check the Google+ API documentation \url{https://developers.google.com/+/api/latest/activities#object.attachments}.
##' @param pl a posting list as retrieved e.g. by \code{\link{harvestPage}}.
##' @return A list containing one data frame per identified attachment type.
##' @export
##' @examples
##' \dontrun{
##' myPosts <- harvestPage("115046504166916768425", ret="list")
##' myPosts.att <- parseAttachments(myPosts)
##' }
parseAttachments <- function(pl) {
  al <- ldply(pl, function(p) unlist(p$object$attachments[[1]]["objectType"]))[,1]
  res <- vector("list", length=length(unique(al)))
  names(res) <- unique(al)
  for (p in pl) {
    if (length(p$object$attachments)>0) { # if no attachment, do nothing
      this.id <- p$id
      this.at <- p$object$attachments[[1]]
      this.ot <- unlist(this.at["objectType"])
      class(this.at) <- c(paste0("a", this.ot), "list")
      res[[this.ot]] <- rbind(res[[this.ot]], parseA(this.at))
      rownames(res[[this.ot]])[length(rownames(res[[this.ot]]))] <- this.id
    }
  }
  return(res)
}

parseA <- function(at) {
  UseMethod(generic = "parseA")
}

parseA.default <- function(at) {
  return(as.data.frame(at))
}

parseA.aarticle <- function(at) {
  name <- unlist(at["displayName"])
  content <- unlist(at["content"])
  url <- unlist(at["url"])
  fiurl <- unlist(at["fullImage"])[1]
  fitype <- unlist(at["fullImage"])[2]
  res <- data.frame(name=ifelse(is.null(name), NA, name), 
                    content=ifelse(is.null(content), NA, content), 
                    url=ifelse(is.null(url), NA, url), 
                    fiurl=ifelse(is.null(fiurl), NA, fiurl),
                    fitype=ifelse(is.null(fitype), NA, fitype))
  rownames(res) <- NULL
  return(res)
}

parseA.aalbum <- function(at) {
  name <- unlist(at["displayName"])
  id <- unlist(at["id"])
  url <- unlist(at["url"])
  nT <- length(at[["thumbnails"]])
  res <- data.frame(name=ifelse(is.null(name), NA, name),
                    id=ifelse(is.null(id), NA, id), 
                    url=ifelse(is.null(url), NA, url), 
                    nT=ifelse(is.null(nT), NA, nT))
  rownames(res) <- NULL
  return(res)
}

parseA.aphoto <- function(at) {
  name <- unlist(at["displayName"])
  id <- unlist(at["id"])
  content <- unlist(at["content"])
  url <- unlist(at["url"])
  imgurl <- unlist(at["image"])[1]
  imgtype <- unlist(at["image"])[2]
  fiurl <- unlist(at["fullImage"])[1]
  fitype <- unlist(at["fullImage"])[2]
  fiheight <- unlist(at["fullImage"])[3]
  fiwidth <- unlist(at["fullImage"])[4]
  res <- data.frame(name=ifelse(is.null(name), NA, name), 
                    id=ifelse(is.null(id), NA, id), 
                    content=ifelse(is.null(content), NA, content), 
                    url=ifelse(is.null(url), NA, url), 
                    imgurl=ifelse(is.null(imgurl), NA, imgurl),
                    imgtype=ifelse(is.null(imgtype), NA, imgtype), 
                    fiurl=ifelse(is.null(fiurl), NA, fiurl), 
                    fitype=ifelse(is.null(fitype), NA, fitype), 
                    fiheight=ifelse(is.null(fiheight), NA, fiheight), 
                    fiwidth=ifelse(is.null(fiwidth), NA, fiwidth))
  rownames(res) <- NULL
  return(res)
}

parseA.avideo <- function(at) {
  name <- unlist(at["displayName"])
  content <- unlist(at["content"])
  url <- unlist(at["url"])
  iurl <- unlist(at["image"])[1]
  itype <- unlist(at["image"])[2]
  iheight <- unlist(at["image"])[3]
  iwidth <- unlist(at["image"])[4]
  eurl <- unlist(at["embed"])[1]
  etype <- unlist(at["embed"])[2]
  res <- data.frame(name=ifelse(is.null(name), NA, name), 
                    content=ifelse(is.null(content), NA, content),
                    url=ifelse(is.null(url), NA, url), 
                    iurl=ifelse(is.null(iurl), NA, iurl), 
                    itype=ifelse(is.null(itype), NA, itype),
                    iheight=ifelse(is.null(iheight), NA, iheight), 
                    iwidth=ifelse(is.null(iwidth), NA, iwidth), 
                    eurl=ifelse(is.null(eurl), NA, eurl), 
                    etype=ifelse(is.null(etype), NA, etype))
  rownames(res) <- NULL
  return(res)
}
