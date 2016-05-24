#' Location Class for Google's Blogger service
#' 
#' This class implements a small subset of the blogger API v3. 
#' The \code{meta} method provides information on the submitted blogposts. The \code{put} method accepts a
#' \code{BlogPostTarget} that can be transfered to Blogger.
#'
#' @examples
#' getSlots("Blogger")
#'
#' @seealso \code{\link{blogger}}, \code{\link{mdreport}}
#'
#' @references \href{https://developers.google.com/blogger/docs/3.0/}{Blogger}
#' @name Blogger-class
#' @rdname Blogger-class
#' @exportClass Blogger
setClass(
    Class="Blogger", 
    representation=representation(
        google.auth="GoogleOAuth2", 
        curl.handle="ANY", # RCurl:::CURLHandle is not exported
        blogid="character", 
        blogtitle="character"
    ),
    contains=c("Location", "Xdata")
)

#' Constructor for Blogger class
#'
#' Instantiates an object and authenticates with google.
#' If the provided \code{oauthfile} parameter points to an existing file,
#' the authentication information is loaded by \code{read.google.oauth2},
#' and the client_id and client_secret information are ignored.
#' If \code{oauthfile} is missing, an initial authentication (directing the
#' user to an website) is performed by \code{google.oauth2}.
#'
#' @param blogurl            URL of the (existing) blog. Defaults to getOption("blogger.blog").
#' @param oauthfile          filename of previously saved authentication information. 
#' @param client_id          client_id. See \code{google.oauth2}
#' @param client_secret      client_secret. See \code{google.oauth2}
#' @param clss               name of the class for convenient inheritance. Defaults to "Blogger".
#' 
#'
#' @return Blogger 
#' @rdname Blogger-class
#' @export
blogger <- function(
    oauthfile=getOption("blogger.oauthfile"),
    client_id=getOption("datamart.client_id"), 
    client_secret=getOption("datamart.client_secret"), 
    blogurl=getOption("blogger.blog"),
    clss="Blogger"
) {
    oa <- if(file.exists(oauthfile)) 
        read.google.oauth2(oauthfile)
    else
        google.oauth2("blogger", client_id=client_id, client_secret=client_secret)
    ch <- RCurl::getCurlHandle()
    
    ret <- RCurl::getURL(
        sprintf("https://www.googleapis.com/blogger/v3/blogs/byurl?url=%s", blogurl),
        .opts = list(
            httpheader = c("Authorization" = sprintf("Bearer %s", query(oa, "access_token"))),
            ssl.verifypeer = TRUE, 
            ssl.verifyhost = TRUE, 
            cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
        ),
        curl=ch
    )
    ret <- RJSONIO::fromJSON(ret)
    if(!("id" %in% names(ret))) stop("Blog '", blogurl, "' not found.")
    
    new(
        clss,
        google.auth=oa, 
        curl.handle=ch, 
        blogid=ret[["id"]], 
        blogtitle=ret[["name"]]
    )
}

#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,Blogger-method
setMethod(
    f="show",
    signature="Blogger",
    definition=function(object) cat(sprintf("<blog '%s' @ Blogger.com>\n", object@blogtitle))
)

#' @docType methods
#' @rdname meta-methods
#' @name meta
#' @aliases meta meta,Blogger-method
#' @export
setMethod(
  f="meta",
  signature=c(self="Blogger"),
  definition=function(self) {
    uri <- paste(
        "https://www.googleapis.com/blogger/v3/blogs/", self@blogid,
        "/posts?",
        "fields=nextPageToken,items(id,title,published,updated)", 
        sep=""
    )
     
    fetchpage <- function(pageToken=NULL) {
        page_uri <- if(!is.null(pageToken)) sprintf("%s&pageToken=%s", uri, pageToken) else uri
        ret <- RCurl::getURL(
            page_uri,
            .encoding = 'UTF-8', 
            .opts = list(
                httpheader = c(Authorization = sprintf("Bearer %s", query(self@google.auth, "access_token"))),
                followlocation = TRUE, 
                ssl.verifypeer = TRUE, 
                ssl.verifyhost = TRUE, 
                cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
            ),
            curl=self@curl.handle
        )
        ret <- RJSONIO::fromJSON(ret)
        entries <- lapply(ret[["items"]], function(node) unlist(node[c("id", "title", "published", "updated")]))
        entries <- as.data.frame(do.call(rbind, entries), stringsAsFactors=FALSE)
        dat <<- rbind(dat, entries) # side effect
        return(ret[["nextPageToken"]])
    }

    dat <- data.frame()
    np <- fetchpage()
    while(!is.null(np)) np <- fetchpage(np)

    #TODO: unescape &#39; etc. in title column
    # test <- "My this &amp; last year&#39;s resolutions"
    # xp <- gregexpr("&#([0-9]+);", test)
    # gsubfn("&#([0-9]+);", function(x) rawToChar(as.raw(as.numeric(x))), test)
    # gsubfn("&([^;]+);", function(x) rawToChar(as.raw(as.numeric(x))), test)
    #> as.integer(charToRaw("\u27"))
    # [1] 39
    # > rawToChar(as.raw(39))
    # [1] "'"
    rownames(dat) <- dat[,"id"]
    dat$published <- strptime(strhead(dat$published, 19), "%Y-%m-%dT%H:%M:%s")
    dat$updated <- strptime(strhead(dat$updated, 19), "%Y-%m-%dT%H:%M:%s")
    idx <- which(colnames(dat)=="id")
    dat <- dat[, -idx]
    return(dat)    
    
  }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,BlogPostTarget,Blogger-method
setMethod(
    f="put",
    signature=c(target="BlogPostTarget", where="Blogger"),
    definition=function(target, where, verbose=TRUE, ...) {
        entry <- list(
            title=target@subject,
            content=target@content,
            labels=Filter(function(x) x!="",target@label) # TODO: this does not work
        )
        
        # overwrite?
        postid <- ""
        if(target@overwrite) {
            m <- meta(where) # this may take some time, it gets all blog entries.
            mvec <- rownames(m)
            names(mvec) <- m[,"title"]
            postid <- mvec[target@subject] 
            if(is.na(postid)) postid <- "" 
        }
        if(verbose) cat("postid='", postid, "'\n")
        
        if(postid=="") {
            if(verbose) cat("New post..\n")
            uri <- sprintf('https://www.googleapis.com/blogger/v3/blogs/%s/posts', where@blogid)
            if(target@draft) uri <- paste(uri, "?isDraft=true", sep="")
            ret <- RCurl::postForm(
                uri,
                # postfields=toJSON(entry),
                curl=where@curl.handle,
                .opts=list(
                    httpheader = c(
                        Authorization = sprintf("Bearer %s", query(where@google.auth, "access_token")),
                        "Content-Type"="application/json"
                    ),
                    postfields=toJSON(entry),
                    customrequest = "POST",
                    ssl.verifypeer = TRUE, 
                    ssl.verifyhost = TRUE, 
                    cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
                ),
                .contentEncodeFun=RCurl::curlPercentEncode
            )        

        } else {
            if(verbose) cat("post update..\n")
            uri <- sprintf('https://www.googleapis.com/blogger/v3/blogs/%s/posts/%s?fetchBody=false', where@blogid, postid)
            if(target@draft) {
                uri <- paste(uri, "&publish=false", sep="")
            }

            ret <- RCurl::postForm(
                uri,
                curl=where@curl.handle,
                .opts=list(
                    httpheader = c(
                        Authorization = sprintf("Bearer %s", query(where@google.auth, "access_token")),
                        "Content-Type"="application/json"
                    ),
                    customrequest = "PUT", # POST does not work on update
                    postfields=toJSON(entry),
                    ssl.verifypeer = TRUE, 
                    ssl.verifyhost = TRUE, 
                    cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")
                ),
                .contentEncodeFun=RCurl::curlPercentEncode
            )        
        }
        
        #extract id
        ret <- RJSONIO::fromJSON(ret)
        if(!"id" %in% names(ret)) stop("posting article failed.")
        return(ret[["id"]])
    }
)


# Project ID:  rstats-datamart-050  
# Project Number:  679853272845

# auth_uri="https://accounts.google.com/o/oauth2/auth",
    # client_secret="w186HM7bQfk98m-wXxkM_zPe",
    # token_uri="https://accounts.google.com/o/oauth2/token",
    # client_email="",
    # redirect_uris":["urn:ietf:wg:oauth:2.0:oob","oob"],
    # client_x509_cert_url="",
    # client_id="679853272845-bud4ot1va793h9hbnib5clpi0luvvr57.apps.googleusercontent.com",
    # auth_provider_x509_cert_url="https://www.googleapis.com/oauth2/v1/certs"
    #
    # code für karsten W.: 4/khCp5EbGzosZW1SxN8Gj7MbGt3RP.Qoo8Q73lhp4dOl05ti8ZT3bJam1QiAI
    
# Here's the OAuth 2.0 scope information for the Blogger API:
# https://www.googleapis.com/auth/blogger
