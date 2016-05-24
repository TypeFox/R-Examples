#
#  Copyright 2014 jSonar Inc
#  All Rights Reserved.
#
#  Licensed under the GNU AFFERO GENERAL PUBLIC LICENSE version 3
#  See http://www.r-project.org/Licenses/AGPL-3
#

SonarConnection <- setClass('SonarConnection', slots=c('url'='character', 'params'='character'));

#' JSON Studio connection
#'
#' Create a connection to a Mongo database through JSON Studio
#' 
#' This function returns a SonarConnection object which can be used with
#' \code{\link{sonarFind}} and \code{\link{sonarAgg}} to query a Mongo
#' database.
#'
#' The parameters for this function are explained in greater detail in the
#' JSON Studio help page \emph{Using the Gateway}.
#'
#' @param url the url where JSON Studio can be accessed
#' @param host the hostname of the Mongo server, as it would be entered from
#'   the JSON Studio login screen
#' @param db the name of the database you intend to access
#' @param port the port number where Mongo is running
#' @param username a username to log in to the database, if necessary
#' @param pwd a password to log in to the database, if necessary
#' @param sdb the name of a database to store JSON Studio-related collections
#' @param ssl TRUE to connect using SSL
#' @param anyCert TRUE to accept any SSL certificate
#' @param krb TRUE to authenticate using Kerberos
#' @param mapCredentials TRUE to map credentials to a functional user account
#'   with which to access data
#' @param secondaryPref TRUE to allow connecting to a secondary of a replica
#'  set and prefer a secondary if the host value passed in is a replica set
#' @return A SonarConnection object to connect to the given Mongo database
#'   through JSON Studio, which can be used with \code{\link{sonarFind}} or
#'   \code{\link{sonarAgg}}.
#'
#' @examples
#' con <- new.SonarConnection('https://localhost:8443', 'localhost', 'test')
#'
#' @export
#' @keywords connection database
#'
#' @family connection
#' @seealso \url{http://jsonstudio.com/wp-content/uploads/2014/04/manual141/_build/html/index.html}
new.SonarConnection <- function(url, host, db, port=27017, username=NULL, pwd=NULL, sdb=NULL, ssl=FALSE, anyCert=FALSE, krb=FALSE, mapCredentials=FALSE, secondaryPref=FALSE)
{
    connectionParams <- list(
        'host' = host,
        'port' = port,
        'db' = db
    );

    # we construct the list this way so null-valued items are implicitly ignored
    connectionParams[['username']] <- username;
    connectionParams[['pwd']] <- pwd;
    connectionParams[['sdb']] <- sdb;
    if (ssl) connectionParams[['ssl']] <- ssl;
    if (anyCert) connectionParams[['anyCert']] <- anyCert;
    if (krb) connectionParams[['krb']] <- krb;
    if (mapCredentials) connectionParams[['mapCredentials']] <- mapCredentials;
    if (secondaryPref) connectionParams[['secondaryPref']] <- secondaryPref;

    connection <- SonarConnection(
        'url' = url,
        'params' = paste(names(connectionParams), connectionParams, sep='=', collapse='&')
    );

    return(connection);
}

sonarGatewayURL <- function(connection, queryName, queryCol, output='csv', type='find', bind=list(), limit=NULL, publishedBy=NULL)
{
    if(! inherits(connection, 'SonarConnection'))
        stop("connection must be a SonarConnection. Create one with new.SonarConnection");

    gatewayParams <- list(
        'output' = output,
        'type' = type,
        'name' = queryName,
        'col' = queryCol
    );
    
    gatewayParams[['limit']] <- limit;
    gatewayParams[['published_by']] <- publishedBy;

    # encode bind parameter types
    for(key in names(bind)) {
        bindParam <- paste('bind', key, sep='.');
        bindValue <- encodeBindValue(bind[[key]]);
        gatewayParams[[bindParam]] <- bindValue;
    }

    queryArgs <- paste(names(gatewayParams), gatewayParams, sep='=', collapse='&');
    allArgs <- paste(connection@params, queryArgs, sep='&');
    gatewayURL <- sprintf("%s/Gateway?%s", connection@url, allArgs);

    return(gatewayURL);
}

encodeBindValue <- function(value)
{
    if(inherits(value, 'character'))
    {
        return(paste('"', value, '"', sep=''));
    }
    else
    {
        return(value);
    }
}
