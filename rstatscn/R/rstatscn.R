library(jsonlite)
library(httr)

statscnbase<-'http://data.stats.gov.cn/easyquery.htm'
rstatscnEnv<-new.env()

#' private function for sec
#' 
#' @return milsec
milSec<-function()
{
	tt=Sys.time()
	ts=format(tt,"%s")
	ts=paste(ts,"000",sep="")
	ts
}

#' the available dbs
#' 
#' the available dbs in the national db
#' @return a data frame with 2 columns , one is the dbcode, another is the db description 
#' @export 
#' @examples 
#'  statscnDbs()
statscnDbs<-function()
{
	dbs <- c("hgnd","hgjd","hgyd","fsnd","fsjd","fsyd","csnd","csyd","gjnd","gjyd","gjydsdj")
	dbnames <- c("national data, yearly","national data,  quaterly","national data, monthly",
		     "province data, yearly","province data, quaterly","province data, monthly",
                     "city data, yearly","city data, monthly", "international data, yearly", 
                     "international data, monthly","3 main countries data, monthly")
	ret=data.frame(dbcode=dbs,description=dbnames)
	return(ret)
}

#' private function for constructing the query parameter for dfwds
#' 
#' @param wdcode string value , one of c("zb","sj","reg") 
#' @param valuecode string value ,  following is the table for available valuecode
#'    zb:   the valudecode can be gotten by statscnQueryZb() function 
#'    sj:   the valudecode can be "2014" for nd db,  "2014C" for jd db.
#'    reg:  the valudecode is the region code fetched by statscnRegions(dbcode) function
#' @return return the queyr string for the http request
genDfwds<-function(wdcode,valuecode)
{
	if( is.na(valuecode) ){
		return("[]")
	}else{
		paste('[{"wdcode":"',wdcode,'","valuecode":"',valuecode,'"}]',sep="")
	}
}
#' private function for check the http status
#' 
#' @param ret the response obj returned by httr package
#' @return return nothing , but if it finds some error , it stop the script
checkHttpStatus<-function(ret)
{
  if (http_status(ret)$category != "Success") {
    stop(sprintf("Bad response from %s", statscnbase))
  }
}
#' private function to convert the returned jason data to a dataframe
#' 
#' @param rawObj the fromJSON output 
#' @param rowcode rowcode in the data frame
#' @param colcode colcode in the data frame
#' @return the contructed data frame
dataJson2df<-function(rawObj,rowcode,colcode)
{
	ret=rawObj
        if(ret[[1]] != 200) {
		stop("Bad response from the statscn server")
	}
        #dataStructure
	#jj is a list
        #jj[[1]] = 200 #return code
        #jj[[2]] is datanode
        #jj[[2]][[1]] is data
        #jj[[2]][[2]] is description
        #jj[[2]][[2]][,"nodes"][[1]] is row description , it is a dataframe
        #jj[[2]][[2]][,"nodes"][[2]] is col description , it is a dataframe
        desList=ret[[2]][[2]][,'nodes']
	rowWdIdx = which(ret[[2]][[2]]$wdcode == rowcode) 
	colWdIdx = which(ret[[2]][[2]]$wdcode == colcode) 
        rowDes=desList[[rowWdIdx]]
        colDes=desList[[colWdIdx]]

        rowNum=nrow(rowDes) 
        colNum=nrow(colDes) 
        rowNames=rowDes[,1]
        colNames=colDes[,1]

        units=rowDes[,'unit']
        units=ifelse(units=="","",paste("(",units,")",sep=""))
        rowNames=paste(rowNames,units,sep="")

        rowCodes=rowDes[,2]
        colCodes=colDes[,2]
        
	#the rowCode and colCode are in the ret[[2]][[1]][,'wds']
        #it is a list , the list length is the same as the data fetched. list[[1]] is for the first data
	#list[[1]] is a dataframe ,  df[1,'valuecode'] is the rowcode , df[2,'valuecode'] is the colcode
        #now we create a dataframe 
        myret=as.data.frame(matrix(rep(NA,rowNum*colNum),nrow=rowNum))
        rownames(myret)=rowCodes
        colnames(myret)=colCodes
        dfdata=ret[[2]][[1]]
        for (k in seq(1,nrow(dfdata))) {
		wddf=dfdata[k,"wds"][[1]]
		myret[wddf[rowWdIdx,'valuecode'],wddf[colWdIdx,'valuecode']] = dfdata[k,'data'][1,'data']
	}
        rownames(myret)=rowNames
        colnames(myret)=colNames
	return(myret)
}

#' the data categories
#' 
#' the sub data categories for the zbid category, dbcode need to be specified, where the dbcode can be fetched by function
#' statscnDbs(). In the returned data frame, the column 'isParent' shows if each sub category is leap category or not
#' @param zbid the father zb/category id , the root id is 'zb'
#' @param dbcode which db will be queried
#' @return the data frame with the sub zbs/categories , if the given zbid is not a Parent zb/category, null list is returned
#' @export 
#' @examples 
#'  statscnQueryZb()
#'  statscnQueryZb('A01',dbcode="hgnd")
statscnQueryZb<-function(zbid="zb",dbcode="hgnd")
{
	curQuery=list(id=zbid,dbcode=dbcode,wdcode="zb",m="getTree")
	yy<-POST(statscnbase,body=curQuery,encode="form")
	assign('lastQuery',curQuery, envir=rstatscnEnv)
        checkHttpStatus(yy)
        jj=fromJSON(content(yy,"text",encoding="utf-8"))
        return(jj)
}
#' the regions in db
#' 
#' the available regions in the specified db, it is used for query the province, city and country code generally
#' @param dbcode the dbcode should be some province db(fs*) , city db(cs*) or internaltional db(gj*)
#' @return the data frame with all the available region codes and names in the db
#' @export 
#' @examples 
#'  statscnRegions('fsnd')
#'  statscnRegions('csnd')
#'  statscnRegions('gjnd')
statscnRegions<-function(dbcode='fsnd')
{
        curQuery<-list(
		m="getOtherWds",
		dbcode=dbcode,
		rowcode="zb",
		colcode="sj",
		wds="[]",
		#dfwds="[]",
		k1=milSec()
	)
        yy<-GET(statscnbase, query=curQuery)
	assign('lastQuery',curQuery, envir=rstatscnEnv)
        checkHttpStatus(yy)
        ret=fromJSON(content(yy,"text",encoding="utf-8"))
        regIndex <- which(ret[[2]]$wdcode == 'reg')
	df <- ret[[2]][,'nodes'][[regIndex]]
	df$sort=NULL
	colnames(df) <- c("regCode","name")
	return(df)
}
#' query data in the statscn db
#' 
#' the main function for querying the statscn database, it will retrieve the data from specified db and orginize the data in a data frame.
#' @param zb the zb/category code to be queried
#' @param dbcode the db code for querying
#' @param rowcode rowcode in the returned data frame
#' @param colcode colcode in the returned data frame
#' @param moreWd more constraint on the data
#'        where the name should be one of c("reg","sj") , which stand for region and sj/time.
#'        the valuecode for reg should be the region code queried by statscnRegions()
#'        the valuecode for sj should be like '2014' for *nd , '2014C' for *jd , '201405' for *yd.
#'        Be noted that , the moreWd name should be different with either rowcode or colcode
#' @return the data frame you are quering
#' @export 
#' @examples 
#' df=statscnQueryData('A0201',dbcode='hgnd')
#' df=statscnQueryData('A0201',dbcode='fsnd',rowcode='zb',colcode='sj',
#'                     moreWd=list(name='reg',value='110000'))
statscnQueryData<-function(zb="A0201",dbcode="hgnd",rowcode='zb',colcode='sj',moreWd=list(name=NA,value=NA))
{
        curQuery<-list(
		m="QueryData",
		dbcode=dbcode,
		rowcode=rowcode,
		colcode=colcode,
		wds=genDfwds(moreWd$name,moreWd$value),
		dfwds=genDfwds("zb",zb),
		k1=milSec()
	)
        yy<-GET(statscnbase, query=curQuery)
	assign('lastQuery',curQuery, envir=rstatscnEnv)
        checkHttpStatus(yy)
        ret=fromJSON(content(yy,"text",encoding="utf-8"))
        return(dataJson2df(ret,curQuery$rowcode,curQuery$colcode))
}
#' fetch the lastN data
#' 
#' fetch the lastN data for the latest query, only affect the number of rows in the returned data.
#' This function can not be used alone , statscnQueryData() has to be called before this function
#' @param n the number of rows to be fetched
#' @return the last n rows data in the latest query
#' @export 
#' @examples 
#' df=statscnQueryData('A0201',dbcode='hgnd')
#' df2=statscnQueryLastN(20)
statscnQueryLastN<-function(n)
{
	wdcode="sj"	
	valuecode=paste("LAST",n,sep="")
	if( is.null(get('lastQuery', envir=rstatscnEnv)) ){
		stop("please call a statscnQueryData for some data firstly")
	}
	curQuery=get('lastQuery', envir=rstatscnEnv)
	if( curQuery$m=="QueryData" ) {
		curQuery$dfwds=genDfwds(wdcode,valuecode)
	}
        yy<-GET(statscnbase, query=curQuery)
	assign('lastQuery',curQuery, envir=rstatscnEnv)
        checkHttpStatus(yy)
        ret=fromJSON(content(yy,"text",encoding="utf-8"))
        return(dataJson2df(ret,curQuery$rowcode,curQuery$colcode))
}

