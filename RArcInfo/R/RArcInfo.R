get.namesofcoverages <- function(directory) 
{
	directory<-as.character(directory)
#	.Call("get_names_of_coverages", as.character(dir))

	if(length(dir(path=directory, pattern="info"))==1)
	{
		fil<-file.info(dir(path=directory, full.names=TRUE))

		covnames<-(dir(path=directory))[fil$isdir]
		covnames<-covnames[covnames!="info"]
	}
	covnames
}
get.tablenames <-function(infodir) 
{
	data<-.Call("get_table_names", as.character(infodir), PACKAGE="RArcInfo")

	#A data frame with all the data
	data.frame(TableName=I(data[[1]]), InfoFile=I(data[[2]]), NFields=data[[3]], RecSize=data[[4]], NRecords=data[[5]], External=I(data[[6]]))
}
get.tablefields <- function(infodir, tablename) 
{
	data<-.Call("get_table_fields", as.character(infodir), as.character(tablename), PACKAGE="RArcInfo")

	#A data frame with all the data
	data.frame(FieldName=I(data[[1]]), FieldType=data[[2]])
}

get.arcdata <- function(datadir, coverage, filename="arc.adf") 
{
	data<-.Call("get_arc_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")

	#a table (dataframe) with the first seven fields is built
	df<-data.frame(ArcId=data[[1]], ArcUserId=data[[2]], FromNode=data[[3]], ToNode=data[[4]], LeftPoly=data[[5]], RightPoly=data[[6]], NVertices=data[[7]])

	list(df, data[[8]])
}


get.bnddata <- function(infodir, tablename) 
	.Call("get_bnd_data", as.character(infodir), as.character(tablename), PACKAGE="RArcInfo")

get.paldata <- function(datadir, coverage, filename="pal.adf") 
{
	data<-.Call("get_pal_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")

	#a table (dataframe) with the first six fields is built
	df<-data.frame(PolygonId=data[[1]], MinX=data[[2]], MinY=data[[3]], MaxX=data[[4]], MaxY=data[[5]], NArcs=data[[6]])

	list(df, data[[7]])
}

get.labdata <- function(datadir, coverage, filename="lab.adf") 
{
	data<-.Call("get_lab_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")
	data.frame(LabelUserID=data[[1]], PolygonID=data[[2]], Coord1X=data[[3]], Coord1Y=data[[4]], Coord2X=data[[5]], Coord2Y=data[[6]], Coord3X=data[[7]], Coord3Y=data[[8]])
}

get.cntdata <- function(datadir, coverage, filename="cnt.adf") 
{
	data<-.Call("get_cnt_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")

	df<-data.frame(PolygonID=data[[1]], CoordX=data[[2]], CoordY=data[[3]], NLabels=data[[4]])

	list(df, data[[5]])
}

get.toldata <- function(datadir, coverage, filename="tol.adf") 
{
	data<-.Call("get_tol_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")
	data.frame(Type=data[[1]], Status=data[[2]], Value=data[[3]])
}


get.txtdata <- function(datadir, coverage, filename="txt.adf") 
{
	data<-.Call("get_txt_data", as.character(datadir), as.character(coverage), as.character(filename), PACKAGE="RArcInfo")

	df<-data.frame(TxtID=data[[1]], UserId=data[[2]], Level=data[[3]], NVerticesLine=data[[4]], NVerticesArrow=data[[5]], Text=data[[6]])

	list(df, data[[7]])
}


get.tabledata <- function(infodir, tablename) 
{
	data<-.Call("get_table_data", as.character(infodir), as.character(tablename), PACKAGE="RArcInfo")

	df<-data.frame(I(data[[1]]))
	l<-length(data)

	if(l>=2)
	{
		for (i in 2:l)
			df<-cbind(df,I(data[[i]]))
	}

	fields<-get.tablefields(infodir, tablename)
	names(data)<-fields[[1]]

	data
}

e00toavc <- function(e00file, avcdir)
{
	.Call("e00toavc", as.character(e00file), as.character(avcdir), PACKAGE="RArcInfo")
}

avctoe00 <- function(avcdir, e00file)
{
	.Call("avctoe00", as.character(avcdir), as.character(e00file), PACKAGE="RArcInfo") 
}
