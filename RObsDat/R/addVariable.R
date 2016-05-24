addVariable <- function(Code, Name, 
			Speciation=rep("Unknown", NROW(Code)), 
			Unit, 
			SampleMedium =rep("Unknown", NROW(Code)),
			ValueType=rep("Unknown", NROW(Code)),
			IsRegular=rep("True", NROW(Code)),
			TimeSupport=rep(0, NROW(Code)),
			TimeUnits=rep("Julian year", NROW(Code)),
			DataType=rep("Unknown", NROW(Code)),
			GeneralCategory=rep("Unknown", NROW(Code)),
			NoDataValue=rep(-999999, NROW(Code))
		       	){

	stopifnot(length(Name) == length(Code))
	stopifnot(length(Speciation) == length(Code))
	stopifnot(length(SampleMedium) == length(Code))
	stopifnot(length(ValueType) == length(Code))
	stopifnot(length(DataType) == length(Code))
	stopifnot(length(GeneralCategory) == length(Code))
	stopifnot(length(Unit) == length(Code))
	stopifnot(length(TimeUnits) == length(Code))
	#Fields with no reference
	stopifnot(length(NoDataValue) == length(Code))
	stopifnot(length(IsRegular) == length(Code))
	stopifnot(length(TimeSupport) == length(Code))


	#check for existing entries
	if(NROW(existing <- getMetadata("Variable",Code=Code, Name=Name, Speciation = Speciation, SampleMedium=SampleMedium, ValueType=ValueType, DataType=DataType, GeneralCategory=GeneralCategory, exact=TRUE))>0){
		warning(paste("Existing Variable entry:", Name, " -- Skiping all imports!"))
		return()
	}
	#Handle CV-Fields: VariableName, Speciation, SampleMedium, ValueType, DataType, GeneralCategory
	VariableNameID <- getID("VariableName",Name)
	SpeciationID <- getID("Speciation",Speciation)
	SampleMediumID <- getID("SampleMedium",SampleMedium)
	ValueTypeID <- getID("ValueType",ValueType)
	DataTypeID <- getID("DataType",DataType)
	GeneralCategoryID <- getID("GeneralCategory",GeneralCategory)
	#Other Fields with foreign key: VariableUnits, TimeUnits
	VariableUnitsID <- getID("Units",Unit)
	TimeUnitsID <- getID("Units",TimeUnits)



	IaddVariable(getOption("odm.handler"), Code=Code, Name=VariableNameID, Speciation=SpeciationID, Unit=VariableUnitsID, SampleMedium=SampleMediumID,ValueType=ValueTypeID, IsRegular=IsRegular, TimeSupport=TimeSupport, TimeUnits=TimeUnitsID, DataType=DataTypeID, GeneralCategory=GeneralCategoryID, NoDataValue=NoDataValue)


}
