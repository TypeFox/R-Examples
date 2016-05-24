#'@name k600.2.kGAS
#'@aliases 
#'k600.2.kGAS
#'k600.2.kGAS.base
#'@title Returns the gas exchange velocity for gas of interest w/ no unit conversions
#'@description 
#'Returns the gas exchange velocity for gas of interest w/ no unit conversions
#'@usage
#'k600.2.kGAS.base(k600,temperature,gas="O2")
#'
#'k600.2.kGAS(ts.data, gas="O2")
#'
#'@param ts.data Object of class data.frame with named columns datetime and k600 and wtr (water temp in deg C). Other columns are ignored
#'@param k600 k600 as vector array of numbers or single number
#'@param temperature Water temperature (deg C) as vector array of numbers or single number
#'@param gas gas for conversion, as string (e.g., 'CO2' or 'O2')
#'@return Numeric value of gas exchange velocity for gas
#'@author Jordan S. Read
#'@seealso \link{k.read} and \link{k.read.base} for functions that calculate k600 estimates
#'@examples
#'  ## single example
#'kO2 <- k600.2.kGAS.base(k600=2.4,temperature=20.4,gas='O2')
#'
#'## Timeseries example
#'#load data
#'data.path = system.file('extdata', package="LakeMetabolizer")
#'sp.data = load.all.data('sparkling', data.path)
#'ts.data = sp.data$data #pull out just the timeseries data
#'
#'
#'#calculate U10 and add it back onto the original 
#'u10 = wind.scale(ts.data)
#'ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
#'ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset  
#'
#'k600 = k.cole(ts.data)
#'ts.data = merge(k600, ts.data)
#'
#'k.gas = k600.2.kGAS(ts.data, 'O2')
#'@export
k600.2.kGAS.base	<-	function(k600,temperature,gas="O2"){
	
	n	<-	0.5
	schmidt	<-	getSchmidt(temperature,gas)
	Sc600	<-	schmidt/600
	
	kGAS	<-	k600*(Sc600^-n)
	return(kGAS)
}

#'@export
k600.2.kGAS = function(ts.data, gas="O2"){
  
	k600 = get.vars(ts.data, 'k600')
  temperature = get.vars(ts.data,'wtr')[,1:2]
  
  kGAS = data.frame(datetime=ts.data$datetime)
  
  kGAS$k.gas = k600.2.kGAS.base(k600[,2], temperature[,2], gas)
  
  return(kGAS)
}
