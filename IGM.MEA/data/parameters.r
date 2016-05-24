parameters <- list(

burst.type = "mi",
s.min=5,
  
perm.n = 1000, ##Might want to specify in parameters or give user option?

elec.min.rate = (1/60),
elec.max.rate = 1000,
well.min.rate = 0,

#parameters for max-interval burst detection
mi.par = list(beg.isi =    0.1,
               end.isi =    0.25,
               min.ibi =    0.8,
               min.durn =   0.05,
               min.spikes = 5),

## Parameters for network.spikes
ns.T = 0.01,    	#time in seconds
ns.N = 3,        #how many coincident electrodes?
sur=100, # num. ms before and after spike to check I think, used in ms

#sahar -added for ver 2.0
# Burst parameters for distribution analysis

# Parameters for inter burst interval distribution analysis
burst.distribution.IBI = list(
  perform = 1,        # 0/1 - run this analysis ?
  min.cases = 900/60, # minimum number of bursts, below which electrode will be ingnored
  # 900 is length of recording divided by 60 sec
  x.lim =         20,  # x max limit for distribution plot
  bins.in.seg =  5,   # how many values to be calculated in each segment of xLim 
  # (overall values will be bins.in.seg * xlim)
  min.values = 0,     # bursts with values below this threshold be ignored with filter.by.min=1
  filter.by.min = 0,  # 0/1 ignore bursts with values below min.values
  per.well = 0),       # 0/1 - perform analysis per well=1 or per electrode=0

# Parameters for burst duration distribution analysis
burst.distribution.durn = list(
  perform = 1,        # 0/1 - run this analysis ?
  min.cases = 900/60, # minimum number of bursts, below which electrode will be ingnored
  # 900 is length of recording divided by 60 sec
  x.lim =         18,  # x max limit for distribution plot
  bins.in.seg =  10,   # how many values to be calculated in each segment of xLim 
  # (overall values will be bins.in.seg * xlim)
  min.values = 0,     # bursts with values below this threshold be ignored with filter.by.min=1
  filter.by.min = 0,  # 0/1 ignore bursts with values below min.values
  per.well = 0),      # 0/1 - perform analysis per well=1 or per electrode=0

# Parameters for inter spike interval within burst distribution analysis
burst.distribution.ISI = list(
  perform = 1,        # 0/1 - run this analysis ?
  min.cases = 900/60, # minimum number of bursts, below which electrode will be ingnored
  x.lim =         0.5,  # x max limit for distribution plot
  bins.in.seg =  100,   # how many values to be calculated in each segment of xLim 
  # (overall values will be bins.in.seg * xlim)
  min.values = 0,     # bursts with values below this threshold be ignored with filter.by.min=1
  filter.by.min = 0,  # 0/1 ignore bursts with values below min.values
  per.well = 0),      # 0/1 - perform analysis per well=1 or per electrode=0

# Parameters for number of spikes in burst distribution analysis
burst.distribution.nSpikes = list(
  perform = 1,        # 0/1 - run this analysis ?
  min.cases = 5, # minimum number of bursts, below which electrode will be ingnored
  x.lim =         200,  # x max limit for distribution plot
  bins.in.seg =  1,   # how many values to be calculated in each segment of xLim 
  # (overall values will be bins.in.seg * xlim)
  min.values = 0,     # bursts with values below this threshold be ignored with filter.by.min=1
  filter.by.min = 0,  # 0/1 ignore bursts with values below min.values
  per.well = 0),      # 0/1 - perform analysis per well=1 or per electrode=0

# Parameters for average spike frequency in burst distribution analysis
burst.distribution.spikeFreq = list(
  perform = 1,        # 0/1 - run this analysis ?
  min.cases = 900/60, # minimum number of bursts, below which electrode will be ingnored
  x.lim =         300,  # x max limit for distribution plot
  bins.in.seg =  1,   # how many values to be calculated in each segment of xLim 
  # (overall values will be bins.in.seg * xlim)
  min.values = 0,     # bursts with values below this threshold be ignored with filter.by.min=1
  filter.by.min = 0,  # 0/1 ignore bursts with values below min.values
  per.well = 0),     # 0/1 - perform analysis per well=1 or per electrode=0

  #network burst parameters
  local_region_min_nAE= 0, #do not change for now
  min_electrodes= 4 ,
  Sigma= c(10,20,50) # a list of window size to be considered 

) # end of parameters list


