
## Limit a bond data set to a certain maturity range 
maturity_range <- function(bonddata,lower,upper) {
  m <- lapply(bonddata,create_maturities_matrix)
  colmax <- function(m) apply(m,2,max)
  m_max <- lapply(m,colmax)

  bonds_in_range <- function(group) names(group[which(
                group>lower & group<upper)])

  # list with ISINs of bonds in range
  isins_range <- lapply(m_max,bonds_in_range)
  index_set <- which(unlist(lapply(isins_range,length))>0)

  # list with positions of bonds in isins_range
  isins_range_pos <- list()

  for (i in seq_along(bonddata)) {
        isins_range_pos[[i]] <- which(bonddata[[i]]$ISIN %in% 
            isins_range[[i]])
  }
  names(isins_range_pos) <- names(bonddata)

  # first part of bonddata for filtering

  N <- which(names(bonddata[[1]]) %in% c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE","PRICE","ACCRUED"))

  first <- function(lst) lst[N]
  filtered <- lapply(bonddata,first)

  bonddata_range <- list()
    for (i in seq_along(bonddata)) {
      bonddata_range[[i]] <- as.list(as.data.frame(filtered[[i]])
                       [isins_range_pos[[i]],])
      # convert to character                       
      bonddata_range[[i]][["ISIN"]] <- as.character(bonddata_range[[i]][["ISIN"]])
      bonddata_range[[i]][["MATURITYDATE"]] <- as.character(bonddata_range[[i]][["MATURITYDATE"]])
      bonddata_range[[i]][["STARTDATE"]] <- as.character(bonddata_range[[i]][["STARTDATE"]])
    }
  names(bonddata_range) <- names(bonddata)

  # list with positions of cashflows in isins_range

  isins_range_pos <- list()
    for (i in seq_along(bonddata)) {
     isins_range_pos[[i]] <- which(bonddata[[i]][["CASHFLOWS"]]
                        [["ISIN"]]%in%isins_range[[i]])
    }
  names(isins_range_pos) <- names(bonddata)

  for (i in seq_along(bonddata)) {
    CASHFLOWS <- as.list(as.data.frame(bonddata[[i]]
             [["CASHFLOWS"]])[isins_range_pos[[i]],])
    CASHFLOWS$ISIN <- as.character(CASHFLOWS$ISIN)
    CASHFLOWS$DATE <- as.character(CASHFLOWS$DATE)             
    bonddata_range[[i]][["CASHFLOWS"]] <- list()
    bonddata_range[[i]][["CASHFLOWS"]] <- CASHFLOWS
  }
  names(bonddata_range) <- names(bonddata)
  bonddata_range

  # add TODAY from bonddata
  for (i in seq_along(bonddata)) {
    bonddata_range[[i]][["TODAY"]] <- bonddata[[i]][["TODAY"]]
  }

  # delete countries where no bonds are available
  bonddata_range <- bonddata_range[index_set]
  bonddata_range
  }

rm_bond <- function(bonddata, group, ISIN) UseMethod("rm_bond")

## remove bonds from a static bonddata set  
rm_bond.couponbonds <- function(bonddata, group, ISIN){
        cf_isin_index <- which(bonddata[[group]]$CASHFLOWS$ISIN %in% ISIN)
 	isin_index <- which(bonddata[[group]]$ISIN %in% ISIN)	
    	bonddata[[group]]$ISIN <-  bonddata[[group]]$ISIN[-isin_index]
    	bonddata[[group]]$MATURITYDATE <- bonddata[[group]]$MATURITYDATE[-isin_index]
    	bonddata[[group]]$ISSUEDATE <- bonddata[[group]]$ISSUEDATE[-isin_index]
    	bonddata[[group]]$COUPONRATE <- bonddata[[group]]$COUPONRATE[-isin_index]
    	bonddata[[group]]$PRICE <- bonddata[[group]]$PRICE[-isin_index]
    	bonddata[[group]]$ACCRUED <- bonddata[[group]]$ACCRUED[-isin_index]
        bonddata[[group]]$CASHFLOWS$ISIN <- bonddata[[group]]$CASHFLOWS$ISIN[-cf_isin_index]
        bonddata[[group]]$CASHFLOWS$CF <- bonddata[[group]]$CASHFLOWS$CF[-cf_isin_index]
        bonddata[[group]]$CASHFLOWS$DATE <- bonddata[[group]]$CASHFLOWS$DATE[-cf_isin_index]
     
	
	class(bonddata) <- "couponbonds"
        bonddata
}

## remove bonds from a dynamic bonddata set 
rm_bond.dyncouponbonds <- function(bonddata, group, ISIN) {
  for (i in seq(length(bonddata))) {
     cf_isin_index <- which(bonddata[[i]][[group]]$CASHFLOWS$ISIN %in% ISIN)
     isin_index <- which(bonddata[[i]][[group]]$ISIN %in% ISIN)
     bonddata[[i]][[group]]$ISIN <- bonddata[[i]][[group]]$ISIN[-isin_index]
     bonddata[[i]][[group]]$MATURITYDATE <- bonddata[[i]][[group]]$MATURITYDATE[-isin_index]
     bonddata[[i]][[group]]$STARTDATE <- bonddata[[i]][[group]]$STARTDATE[-isin_index]
     bonddata[[i]][[group]]$COUPONRATE <- bonddata[[i]][[group]]$COUPONRATE[-isin_index]
     bonddata[[i]][[group]]$PRICE <- bonddata[[i]][[group]]$PRICE[-isin_index]
     bonddata[[i]][[group]]$ACCRUED <- bonddata[[i]][[group]]$ACCRUED[-isin_index]
     bonddata[[i]][[group]]$CASHFLOWS$ISIN <- bonddata[[i]][[group]]$CASHFLOWS$ISIN[-cf_isin_index]
     bonddata[[i]][[group]]$CASHFLOWS$CF <- bonddata[[i]][[group]]$CASHFLOWS$CF[-cf_isin_index]
     bonddata[[i]][[group]]$CASHFLOWS$DATE <- bonddata[[i]][[group]]$CASHFLOWS$DATE[-cf_isin_index]
  }
     class(bonddata) <- "dyncouponbonds"
     bonddata
}



## preprocess a static bonddataset,i.e., sort data, calculate cashflows, maturity matrix, yields 
prepro_bond <- function(group,
           bonddata,
           matrange="all"){

  # select given group from bonddata
  bonddata <- bonddata[group]

  # select data according to chosen maturity range
  if (length(matrange)==1) {bonddata <- bonddata }else
   {bonddata <- maturity_range(bonddata,matrange[1],matrange[2]) }

  # number of groups 
  n_group <- length(bonddata) 
  
  # group sequence
  sgroup <- seq(n_group)
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE),
                 sgroup,SIMPLIFY=FALSE)
  
  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE),
                sgroup,SIMPLIFY=FALSE)
  
  # calculate dirty prices
  p <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)

  
  # extract accrued interest
  ac <- mapply(function(k) bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)

  # assign ISIN 
  for(k in sgroup) {
                    names(p[[k]]) <- bonddata[[k]]$ISIN
                    names(ac[[k]]) <- bonddata[[k]]$ISIN
                 
                  }
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  p <- mapply(function(k) p[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  ac <- mapply(function(k) ac[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
 
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]]),
                   sgroup,SIMPLIFY=FALSE)
  # calculate duration   
  duration <- mapply(function(k) duration(cf_p[[k]],m_p[[k]],y[[k]][,2]),
                   sgroup,SIMPLIFY=FALSE)

  res <- list(n_group=n_group,sgroup=sgroup,cf=cf,cf_p=cf_p,m=m,m_p=m_p,p=p,ac=ac,y=y,duration=duration,timestamp=bonddata[[1]]$TODAY)
  res
}

## Bonddata postprocessing function                      

postpro_bond <- function(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda){
  
 # theoretical bond prices with estimated parameters
 phat <- mapply(function(k) bond_prices(method,opt_result[[k]]$par,
       m[[k]],cf[[k]],lambda)$bond_prices,sgroup,SIMPLIFY=FALSE)

  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
 
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"

  
  # maturity interval
  t <- seq(round(min(mapply(function(i) min(y[[i]][,1]), sgroup)),2),
                           ceiling(max(mapply(function(i) max(y[[i]][,1]), sgroup))),0.01)
  
  
  # calculate zero coupon yield curves  
  zcy_curves <- mapply(function(k) cbind(t,spotrates(method,opt_result[[k]]$par,t,lambda)/100),sgroup,SIMPLIFY=FALSE)
 
  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
                                      
  # calculate spread curves              	 
   if(n_group != 1) {  
   s_curves <- mapply(function(k) cbind(t,(zcy_curves[[k]][,2] - zcy_curves[[1]][,2])),sgroup,
   					SIMPLIFY=FALSE)
    } else s_curves = "none"
   
   for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
   class(s_curves) <- "s_curves"
    
  # calculate extrapolation point                        
  expoints <- mapply(function(k) which(zcy_curves[[k]][,1] > 
                 mapply(function(i) max(y[[i]][,1]), seq(n_group))[k])[1],sgroup, SIMPLIFY=FALSE )  
        
  # calculate forward rate curves 
  fwr_curves <-  mapply(function(k) cbind(t,forwardrates(method,opt_result[[k]]$par,t,lambda)/100),sgroup,SIMPLIFY=FALSE)                   
                                        
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"

  
  # calculate discount factor curves 
  df_curves <- mapply(function(k) cbind(zcy_curves[[k]][,1],exp(-zcy_curves[[k]][,1]*
                                        zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
   for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
   class(df_curves) <- "df_curves"

  res <- list(phat=phat,perrors=perrors,yhat=yhat,yerrors=yerrors,t=t,zcy_curves=zcy_curves,
              s_curves=s_curves,expoints=expoints,fwr_curves=fwr_curves,df_curves=df_curves,opt_result=opt_result)
  res

}
