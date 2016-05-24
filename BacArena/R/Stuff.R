globalVariables(c("Ec_core"))


#color dictionary of 269 maximally distinct colors from all previous colors 
colpal1 <- c("#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B")

# 20 optimally distinct colors
colpal2 <- c("#C48736", "#CE54D1", "#96CED5", "#76D73C", "#403552", "#D4477D", "#5A7E36", "#D19EC4", "#CBC594", "#722A2D", "#D0CD47", "#CF4A31", "#7B6FD0", "#597873", "#6CD3A7", "#484125", "#C17E73", "#688EC1",  "#844081", "#7DD06F")

# K. Kelly (1965): Twenty-two colors of maximum contrast. // Color Eng., 3(6)
colpal3 = c(
  "#FFB300", # Vivid Yellow
  "#803E75", # Strong Purple
  "#FF6800", # Vivid Orange
  "#A6BDD7", # Very Light Blue
  "#C10020", # Vivid Red
  "#CEA262", # Grayish Yellow
  "#817066", # Medium Gray
  "#007D34", # Vivid Green
  "#F6768E", # Strong Purplish Pink
  "#00538A", # Strong Blue
  "#FF7A5C", # Strong Yellowish Pink
  "#53377A", # Strong Violet
  "#FF8E00", # Vivid Orange Yellow
  "#B32851", # Strong Purplish Red
  "#F4C800", # Vivid Greenish Yellow
  "#7F180D", # Strong Reddish Brown
  "#93AA00", # Vivid Yellowish Green
  "#593315", # Deep Yellowish Brown
  "#F13A13", # Vivid Reddish Orange
  "#232C16" # Dark Olive Green
)

# Diffusion pde solver function
Diff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC
    return (list(dCONC))
  })
}

# advection + diffusion + conservation
AdvecDiffConserv2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    vgrid <- setup.prop.2D(value = 1, y.value=0, grid = gridgeometry.grid2D)
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     v.grid = vgrid, flux.x.down=0, flux.x.up=0, flux.y.up=0, flux.y.down=0)$dC
    return (list(dCONC))
  })
}

AdvecDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    vgrid <- setup.prop.2D(value = 1, y.value=0, grid = gridgeometry.grid2D)
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    #dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, v.grid = vgrid, C.x.down=rep(0.1, gridgeometry.grid2D$x.N))$dC
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     v.grid = vgrid, 
                     C.y.down=rep(1, gridgeometry.grid2D$x.N),
                     flux.x.down=1, flux.x.up=-1, flux.y.up=0, flux.y.down=0)$dC
    return (list(dCONC))
  })
}

ConstBoundDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     C.y.down=rep(boundS, gridgeometry.grid2D$y.N),
                     C.y.up=rep(boundS, gridgeometry.grid2D$y.N),
                     C.x.down=rep(boundS, gridgeometry.grid2D$x.N),
                     C.x.up=rep(boundS, gridgeometry.grid2D$x.N))$dC
    return (list(dCONC))
  })
}

InfluxBoundDiff2d <- function (t, y, parms)  {
  # geometry values are in parms
  with (as.list(parms), {
    CONC  <- matrix(nrow = gridgeometry.grid2D$x.N, ncol = gridgeometry.grid2D$y.N, data = y)
    dCONC <- tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid, 
                     flux.y.down=rep(-boundS, gridgeometry.grid2D$y.N),
                     flux.y.up=rep(boundS, gridgeometry.grid2D$y.N),
                     flux.x.down=c(0,rep(-boundS, gridgeometry.grid2D$x.N-2),0),  # reduce edge bias
                     flux.x.up=c(0, rep(boundS, gridgeometry.grid2D$x.N-2),0))$dC # "
    return (list(dCONC))
  })
}





# function estimates lrw array size paremeter needed to solve stiffy diffusion pde with solver lsodes
estimate_lrw <- function(grid_n, grid_m){
  x=c(10*10, 25*25, 51*51, 61*61, 71*71, 81*81, 91*91, 101*101)
  y=c(3901, 29911, 160000, 230000, 330000, 430000, 580000, 710000)
  lm <- lm(y~x)
  #summary(lm)
  #plot(x,y)
  #abline(coef(lm))
  #abline(coef=c(0, lm$coefficients[2]))
  lrw <- as.numeric(lm$coefficients[2]*grid_n*grid_m + grid_n*grid_m*100)
  #lrw <- ((grid_n*grid_m)*18.5 + 20)*10 -> alternative function
  return(lrw)
}

#' @title Start simulation
#'
#' @description The function \code{openArena} can be used to start a default simulation.
#' @export
#' @rdname openArena
#' @importFrom utils data
#'
#' @return Returns an object of class \code{Eval} which can be used for subsequent analysis steps.
#' @examples
#' \donttest{ 
#' sim <- openArena()
#' evalArena(sim, time=7, phencol = TRUE, 
#'           plot_items=c("Population", "EX_o2(e)", "EX_for(e)",
#'           "EX_glc(e)", "EX_for(e)"))
#'}
openArena <- function(){
  data(Ec_core, envir = environment())
  bac = Bac(model=Ec_core, type="E. coli")
  arena = Arena(n=50, m=50, stir=F, Lx=0.0125, Ly=0.0125)
  addOrg(arena, bac, amount=50)
  addSubs(arena, smax=0, difspeed=6.7e-6)
  changeSub(arena, smax=0.05, unit="mM", 
            mediac = c("EX_glc(e)", "EX_o2(e)", "EX_h2o(e)", "EX_pi(e)", "EX_nh4(e)", "EX_h(e)"))
  sim <- simEnv(arena, time=10)  
  plotCurves2(sim, legendpos="left")
  return(sim)
}

reset_screen <- function(){
  par(mfrow=c(1,1))
}
