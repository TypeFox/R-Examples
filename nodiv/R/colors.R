

int.color <- function(col, n, alpha){
  if (n == length(col) & alpha == 1) 
    return(col)
  rgb.col <- t(col2rgb(col))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, , length(col))
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    hold <- spline(x, rgb.col[, k], xout = xg)$y
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  if (alpha == 1) {
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  }
  else {
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255, 
        alpha = alpha)
  }
}


choose.colors <- function(vec, zlim = NULL, coltype = c("auto", "ramp", "monochrome", "divergent", "individual"))
{
  
  coltype = match.arg(coltype)
  
  if(is.character(vec))
    vec <- as.factor(vec)
  
  if(is.factor(vec))
    zlim = c(1,length(levels(vec)))
  
  if(is.null(zlim)){
    zlim <- range(vec, na.rm = TRUE)
    
    if(zlim[1] * zlim[2] < 0)
      zlim <- c(-max(abs(vec), na.rm = TRUE), max(abs(vec), na.rm = TRUE))
  }
  
  if(coltype == "auto"){
    if(is.factor(vec))
      coltype <- "individual" else {
        
        if(zlim[1] * zlim[2] < 0){
          if(sum(zlim) == 0){
            coltype <- "divergent" 
          } else {
            warning("ramp color scheme chosen - if you want divergent colors use zlims symmetric around 0")
            coltype <- "ramp"
          }
          
        } else coltype <- "ramp" 
      }
  }
  
  if(coltype == "individual"){
    n <- length(unique(vec))
    if(requireNamespace("RColorBrewer")){
      if(n <= 9)
        ret <- RColorBrewer::brewer.pal(n, "Set1") else
          if(n <= 12)
            ret <- RColorBrewer::brewer.pal(n, "Set3") else
              if(n <= 21)
                ret <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"))[1:n] else 
                  if(requireNamespace("colorspace"))
                    ret <- colorspace::rainbow_hcl(n) else
                      ret <- rainbow(n)
    } else {
      if(n <= 8)
        ret <- palette()[1:n] else
        {
          if(requireNamespace("colorspace"))
            ret <- colorspace::rainbow_hcl(n) else 
              ret <- rainbow(n)
        }
      
    }
  }
  
  if(coltype == "ramp"){
    ret <- parula() 
  # ret <- HMrainbow 
    #if(requireNamespace("fields"))
    #  ret <- fields::tim.colors(64) else
    #    ret <- rev(terrain.colors(64))
  }
  
  if(coltype == "divergent"){
    #if(requireNamespace("RColorBrewer"))
      #ret <- RColorBrewer::brewer.pal(11, "RdBu") else
        #if(requireNamespace("colorspace"))
          #ret <- colorspace::diverge_hcl(64) else
            ret <- redblue(64)
  }
  
  if(coltype == "monochrome"){
    ret <- blackbodycol()
    #if(requireNamespace("RColorBrewer"))
    #  ret <- RColorBrewer::brewer.pal(9, "YlOrRd") else
    #    if(requireNamespace("colorspace"))
    #      ret <- colorspace::sequential_hcl(64) else
    #        ret <- rev(heat.colors(64))
  }
  
  if(length(unique(na.omit(vec))) == 2)
    ret <- c("darkgrey", "red")
  
  attr(ret, "zlim") <- zlim
  
  return(ret)
}

custom_palette <- function(colname = c("parula", "jet", "blackbody", "HMblueyellow", 
                                       "HMrainbow", "HMlinear_optimal", "HMoptimal_scale", 
                                       "cube1", "cubeyf1", "redblue", "moreland"), n = NULL, alpha = 1){
  colname <- match.arg(colname)
  if(is.null(n)){
    if(colname %in% c("parula", "jet"))
      n <- 64 else 
        if(colname %in% c("RedBlue")) n <- 12 else 
          if(colname %in% c("moreland")) n <- 33 else 
            n <- 256
  }
  switch(colname,
         parula = parula(n, alpha),
         jet = jetcol(n, alpha),
         blackbody = blackbodycol(n, alpha),
         HMblueyellow = bluetoyellow(n, alpha),
         HMrainbow = HMrainbow(n, alpha),
         HMlinear_optimal = Lin_optimcol(n, alpha),
         HMoptimal_scale = optim_scalecol(n, alpha),
         cube1 = cube1col(n, alpha),
         cubeyf1 = cubeyf1col(n, alpha),
         redblue = redblue(n, alpha),
         moreland = moreland(n, alpha)
         )
}

create.cols <- function(vec, col, zlim, coltype = c("auto", "ramp", "monochrome", "divergent", "individual"))
{
  coltype = match.arg(coltype)
  
  if(missing(zlim)){
    zlim <- range(vec, na.rm = TRUE)
    if(zlim[1] * zlim[2] < 0)
      zlim <- c(-max(abs(vec), na.rm = TRUE), max(abs(vec), na.rm = TRUE))
  }
  
  if(missing(col)) 
    col <- choose.colors(vec, zlim, coltype)
  
  if(length(col) == length(unique(vec)))
    return(col[match(vec, sort(unique(vec)))])
  
  if(min(vec, na.rm = TRUE) == 0 & identical(vec, floor(vec)) & length(col) > 3) col <- c("grey", col)
  
  vec = vec - zlim[1]
  vec = floor(vec * (length(col)-1)/(zlim[2]- zlim[1]))+1
  return(col[vec])
}

redblue <- function(n = 12, alpha = 1){
  orig <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", 
            "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", 
            "#053061")
  int.color(rev(orig), n, alpha)
}

moreland <- function(n = 33, alpha = 1){
  orig <- c("#3B4CC0", "#445ACC", "#4D68D7", "#5775E1", "#6282EA", 
            "#6C8EF1", "#779AF7", "#82A5FB", "#8DB0FE", "#98B9FF", 
            "#A3C2FF", "#AEC9FD", "#B8D0F9", "#C2D5F4", "#CCD9EE", 
            "#D5DBE6", "#DDDDDD", "#E5D8D1", "#ECD3C5", "#F1CCB9", 
            "#F5C4AD", "#F7BBA0", "#F7B194", "#F7A687", "#F49A7B", 
            "#F18D6F", "#EC7F63", "#E57058", "#DE604D", "#D55042", 
            "#CB3E38", "#C0282F", "#B40426")
  int.color(orig, n, alpha)
}

parula <- function(n = 64, alpha = 1){
  orig <- c("#352A87", "#363093", "#3637A0", "#353DAD", "#3243BA", 
            "#2C4AC7", "#2053D4", "#0F5CDD", "#0363E1", "#0268E1", 
            "#046DE0", "#0871DE", "#0D75DC", "#1079DA", "#127DD8", 
            "#1481D6", "#1485D4", "#1389D3", "#108ED2", "#0C93D2", 
            "#0998D1", "#079CCF", "#06A0CD", "#06A4CA", "#06A7C6", 
            "#07A9C2", "#0AACBE", "#0FAEB9", "#15B1B4", "#1DB3AF", 
            "#25B5A9", "#2EB7A4", "#38B99E", "#42BB98", "#4DBC92", 
            "#59BD8C", "#65BE86", "#71BF80", "#7CBF7B", "#87BF77", 
            "#92BF73", "#9CBF6F", "#A5BE6B", "#AEBE67", "#B7BD64", 
            "#C0BC60", "#C8BC5D", "#D1BB59", "#D9BA56", "#E1B952", 
            "#E9B94E", "#F1B94A", "#F8BB44", "#FDBE3D", "#FFC337", 
            "#FEC832", "#FCCE2E", "#FAD32A", "#F7D826", "#F5DE21", 
            "#F5E41D", "#F5EB18", "#F6F313", "#F9FB0E")
  int.color(orig, n, alpha)
}

bluetoyellow <- function(n = 256, alpha = 1){
  orig <- c("#0707FE", "#1717FC", "#1E1EFA", "#2424F8", "#2828F7",
            "#2C2CF5", "#2F2FF3", "#3232F2", "#3434F0", "#3737EF",
            "#3939EE", "#3B3BEC", "#3D3DEB", "#3F3FEA", "#4141E9",
            "#4242E7", "#4444E6", "#4545E5", "#4747E4", "#4848E3", 
            "#4A4AE2", "#4B4BE1", "#4C4CE1", "#4E4EE0", "#4F4FDF",
            "#5050DE", "#5151DD", "#5252DD", "#5454DC", "#5555DB",
            "#5656DA", "#5757DA", "#5858D9", "#5959D8", "#5A5AD8", 
            "#5B5BD7", "#5C5CD6", "#5D5DD6", "#5E5ED5", "#5F5FD5", 
            "#6060D4", "#6161D4", "#6262D3", "#6262D2", "#6363D2", 
            "#6464D1", "#6565D1", "#6666D0", "#6767D0", "#6868D0", 
            "#6969CF", "#6969CF", "#6A6ACE", "#6B6BCE", "#6C6CCD", 
            "#6D6DCD", "#6E6ECC", "#6E6ECC", "#6F6FCC", "#7070CB", 
            "#7171CB", "#7272CA", "#7272CA", "#7373CA", "#7474C9", 
            "#7575C9", "#7676C8", "#7676C8", "#7777C8", "#7878C7", 
            "#7979C7", "#7979C7", "#7A7AC6", "#7B7BC6", "#7C7CC6", 
            "#7C7CC5", "#7D7DC5", "#7E7EC5", "#7F7FC4", "#8080C4", 
            "#8080C3", "#8181C3", "#8282C3", "#8282C2", "#8383C2", 
            "#8484C2", "#8585C1", "#8585C1", "#8686C1", "#8787C0", 
            "#8888C0", "#8888C0", "#8989BF", "#8A8ABF", "#8B8BBF", 
            "#8B8BBE", "#8C8CBE", "#8D8DBE", "#8E8EBD", "#8E8EBD", 
            "#8F8FBD", "#9090BC", "#9090BC", "#9191BC", "#9292BB", 
            "#9393BB", "#9393BB", "#9494BA", "#9595BA", "#9595BA", 
            "#9696B9", "#9797B9", "#9898B9", "#9898B8", "#9999B8", 
            "#9A9AB8", "#9A9AB7", "#9B9BB7", "#9C9CB6", "#9D9DB6", 
            "#9D9DB6", "#9E9EB5", "#9F9FB5", "#9F9FB5", "#A0A0B4", 
            "#A1A1B4", "#A2A2B4", "#A2A2B3", "#A3A3B3", "#A4A4B2", 
            "#A4A4B2", "#A5A5B2", "#A6A6B1", "#A7A7B1", "#A7A7B0", 
            "#A8A8B0", "#A9A9B0", "#A9A9AF", "#AAAAAF", "#ABABAE", 
            "#ACACAE", "#ACACAD", "#ADADAD", "#AEAEAD", "#AEAEAC", 
            "#AFAFAC", "#B0B0AB", "#B1B1AB", "#B1B1AA", "#B2B2AA", 
            "#B3B3A9", "#B3B3A9", "#B4B4A8", "#B5B5A8", "#B5B5A7", 
            "#B6B6A7", "#B7B7A6", "#B8B8A6", "#B8B8A5", "#B9B9A5", 
            "#BABAA4", "#BABAA4", "#BBBBA3", "#BCBCA3", "#BDBDA2", 
            "#BDBDA2", "#BEBEA1", "#BFBFA1", "#BFBFA0", "#C0C09F", 
            "#C1C19F", "#C2C29E", "#C2C29E", "#C3C39D", "#C4C49D",
            "#C4C49C", "#C5C59B", "#C6C69B", "#C7C79A", "#C7C799",
            "#C8C899", "#C9C998", "#C9C997", "#CACA97", "#CBCB96", 
            "#CCCC95", "#CCCC95", "#CDCD94", "#CECE93", "#CECE92", 
            "#CFCF92", "#D0D091", "#D1D190", "#D1D18F", "#D2D28F", 
            "#D3D38E", "#D3D38D", "#D4D48C", "#D5D58B", "#D6D68A", 
            "#D6D68A", "#D7D789", "#D8D888", "#D8D887", "#D9D986", 
            "#DADA85", "#DBDB84", "#DBDB83", "#DCDC82", "#DDDD81", 
            "#DDDD80", "#DEDE7F", "#DFDF7E", "#E0E07D", "#E0E07C",
            "#E1E17B", "#E2E27A", "#E2E279", "#E3E377", "#E4E476", 
            "#E5E575", "#E5E574", "#E6E672", "#E7E771", "#E8E870", 
            "#E8E86E", "#E9E96D", "#EAEA6B", "#EAEA6A", "#EBEB68", 
            "#ECEC67", "#EDED65", "#EDED64", "#EEEE62", "#EFEF60", 
            "#EFEF5E", "#F0F05C", "#F1F15B", "#F2F259", "#F2F256", 
            "#F3F354", "#F4F452", "#F5F550", "#F5F54D", "#F6F64A", 
            "#F7F748", "#F7F745", "#F8F841", "#F9F93E", "#FAFA3A", 
            "#FAFA36", "#FBFB31", "#FCFC2C", "#FDFD25", "#FDFD1C", 
            "#FEFE0D")
    int.color(orig, n, alpha)
}
  
jetcol <- function(n = 64, alpha = 1){
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
            "#AF0000", "#9F0000", "#8F0000", "#800000")
    int.color(orig, n, alpha)
}

cube1col <- function(n = 255, alpha = 1){
  orig <- c("#740081", "#760085", "#770088", "#78008B", "#79008E",
            "#7A0091", "#7B0094", "#7C0097", "#7D009A", "#7E009D", 
            "#7F01A0", "#8002A3", "#8004A6", "#8106A9", "#8208AC", 
            "#830BB0", "#830DB3", "#8410B6", "#8512B9", "#8514BC", 
            "#8517BF", "#851AC2", "#851DC5", "#8520C8", "#8523CB", 
            "#8526CD", "#8429D0", "#832CD2", "#832FD4", "#8232D6", 
            "#8135D8", "#8038DA", "#803BDC", "#7F3EDE", "#7E40E0", 
            "#7E42E2", "#7D45E4", "#7C47E7", "#7C49E9", "#7B4BEB", 
            "#794DED", "#784FF0", "#7751F2", "#7653F4", "#7555F6", 
            "#7457F7", "#7359F9", "#725BFA", "#715DFC", "#705EFC", 
            "#6F60FD", "#6E62FD", "#6D64FD", "#6C66FD", "#6B68FD", 
            "#6A6BFC", "#696DFC", "#686FFB", "#6771FB", "#6673FA", 
            "#6675F9", "#6577F8", "#6479F7", "#637BF7", "#627DF6", 
            "#617EF5", "#6080F4", "#5F82F3", "#5E84F2", "#5D86F1", 
            "#5C87F0", "#5B89EF", "#5A8AEE", "#598CEC", "#588EEB", 
            "#578FEA", "#5691E8", "#5592E6", "#5394E5", "#5295E3", 
            "#5197E2", "#5098E0", "#4F99DE", "#4E9BDD", "#4D9CDB", 
            "#4C9ED9", "#4A9FD7", "#49A1D6", "#48A2D4", "#46A4D2", 
            "#45A5D0", "#43A7CE", "#42A8CC", "#40AACA", "#3FABC9", 
            "#3DADC7", "#3CAEC5", "#3BAFC2", "#3AB1C0", "#39B2BE", 
            "#39B3BC", "#38B5BA", "#38B6B8", "#38B7B6", "#38B8B4", 
            "#39B9B2", "#39BBB0", "#3ABCAE", "#3ABDAB", "#3BBEA9", 
            "#3CBFA7", "#3DC0A5", "#3EC1A3", "#3EC2A0", "#3FC39E", 
            "#40C49C", "#41C59A", "#42C697", "#43C795", "#43C893", 
            "#44C991", "#45CA8E", "#45CB8C", "#46CB8A", "#47CC87",
            "#47CD85", "#48CE83", "#49CF81", "#49D07E", "#4AD07C",
            "#4BD17A", "#4BD277", "#4CD375", "#4DD373", "#4DD471", 
            "#4ED56F", "#4FD66C", "#50D76A", "#50D868", "#51D865", 
            "#52D962", "#52DA60", "#53DB5D", "#54DC5A", "#54DD57", 
            "#55DE55", "#56DE52", "#57DF50", "#57E04E", "#58E14C", 
            "#59E14B", "#5AE24A", "#5BE349", "#5CE349", "#5EE449", 
            "#5FE549", "#61E549", "#63E64A", "#65E64A", "#68E74A", 
            "#6AE74B", "#6DE84B", "#6FE84C", "#72E94C", "#75E94D", 
            "#78EA4E", "#7AEA4E", "#7DEA4F", "#80EB4F", "#82EB50", 
            "#85EB50", "#87EB50", "#89EB51", "#8CEB51", "#8EEB52", 
            "#91EB52", "#93EB52", "#96EC53", "#98EC53", "#9BEC54", 
            "#9DEC54", "#A0EC54", "#A2EC55", "#A5EC55", "#A7EC55", 
            "#A9EC56", "#ABEC56", "#ADEC56", "#AFEC57", "#B1EC57", 
            "#B4EC57", "#B6EC57", "#B8EC58", "#B9EC58", "#BBEC58", 
            "#BDEC58", "#BFEC59", "#C1EC59", "#C3EC59", "#C4EC59", 
            "#C6EC59", "#C8EC59", "#C9EC5A", "#CBEC5A", "#CCEC5A", 
            "#CDEC5A", "#CFEC5A", "#D0EB5A", "#D1EA5B", "#D2EA5B", 
            "#D3E95B", "#D4E85B", "#D5E65B", "#D6E55B", "#D7E45B", 
            "#D8E25B", "#D9E15B", "#DAE05C", "#DBDE5C", "#DCDD5C", 
            "#DDDB5C", "#DEDA5C", "#DFD95C", "#E0D75C", "#E2D65C", 
            "#E3D55D", "#E5D35D", "#E6D25D", "#E7D05D", "#E9CE5D", 
            "#EACD5D", "#ECCB5D", "#EDC95E", "#EEC85E", "#EFC65E", 
            "#F0C45E", "#F1C25E", "#F2C05E", "#F3BE5E", "#F3BC5E", 
            "#F4BA5E", "#F4B85E", "#F5B65E", "#F5B45E", "#F6B25E", 
            "#F6B05D", "#F7AD5D", "#F7AB5D", "#F8A85D", "#F8A65D", 
            "#F8A35C", "#F9A15C", "#F99E5C", "#F99C5C", "#F9995B", 
            "#F9965B")
  int.color(orig, n, alpha)
}

cubeyf1col <- function(n = 255, alpha = 1){
  orig <- c("#7B0290", "#7B0392", "#7C0495", "#7D0597", "#7E0699", "#7F079C", "#80089E", "#8108A0", "#810AA2", "#820BA5", "#830CA8", "#830CAB", "#840CAE", "#850DB1", "#850FB3", "#8611B5", "#8613B8", "#8615BA", "#8618BC", "#861BBE", "#861DC0", "#8620C2", "#8622C4", "#8624C6", "#8527C8", "#8529CA", "#852BCD", "#852DCE", "#842FD0", "#8431D2", "#8433D4", "#8335D6", "#8337D8", "#8239DA", "#823ADB", "#813CDD", "#813EDF", "#8040E0", "#8042E2", "#7F44E4", "#7F45E6", "#7E47E8", "#7D48EA", "#7C4AEB", "#7A4DEE", "#794FEF", "#7850F1", "#7751F2", "#7653F4", "#7554F5", "#7556F7", "#7457F8", "#7359FA", "#725AFB", "#715BFB", "#715CFC", "#705EFC", "#6F5FFC", "#6E61FD", "#6E62FD", "#6D64FD", "#6C65FD", "#6C67FE", "#6B68FE", "#6A6AFE", "#696BFF", "#696CFF", "#686EFF", "#6870FE", "#6871FE", "#6773FD", "#6774FD", "#6676FC", "#6577FB", "#6478FA", "#637AF9", "#627CF8", "#617EF7", "#6081F5", "#6082F4", "#5F83F3", "#5E84F2", "#5D86F1", "#5C87F0", "#5B88EF", "#5B89EE", "#5A8AED", "#588BEC", "#588DEA", "#578EE9", "#568FE8", "#5590E7", "#5491E6", "#5493E4", "#5294E3", "#5296E1", "#5197E0", "#5098DF", "#4F99DE", "#4E9ADD", "#4E9CDC", "#4D9DDA", "#4B9ED9", "#4B9FD8", "#4AA0D7", "#49A1D6", "#48A2D4", "#47A4D3", "#46A5D2", "#45A7D0", "#43A9CD", "#42AACB", "#42ABC9", "#41ACC8", "#3FADC7", "#3FAEC5", "#3EAFC4", "#3DB0C3", "#3CB1C1", "#3CB2C0", "#3BB3BE", "#3AB4BC", "#39B5BB", "#38B6B9", "#38B7B8", "#37B8B6", "#37B8B5", "#38B9B3", "#39BAB1", "#39BBB0", "#3ABCAF", "#3ABCAD", "#3BBDAB", "#3BBEAA", "#3CBFA8", "#3DC0A6", "#3DC0A4", "#3EC1A3", "#3EC1A1", "#3FC29F", "#3FC39E", "#40C49C", "#41C599", "#42C696", "#43C794", "#43C893", "#44C991", "#45CA8F", "#45CA8E", "#46CB8C", "#46CB8A", "#46CC88", "#46CD86", "#47CD85", "#48CE83", "#48CF81", "#49D080", "#49D07E", "#4AD17C", "#4AD27A", "#4BD278", "#4CD376", "#4CD374", "#4DD472", "#4ED470", "#4ED56F", "#4FD66D", "#4FD66B", "#50D76A", "#50D768", "#51D866", "#52D865", "#52D964", "#53DA62", "#54DB5F", "#55DC5C", "#56DD5A", "#57DD58", "#58DE57", "#58DE55", "#59DF54", "#59DF52", "#59E051", "#5AE150", "#5AE14F", "#5BE24E", "#5BE24C", "#5CE34A", "#5DE349", "#5DE447", "#5FE447", "#60E547", "#62E548", "#63E648", "#65E648", "#67E748", "#68E749", "#6AE74A", "#6CE74A", "#6DE84B", "#6FE84B", "#71E84C", "#73E94D", "#75E94D", "#76E94D", "#78E94D", "#7BE94E", "#7EEA4F", "#82EA50", "#84EB50", "#86EB50", "#88EB51", "#8AEB51", "#8CEB52", "#8EEB52", "#91EB52", "#93EC53", "#95EC53", "#97EC54", "#9AEC54", "#9CEC54", "#9FEC55", "#A1EC55", "#A3EC55", "#A4EC55", "#A4EC55", "#A6EC55", "#A9EC56", "#ABEC56", "#ACEC56", "#ADEC56", "#AEEC56", "#B0EC57", "#B3EC57", "#B4EC57", "#B4EC57", "#B6EC57", "#B7EC58", "#B8EC58", "#B8EC58", "#BBEC58", "#BDEC58", "#BFEC59", "#C0EC59", "#C1EC59", "#C2EC59", "#C5EC59", "#C7EC5A", "#C8EC5A", "#C9EC5A", "#CAEC5A", "#CCEC5A", "#CCEC5A", "#CDEC5A", "#CFEC5B", "#D1EB5B")
  int.color(orig, n, alpha)
}

blackbodycol <- function(n = 255, alpha = 1){
  orig <- c("#000000", "#230000", "#340000", "#3C0000", "#3F0100", "#400200", "#440500", "#450600", "#480800", "#4A0A00", "#4D0C00", "#4E0E00", "#511000", "#531100", "#551300", "#561400", "#591600", "#5B1800", "#5C1900", "#5E1A00", "#5F1C00", "#621E00", "#641F00", "#662100", "#672200", "#692300", "#6A2400", "#6C2600", "#6D2700", "#6F2800", "#702A00", "#722B00", "#732C00", "#752D00", "#772F00", "#772F00", "#783000", "#7A3100", "#7B3300", "#7D3400", "#7D3400", "#7E3500", "#803600", "#813800", "#813800", "#833900", "#843A00", "#863B00", "#863B00", "#883D00", "#893E00", "#893E00", "#8B3F00", "#8B3F00", "#8C4100", "#8E4200", "#8E4200", "#8F4300", "#8F4300", "#914400", "#914400", "#924600", "#924600", "#944700", "#944700", "#954800", "#954800", "#974900", "#974900", "#994B00", "#994B00", "#9A4C00", "#9A4C00", "#9A4C00", "#9C4D00", "#9C4D00", "#9D4F00", "#9D4F00", "#9F5000", "#9F5000", "#9F5000", "#A05100", "#A05100", "#A25200", "#A25200", "#A35400", "#A35400", "#A55500", "#A55500", "#A65600", "#A65600", "#A65600", "#A85700", "#A85700", "#AA5900", "#AA5900", "#AB5A00", "#AB5A00", "#AD5B00", "#AD5B00", "#AE5D00", "#AE5D00", "#B05E00", "#B05E00", "#B15F00", "#B15F00", "#B36000", "#B36000", "#B46200", "#B66300", "#B66300", "#B76400", "#B76400", "#B96600", "#B96600", "#BB6700", "#BB6700", "#BC6800", "#BC6800", "#BE6900", "#BF6B00", "#BF6B00", "#C16C00", "#C16C00", "#C26D00", "#C46E00", "#C46E00", "#C57000", "#C57000", "#C77100", "#C87200", "#C87200", "#CA7400", "#CA7400", "#CC7500", "#CD7600", "#CD7600", "#CF7700", "#D07900", "#D07900", "#D27A00", "#D37B00", "#D37B00", "#D57C00", "#D67E00", "#D67E00", "#D87F00", "#D98000", "#D98000", "#DB8200", "#DD8300", "#DD8300", "#DE8400", "#E08500", "#E08500", "#E18700", "#E38800", "#E38800", "#E48900", "#E68A00", "#E68A00", "#E78C00", "#E98D00", "#E98D00", "#EA8E00", "#EC9000", "#EC9000", "#EE9100", "#EF9200", "#F19300", "#F19300", "#F29500", "#F49600", "#F49600", "#F59700", "#F79900", "#F79900", "#F89A00", "#FA9B00", "#FB9C00", "#FB9C00", "#FD9E00", "#FF9F00", "#FF9F00", "#FFA000", "#FFA100", "#FFA300", "#FFA300", "#FFA400", "#FFA500", "#FFA700", "#FFA700", "#FFA800", "#FFA900", "#FFA900", "#FFAA00", "#FFAC00", "#FFAD00", "#FFAD00", "#FFAE00", "#FFAF00", "#FFB100", "#FFB200", "#FFB300", "#FFB500", "#FFB500", "#FFB600", "#FFB700", "#FFB800", "#FFBB07", "#FFBC0A", "#FFBD0E", "#FFBF12", "#FFC015", "#FFC119", "#FFC31D", "#FFC524", "#FFC628", "#FFC82B", "#FFCA33", "#FFCC36", "#FFCE3D", "#FFCF41", "#FFD248", "#FFD34C", "#FFD653", "#FFD85B", "#FFDB62", "#FFDD69", "#FFDF6D", "#FFE174", "#FFE47B", "#FFE886", "#FFEA8E", "#FFED95", "#FFEF9C", "#FFF0A0", "#FFF3A7", "#FFF6AE", "#FFF8B6", "#FFF9B9", "#FFFCC1", "#FFFDC4", "#FFFFCC", "#FFFFCF", "#FFFFD3", "#FFFFDA", "#FFFFDE", "#FFFFE1", "#FFFFE5", "#FFFFE9", "#FFFFEC", "#FFFFF0", "#FFFFF4", "#FFFFF7", "#FFFFFF")
  int.color(orig, n, alpha)
}
  
HMrainbow <- function(n = 255, alpha = 1){
  orig <- c("#000000", "#2D0024", "#38002E", "#3C0031", "#430036", "#46003B", "#47003D", "#4B0044", "#4A0049", "#4A004D", "#490051", "#470057", "#45015A", "#44025E", "#420361", "#3F0666", "#3D076A", "#3A0A6D", "#380C71", "#350F74", "#301277", "#2F1479", "#2C177C", "#291B80", "#281C81", "#252084", "#222486", "#1D2B89", "#19348A", "#18398B", "#183E8D", "#18408E", "#17418E", "#17458F", "#17478E", "#17478E", "#17498E", "#174B8E", "#174B8E", "#174E8E", "#17508E", "#17508E", "#17528D", "#17558D", "#17558D", "#17578C", "#17578C", "#185A8C", "#185A8C", "#185D8B", "#185D8B", "#185D8B", "#185D8B", "#18618B", "#18618B", "#19658A", "#19658A", "#196889", "#196889", "#196889", "#1A6C89", "#1A6C89", "#1B6F88", "#1B6F88", "#1B6F88", "#1B7387", "#1B7387", "#1C7686", "#1C7686", "#1D7A85", "#1D7A85", "#1D7A85", "#1D7A85", "#1D7D84", "#1D7D84", "#1E8083", "#1E8083", "#1F8382", "#1F8382", "#1F8382", "#208680", "#208680", "#21897F", "#21897F", "#21897F", "#228C7D", "#228C7D", "#238E7B", "#238E7B", "#249179", "#249179", "#249179", "#259376", "#259376", "#269674", "#269674", "#289871", "#289871", "#299A6F", "#299A6F", "#2A9C6C", "#2A9C6C", "#2B9E6A", "#2B9E6A", "#2B9E6A", "#2DA068", "#2DA068", "#2EA265", "#2EA265", "#30A463", "#30A463", "#32A661", "#32A661", "#33A85F", "#35AA5D", "#35AA5D", "#35AA5D", "#37AC5B", "#37AC5B", "#39AE58", "#39AE58", "#3BAF56", "#3EB154", "#40B252", "#40B252", "#43B450", "#43B450", "#45B54F", "#48B74D", "#48B74D", "#48B74D", "#4BB84C", "#4DBA4A", "#50BB49", "#53BD48", "#57BE48", "#5BBF47", "#5FC046", "#63C146", "#67C246", "#6BC346", "#6FC446", "#6FC446", "#73C446", "#77C546", "#7BC546", "#82C647", "#85C747", "#89C748", "#8CC748", "#8FC749", "#8FC749", "#93C749", "#96C74A", "#99C74A", "#9CC74B", "#A0C84C", "#A7C84E", "#AAC84F", "#ADC84F", "#ADC84F", "#B1C850", "#B4C851", "#B7C752", "#BAC752", "#BEC753", "#C4C755", "#C7C655", "#C7C655", "#CBC656", "#CEC557", "#D4C559", "#D7C45A", "#DAC35B", "#E0C25E", "#E0C25E", "#E6C160", "#E9C062", "#ECBE64", "#EEBD68", "#F0BC6A", "#F0BC6A", "#F2BB6E", "#F4B972", "#F5B874", "#F7B778", "#F8B67B", "#F8B67B", "#FAB57D", "#FBB480", "#FCB482", "#FDB485", "#FDB485", "#FEB486", "#FEB38A", "#FFB38E", "#FFB391", "#FFB391", "#FFB398", "#FFB4A1", "#FFB4A4", "#FFB4A7", "#FFB4A7", "#FFB5A9", "#FFB5AA", "#FFB6AD", "#FFB7B0", "#FFB7B0", "#FFB8B3", "#FFB9B3", "#FFB9B6", "#FFBAB6", "#FFBAB6", "#FFBBB9", "#FFBCB9", "#FFBDBC", "#FFBDBC", "#FFBEBC", "#FFBFBF", "#FFC0BF", "#FFC2C2", "#FFC2C2", "#FFC5C5", "#FFC6C6", "#FFC8C8", "#FFC9C9", "#FFC9C9", "#FFCACA", "#FFCBCB", "#FFCDCD", "#FFCECE", "#FFCECE", "#FFD0D0", "#FFD1D1", "#FFD3D3", "#FFD7D7", "#FFD8D8", "#FFD8D8", "#FFDADA", "#FFDBDB", "#FFDDDD", "#FFDFDF", "#FFE2E2", "#FFE4E4", "#FFE6E6", "#FFE6E6", "#FFE8E8", "#FFEBEB", "#FFEDED", "#FFF0F0", "#FFF3F3", "#FFF6F6", "#FFF9F9", "#FFFBFB", "#FFFDFD", "#FFFFFF")
  int.color(orig, n, alpha)
}  
  
  
optim_scalecol <- function(n = 255, alpha = 1){
  orig <- c("#000000", "#040000", "#080000", "#0C0000", "#100000", "#140000", "#180000", "#1C0000", "#200000", "#240000", "#280000", "#2C0000", "#300000", "#340000", "#380000", "#3C0000", "#400000", "#440000", "#480000", "#4C0000", "#500000", "#540000", "#580000", "#5C0000", "#600000", "#640000", "#680000", "#6C0000", "#700000", "#740000", "#780300", "#780500", "#780800", "#780A00", "#780D00", "#780F00", "#781200", "#781400", "#781700", "#781900", "#781C00", "#781E00", "#782100", "#782300", "#782600", "#782800", "#782B00", "#782D00", "#783000", "#783200", "#783500", "#783700", "#783A00", "#783C00", "#783F00", "#784100", "#784400", "#784600", "#784900", "#784B00", "#784E00", "#785000", "#785300", "#785500", "#785800", "#785A00", "#785D00", "#785F00", "#786200", "#786400", "#786700", "#786900", "#786C00", "#786E00", "#787100", "#787300", "#787600", "#787800", "#787B00", "#787D03", "#788005", "#788208", "#78850A", "#78870D", "#788A0F", "#788C12", "#788F14", "#789117", "#799419", "#7A961C", "#7B991E", "#7C9B21", "#7D9E23", "#7EA026", "#7FA328", "#80A52B", "#81A82D", "#82AA30", "#83AD32", "#84AF35", "#84B237", "#84B43A", "#84B73C", "#84B93F", "#84BC41", "#84BE44", "#84C146", "#84C349", "#84C64B", "#84C84E", "#84CB50", "#84CD53", "#84D055", "#84D258", "#84D55A", "#84D75D", "#84DA5F", "#84DC62", "#84DF64", "#84E167", "#84E469", "#84E66C", "#84E96E", "#84EB71", "#84EE73", "#84F076", "#84F378", "#84F57B", "#84F87D", "#84FA80", "#84FD82", "#84FF85", "#84FF87", "#84FF8A", "#84FF8C", "#84FF8F", "#84FF91", "#84FF94", "#84FF96", "#84FF99", "#84FF9B", "#84FF9E", "#84FFA0", "#84FFA3", "#84FFA5", "#84FFA8", "#84FFAA", "#84FFAD", "#84FFAF", "#84FFB2", "#84FFB4", "#84FFB7", "#84FFB9", "#84FFBC", "#84FFBE", "#84FFC1", "#84FFC3", "#84FFC6", "#84FFC8", "#84FFCB", "#84FFCD", "#84FFD0", "#84FFD2", "#84FFD5", "#84FFD7", "#84FFDA", "#84FFDC", "#84FFDF", "#84FFE1", "#84FFE4", "#84FFE6", "#84FFE9", "#84FFEB", "#84FFEE", "#84FFF0", "#84FFF3", "#84FFF5", "#84FFF8", "#84FFFA", "#84FFFD", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#84FFFF", "#87FFFF", "#89FFFF", "#8CFFFF", "#8FFFFF", "#91FFFF", "#94FFFF", "#96FFFF", "#99FFFF", "#9BFFFF", "#9EFFFF", "#A0FFFF", "#A3FFFF", "#A5FFFF", "#A8FFFF", "#AAFFFF", "#ADFFFF", "#AFFFFF", "#B2FFFF", "#B4FFFF", "#B7FFFF", "#B9FFFF", "#BCFFFF", "#BEFFFF", "#C1FFFF", "#C3FFFF", "#C6FFFF", "#C8FFFF", "#CBFFFF", "#CDFFFF", "#D0FFFF", "#D2FFFF", "#D5FFFF", "#D7FFFF", "#DAFFFF", "#DCFFFF", "#DFFFFF", "#E1FFFF", "#E4FFFF", "#E6FFFF", "#E9FFFF", "#EBFFFF", "#EEFFFF", "#F0FFFF", "#F3FFFF", "#F5FFFF", "#F8FFFF", "#FAFFFF", "#FDFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")
  int.color(orig, n, alpha)
}

Lin_optimcol <- function(n = 255, alpha = 1){
  orig <- c("#000000", "#000000", "#000000", "#010000", "#020000", "#020000", "#030000", "#030000", "#040000", "#050000", "#050000", "#060000", "#070000", "#070000", "#080000", "#090000", "#090000", "#0A0000", "#0B0000", "#0C0000", "#0D0000", "#0E0000", "#0F0000", "#100000", "#110000", "#120000", "#130000", "#140000", "#150000", "#160000", "#170000", "#190000", "#1A0000", "#1B0000", "#1C0000", "#1E0000", "#1F0000", "#210000", "#220000", "#230000", "#250000", "#270000", "#280000", "#2B0000", "#2D0000", "#2E0000", "#310000", "#330000", "#350000", "#360000", "#380000", "#3A0000", "#3C0000", "#3E0000", "#400000", "#430000", "#450000", "#470000", "#4A0000", "#4C0000", "#500000", "#510000", "#540000", "#560000", "#590000", "#5C0000", "#5E0000", "#610000", "#640000", "#670000", "#6A0000", "#6D0000", "#700000", "#730000", "#750000", "#7A0000", "#7E0000", "#800000", "#830000", "#870000", "#870000", "#870100", "#870200", "#870300", "#870400", "#870600", "#870600", "#870800", "#870900", "#870A00", "#870B00", "#870D00", "#870D00", "#870F00", "#871100", "#871100", "#871300", "#871500", "#871600", "#871700", "#871900", "#871A00", "#871B00", "#871D00", "#871F00", "#872000", "#872100", "#872300", "#872400", "#872600", "#872800", "#872A00", "#872C00", "#872E00", "#872F00", "#873100", "#873300", "#873400", "#873600", "#873800", "#873900", "#873B00", "#873E00", "#873F00", "#874100", "#874300", "#874500", "#874800", "#874900", "#874C00", "#874E00", "#875000", "#875200", "#875400", "#875700", "#875800", "#875A00", "#875D00", "#875F00", "#876200", "#876500", "#876700", "#876A00", "#876B00", "#876E00", "#877100", "#877300", "#877600", "#877900", "#877C00", "#877F00", "#878100", "#878500", "#878700", "#878A00", "#878D00", "#879000", "#879400", "#879600", "#879B00", "#879D00", "#87A000", "#87A300", "#87A600", "#87AA00", "#87AE00", "#87B100", "#87B400", "#87B800", "#87BC00", "#87C000", "#87C300", "#87C800", "#87CB00", "#87CD00", "#87D200", "#87D600", "#87DA00", "#87DE00", "#87E200", "#87E700", "#87EC00", "#87EF00", "#87F400", "#87F900", "#87FE00", "#87FF01", "#87FF05", "#87FF0A", "#87FF0F", "#87FF14", "#87FF17", "#87FF1C", "#87FF21", "#87FF26", "#87FF2B", "#87FF2D", "#87FF31", "#87FF36", "#87FF3B", "#87FF41", "#87FF46", "#87FF4A", "#87FF50", "#87FF54", "#87FF5A", "#87FF5F", "#87FF62", "#87FF68", "#87FF6E", "#87FF74", "#87FF78", "#87FF7D", "#87FF83", "#87FF89", "#87FF90", "#87FF95", "#87FF9A", "#87FF9E", "#87FFA5", "#87FFAC", "#87FFB3", "#87FFBA", "#87FFBF", "#87FFC6", "#87FFCB", "#87FFD3", "#87FFD8", "#87FFE0", "#87FFE8", "#87FFF0", "#87FFF8", "#87FFFE", "#87FFFF", "#8CFFFF", "#92FFFF", "#99FFFF", "#9CFFFF", "#A1FFFF", "#A8FFFF", "#ACFFFF", "#B1FFFF", "#B6FFFF", "#BDFFFF", "#C0FFFF", "#C7FFFF", "#CCFFFF", "#D2FFFF", "#D7FFFF", "#DCFFFF", "#E1FFFF", "#E8FFFF", "#ECFFFF", "#F0FFFF", "#F8FFFF", "#FFFFFF")
  int.color(orig, n, alpha)
}