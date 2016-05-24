`jpng` <-
function(file="tmp", width = 8, height = 8, P=NULL,  bg = "white")
{
  if(missing(file)) { file="Jrast" }
        if(missing(P)) { P = NULL }
        if(missing(width)) { width = NULL }
        if(missing(height)) { height = NULL }

        if(!is.null(width) & !is.null(height) )
            {
                P = c(width, height)
            }
        
        if(is.null(P) & is.null(width) & is.null(height)   )
            {
                P = par('din')
                P = round(P, digits=2)
                width= P[1]
                height =P[2]
            }
        if(!is.null(P))
            {
                width= P[1]
                height =P[2]
            }
        if(length( grep('.png', file)) > 0)
            {
                psname = file
            }
        else
            {
                psname = RPMG::local.file(file, "png")
            }
        ## P = round(par('pin'))
        png(filename = psname,
            width =width , height =height , units = "in", res=300,
            pointsize = 12, bg = bg )
        return(psname)
  
}

jpdf<-function(file="tmp", width = 8, height = 8, P=NULL)
{
  if(missing(file)) { file="Jpdf" }
        if(missing(P)) { P = NULL }
        if(missing(width)) { width = NULL }
        if(missing(height)) { height = NULL }

        if(!is.null(width) & !is.null(height) )
            {
                P = c(width, height)
            }
        
        if(is.null(P) & is.null(width) & is.null(height)   )
            {
                P = par('din')
                P = round(P, digits=2)
                width= P[1]
                height =P[2]
            }
        if(!is.null(P))
            {
                width= P[1]
                height =P[2]
            }
        if(length( grep('.pdf', file)) > 0)
            {
                psname = file
            }
        else
            {
                psname = RPMG::local.file(file, "pdf")
            }
        ## P = round(par('pin'))
        pdf(file = psname,
            width =width , height =height , paper = "special",
          onefile = TRUE )
        return(psname)
  
}
