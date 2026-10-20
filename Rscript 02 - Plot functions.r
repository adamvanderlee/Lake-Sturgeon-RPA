#        1         2         3         4         5         6         7         8
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# PLOT FUNCTIONS
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# PAIRS PLOT FUNCTION
#-------------------------------------------------------------------------------

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

## computes variance inflation factors (VIF) and reduced model
## to only params with VIF<cutoff
corVIF = function(x,y,cutoff=5,cor="FALSE"){
     require(car)
     if(cor){
         cat("Correlation",":","\n")
         print(cor(x,use="complete.obs"))
     }    
     dat <- x
     count = 0
     flag <- TRUE
     cat("\n","VIF :","\n")
     while(flag){
         count = count+1
         fit <- lm( y ~ ., data=dat)
         vfit <- vif(fit)
         cat("\n","Interation",count,":","\n")
         print(vfit)
         if(max(vfit) > cutoff){
             dat <- dat[, -which.max(vfit)]
         } else {
             flag <- FALSE
         }
     }
}

#-------------------------------------------------------------------------------
# PAUSE PLOT
#-------------------------------------------------------------------------------
# stops between successive plots
keybd=function(key) {
      "Done"
}
pausePlot=function(wait=0) {
   if (wait==0) {
     getGraphicsEvent("Any Key To Continue",onMouseDown = NULL,
                       onMouseMove = NULL,onKeybd = keybd)
   }
   else {
     start=Sys.time()
     stop=Sys.time()+wait
     now=Sys.time()
     while (now < stop) {
       now=Sys.time()
     }
   }
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Multiple plot function
#
# FROM : http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Deterines default ggplot colours for n factors
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

ggplot_colour <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#-------------------------------------------------------------------------------
