#' Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @title Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0    Initial value.
#' @param r    Growth rate parameter.
#' @param k    Carrying capacity.
#' @param N    Length of the time series.
#' @param type    One of: "driving" (default), "damping", "logistic", "vanGeert1991".
#'
#' @return A timeseries object of length N.
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#' @seealso \code{\link{growth.ac.cond}}
#'
#' @examples
#' # The logistic map in the chaotic regime
#' growth.ac(Y0 = 0.01, r = 4, type = "logistic")
growth.ac <- function(Y0 = 0.01, r = 1, k = 1, N = 100, type = c("driving", "damping", "logistic", "vanGeert")[1]){
    # Create a vector Y of length N, which has value Y0 at Y[1]
    if(N>1){
        Y <- as.numeric(c(Y0, rep(NA,N-2)))
        # Conditional on the value of type ...
        switch(type,
               # Iterate N steps of the difference function with values passed for Y0, k and r.
               driving  = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] ),
               damping  = k + sapply(seq_along(Y), function(t) Y[[t+1]] <<- - r * Y[t]^2 / k),
               logistic = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] * ((k - Y[t]) / k)),
               vanGeert = sapply(seq_along(Y), function(t) Y[[t+1]] <<- Y[t] * (1 + r - r * Y[t] / k))
        )}
    return(ts(Y))
}


#' @title Un-initialise It: Unload and/or uninstall R packages
#' @description \code{un.IT} will check if the Packages in the list argument \code{loose} are installed on the system and unload them. If \code{unT=TRUE} it will first unload the packages if they are loaded, and then proceed to uninstall them.
#' @param loose    A vector of package names to be unloaded.
#' @param unT    Logical. If \code{TRUE} (default), packages in \code{loose} wil be un-installed if they are available on the system.
#'
#' @export
#'
#'@author Fred Hasselman
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}
#'
#' @examples
#' \dontnrun{un.IT(loose = c("reshape2", "plyr", "dplyr"), unT = FALSE)}
un.IT <- function(loose,unT=FALSE){
    dp <- .packages()
    if(any(loose %in% dp)){
        for(looseLib in loose[(loose %in% dp)]){detach(paste0("package:",looseLib), unload=TRUE,character.only=TRUE)}
    }
    rm(dp)
    if(unT==TRUE){
        dp <- .packages(all.available=TRUE)
        if(any(loose %in% dp)){remove.packages(loose[(loose %in% dp)])}
    }
}
#' Conditional Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0 Initial value
#' @param r Growth rate parameter
#' @param k Carrying capacity
#' @param cond Conditional rules passed as a data.frame of the form: cbind.data.frame(Y = ..., par = ..., val = ...)
#' @param N Length of the time series
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#' @seealso \code{\link{growth.ac.cond}}
#'
#' @examples
#' # Plot with the default settings
#' xyplot(growth.ac.cond())
#'
#' # The function such that it can take a set of conditional rules and apply them sequentially during the iterations.
#' # The conditional rules are passed as a `data.frame`
#' (cond <- cbind.data.frame(Y = c(0.2, 0.6), par = c("r", "r"), val = c(0.5, 0.1)))
#' xyplot(growth.ac.cond(cond=cond))
#'
#' # Combine a change of `r` and a change of `k`
#' (cond <- cbind.data.frame(Y = c(0.2, 1.99), par = c("r", "k"), val = c(0.5, 3)))
#' xyplot(growth.ac.cond(cond=cond))
#'
#' # A fantasy growth process
#' (cond <- cbind.data.frame(Y = c(0.1, 1.99, 1.999, 2.5, 2.9), par = c("r", "k", "r", "r","k"), val = c(0.3, 3, 0.9, 0.1, 1.3)))
#' xyplot(growth.ac.cond(cond=cond))
growth.ac.cond <- function(Y0 = 0.01, r = 0.1, k = 2, cond = cbind.data.frame(Y = 0.2, par = "r", val = 2), N = 100){
    # Create a vector Y of length N, which has value Y0 at Y[1]
    Y <- c(Y0, rep(NA, N-1))
    # Iterate N steps of the difference equation with values passed for Y0, k and r.
    cnt <- 1
    for(t in seq_along(Y)){
        # Check if the current value of Y is greater than the threshold for the current conditional rule in cond
        if(Y[t] > cond$Y[cnt]){
            # If the threshold is surpassed, change the parameter settings by evaluating: cond$par = cond$val
            eval(parse(text = paste(cond$par[cnt], "=", cond$val[cnt])))
            # Update the counter if there is another conditional rule in cond
            if(cnt < nrow(cond)){cnt <- cnt + 1}
        }
        # Van Geert growth model
        Y[[t+1]] <- Y[t] * (1 + r - r * Y[t] / k)
    }
    return(ts(Y))
}


# PLOTS -------------------------------------------------------------------
disp <- function(message='Hello world!', header = TRUE, footer = TRUE){

    mWidth <- max(laply(message,nchar))

    if(is.character(header)){
        hWidth <- max(laply(header,nchar))
        mWidth <- max(hWidth,mWidth)
    }

    dmessage <- list()
    for(m in 1:length(message)){
       # b <- floor((mWidth-nchar(message[m]))/2)
        e <- mWidth-nchar(message[m])
        dmessage[[m]] <- paste0('§ ',message[m]) #,paste0(rep(' ',e),collapse=""),'\n\t')
                                #paste0('§ ',paste0(rep(" ",mWidth),collapse=""),' §'))
    }
    # if(m > 1){dmessage[[m]] <- paste0(dmessage[[m]],}

   # mWidth <- max(laply(dmessage, nchar))
    banner <- paste0(rep('~', mWidth), collapse = "")
    if(is.character(header)){
        b <- floor((nchar(banner)-nchar(header))/2)
        e <- ceiling((nchar(banner)-nchar(header))/2)
            leader <- paste0('\n\t',paste0(rep('~',b),collapse=""),header,paste0(rep('~',e),collapse=""))
        }
    if(header == TRUE){
            leader <- banner
        }
    if(header == FALSE){
            leader <- paste0('§') #,paste0(rep(" ",nchar(banner)-2),collapse="")) #,'§')
        }

    if(footer){
            cat(paste0('\n\t',leader,'\n\t',dmessage,'\n\t',banner,'\n'))
        } else {
            cat(paste0('\n\t',leader,'\n\t',dmessage))
        }
}

gg.theme <- function(type=c("clean","noax")[1],useArial = F, afmPATH="~/Dropbox"){
    require(ggplot2)
    if(useArial){
        set.Arial(afmPATH)
        bf_font="Arial"
    } else {bf_font="Helvetica"}

    switch(type,
           clean = theme_bw(base_size = 16, base_family=bf_font) +
               theme(axis.text.x     = element_text(size = 14),
                     axis.title.y    = element_text(vjust = +1.5),
                     panel.grid.major  = element_blank(),
                     panel.grid.minor  = element_blank(),
                     legend.background = element_blank(),
                     legend.key = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.line  = element_line(colour = "black")),

           noax = theme(line = element_blank(),
                        text  = element_blank(),
                        title = element_blank(),
                        plot.background = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank())
    )
}

plotHolder <- function(useArial = F,afmPATH="~/Dropbox"){
    require(ggplot2)
    ggplot() +
        geom_blank(aes(1,1)) +
        theme(line = element_blank(),
              text  = element_blank(),
              title = element_blank(),
              plot.background = element_blank(),
              #           panel.grid.major = element_blank(),
              #           panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()
              #           axis.title.x = element_blank(),
              #           axis.title.y = element_blank(),
              #           axis.text.x = element_blank(),
              #           axis.text.y = element_blank(),
              #           axis.ticks = element_blank()
        )
}

set.Arial <- function(afmPATH="~/Dropbox"){
    # Set up PDF device on MAC OSX to use Arial as a font in Graphs
    if(nchar(afmPATH>0)){
        if(file.exists(paste0(afmPATH,"/Arial.afm"))){
            Arial <- Type1Font("Arial",
                               c(paste(afmPATH,"/Arial.afm",sep=""),
                                 paste(afmPATH,"/Arial Bold.afm",sep=""),
                                 paste(afmPATH,"/Arial Italic.afm",sep=""),
                                 paste(afmPATH,"/Arial Bold Italic.afm",sep="")))
            if(!"Arial" %in% names(pdfFonts())){pdfFonts(Arial=Arial)}
            if(!"Arial" %in% names(postscriptFonts())){postscriptFonts(Arial=Arial)}
            return()
        } else {disp(header='useArial=TRUE',message='The directory did not contain the *.afm version of the Arial font family')}
    } else {disp(header='useArial=TRUE',message='Please provide the path to the *.afm version of the Arial font family')}
}


# GRAPH PLOTTING ---------------------------------------------------------------
graph2svg <- function(TDM,pname){

  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  g <- simplify(g)

  # Remove vertices that were used in the search query
  Vrem <- which(V(g)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  g <- (g - V(g)$name[Vrem])

  # Set colors and sizes for vertices
  V(g)$degree <- degree(g)
  rev         <- scaleRange(log1p(V(g)$degree))
  rev[rev<=0.3]<-0.3

  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
  V(g)$size        <- 10*scaleRange(V(g)$degree)
  V(g)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, rev)
  V(g)$label.cex   <- scaleRange(V(g)$degree)+.1

  # set edge width and color
  rew <- scaleRange(E(g)$weight)
  rew[rew<=0.3]<-0.3

  E(g)$width <- 2*scaleRange(E(g)$weight)
  E(g)$color <- rgb(.5, .5, 0, rew)
  set.seed(958)

  svg(paste(pname,sep=""),width=8,height=8)
  plot(g, layout=layout.fruchterman.reingold(g))
  dev.off()

  return(g)
}

# Plot vertex neighbourhood
hoodGraph2svg <- function(TDM,Vname,pname){

   # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  ig <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  ig <- simplify(ig)

  # Remove vertices that were used in the search query
  Vrem <- which(V(ig)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  ig <- (ig - V(ig)$name[Vrem])

  # This is a deletion specific for the Neighbourhood graphs
  Vrem <- which(V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
  ig   <- ig - V(ig)$name[Vrem]

  idx <- which(V(ig)$name==Vname)
  sg  <- graph.neighborhood(ig, order = 1, nodes=V(ig)[idx], mode = 'all')[[1]]

  # set colors and sizes for vertices
  V(sg)$degree <- degree(sg)

  rev<-scaleRange(log1p(V(sg)$degree))
  rev[rev<=0.3]<-0.3

  V(sg)$color <- rgb(scaleRange(V(sg)$degree), 1-scaleRange(log1p(V(sg)$degree*V(sg)$degree)),  0, rev)

  V(sg)$size        <- 35*scaleRange(V(sg)$degree)
  V(sg)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(sg)$label       <- V(sg)$name
  V(sg)$label.color <- rgb(0, 0, 0, rev)
  V(sg)$label.cex   <- scaleRange(V(sg)$degree)

  # set edge width and color
  rew<-scaleRange(E(sg)$weight)
  rew[rew<=0.3]<-0.3

  E(sg)$width <- 6*scaleRange(E(sg)$weight)
  E(sg)$color <- rgb(.5, .5, 0, rew)

  idV <- which(V(sg)$name==Vname)
  idE <- incident(sg,V(sg)[[idV]])
  E(sg)$color[idE] <- rgb(0, 0, 1 ,0.8)

  set.seed(958)

  idx <- which(V(sg)$name==Vname)
  svg(paste(pname,sep=""),width=8,height=8)
  plot(sg,layout=layout.star(sg,center=V(sg)[idx]))
  dev.off()

  return(sg)
}

# Conversion formulas for self-affinity parameter estimates (sap) to Dimension (fd) suggested in Hasselman (2013)
# PSD slope (if signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure)
psd2fd <- function(sap){return(fd <- 3/2 + ((14/33)*tanh(sap*log(1+sqrt(2)))) )}
# DFA slope (this is H in DFA)
dfa2fd <- function(sap){return(fd <- 2-(tanh(log(3)*sap)) )}
# SDA slope (simple 1-sap, but note that for some signals different PSD slope values project to 1 SDA slope)
sda2fd <- function(sap){return(fd <- 1-sap)}


sliceTS<-function(TSmat,epochSz=1) {
  # Slice columns of TSmat in epochs of size = epochSz
  init(c("plyr"))

  N<-dim(TSmat)
  return(llply(seq(1,N[1],epochSz),function(i) TSmat[i:min(i+epochSz-1,N[1]),1:N[2]]))
}

fltrIT <- function(TS,f){
  # Apply filtfilt to TS using f (filter settings)
  require("signal")

  return(filtfilt(f=f,x=TS))

}

SWtest0 <- function(g){
  Nreps <- 10;
  histr  <- vector("integer",Nreps)
  target<- round(mean(degree(g)))
  now   <- target/2
  for(i in 1:Nreps){
    gt      <- watts.strogatz.game(dim=1, size=length(degree(g)), nei=now, 0)
    histr[i] <- round(mean(degree(gt)))
    ifelse(histr[i] %in% histr,break,{
      ifelse(histr[i]>target,{now<-now-1},{
        ifelse(histr[i]<target,{now<-now+1},{
          break})
      })
    })
  }
  return(gt)
}


# SWtestV <- function(g,N){
#  return(list(cp=transitivity(g,type="global"),cpR=transitivity(rewire(g,mode=c("simple"),niter=N),type="global"),lp=average.path.length(g), lpR=average.path.length(rewire(g,mode=c("simple"),niter=N))))
# }

SWtestE <- function(g,p=1,N=20){
  values <- matrix(nrow=N,ncol=6,dimnames=list(c(1:N),c("cp","cpR","cp0","lp","lpR","lp0")))

  for(n in 1:N) {
    gt<-SWtest0(g)
    values[n,] <- c(transitivity(g,type="localaverage"),transitivity(rewire.edges(g,prob=p),type="localaverage"),transitivity(gt,type="localaverage"),average.path.length(g),average.path.length(rewire.edges(g,prob=p)),average.path.length(gt))}
  values[n,values[n,]==0] <- NA #values[n,values[n,]==0]+1e-8}

  values   <- cbind(values,(values[,1]/values[,2])/(values[,4]/values[,5]),(values[,1]/values[,3]),(values[,4]/values[,6]),((values[,1]/values[,3])/values[,2])/((values[,4]/values[,6])/values[,5]))
  valuesSD <- data.frame(matrix(apply(values[,1:10],2,sd,na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  valuesAV <- data.frame(matrix(colMeans(values[,1:10],na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  return(list(valuesAV=valuesAV,valuesSD=valuesSD,valuesSE=valuesSD/sqrt(N)))
}

PLFsmall <- function(g){
  reload <- FALSE
  if("signal" %in% .packages()){
    warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
    reload <- TRUE
    detach("package:signal", unload=TRUE)}

  if(length(V(g))>100){warning("Vertices > 100, no need to use PLFsmall, use a binning procedure");break}
  d<-degree(g,mode="all")
  y<-hist(d,breaks=0:length(V(g)),plot=F)$density
  y<-y[y>0]
  if(length(y)==2){warning("Caution... Log-Log slope is a bridge (2 points)")}
  if(length(y)<2){warning("Less than 2 points in Log-Log regression... aborting");break}
  alpha=coef(lm(log(y) ~ poly(log(1:length(y)), degree=1), na.action="na.exclude") )[2]

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(alpha)
}

FDrel <- function(g){
  d<-degree(g,mode="all")
  nbreaks <- round(length(V(g))/2)-1
  y<-hist(d,breaks=nbreaks,plot=F)$density
  y<-y[y>0]
  return(FD <- -sum(y*log2(y))/-(log2(1/length(y))))
}

# scale.R -------------------------------------------------------------

# # Three uses:
# #
# # 1. scale.RANGE(x)             Scale x to data range: min(x.out)==0;      max(x.out)==1
# # 2. scale.RANGE(x,mn,mx)       Scale x to arg. range: min(x.out)==mn==0;  max(x.out)==mx==1
# # 3. scale.RANGE(x,mn,mx,lo,hi) Scale x to arg. range: min(x.out)==mn==lo; max(x.out)==mx==hi
# #
# # Examples:
# #
# # Works on numeric objects
# somenumbers <- cbind(c(-5,100,sqrt(2)),c(exp(1),0,-pi))
#
# scale.RANGE(somenumbers)
# scale.RANGE(somenumbers,mn=-100)
#
# # Values < mn will return < lo (default=0)
# # Values > mx will return > hi (default=1)
# scale.RANGE(somenumbers,mn=-1,mx=99)
#
# scale.RANGE(somenumbers,lo=-1,hi=1)
# scale.RANGE(somenumbers,mn=-10,mx=101,lo=-1,hi=4)

scale.R <- function(x,mn=min(x,na.rm=T),mx=max(x,na.rm=T),lo=0,hi=1){
    x <- as.data.frame(x)
    u <- x
    for(i in 1:dim(x)[2]){
        mn=min(x[,i],na.rm=T)
        mx=max(x[,i],na.rm=T)
        if(mn>=mx){warning("Minimum (mn) >= maximum (mx).")}
        if(lo>=hi){warning("Lowest scale value (lo) >= highest scale value (hi).")}
        ifelse(mn==mx,{u[,i]<-rep(mx,length(x[,i]))},{
            u[,i]<-(((x[i]-mn)*(hi-lo))/(mx-mn))+lo
            id<-complete.cases(u[,i])
            u[!id,i]<-0
        })
    }
    return(u)
}

Rmd2htmlWP <- function(infile, outfile, sup = T) {
    require(markdown)
    require(knitr)
    mdOpt <- markdownHTMLOptions(default = T)
    mdOpt <- mdOpt[mdOpt != "mathjax"]
    mdExt <- markdownExtensions()
    mdExt <- mdExt[mdExt != "latex_math"]
    if (sup == T) {
        mdExt <- mdExt[mdExt != "superscript"]
    }
    knit2html(input = infile, output = outfile, options = c(mdOpt), extensions = c(mdExt))
}

# MULTIPLOT FUNCTION ------------------------------------------------------------------------------------------------------------------
#
# [copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ]
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multi.PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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

# TRY … CATCH -------------------------------------------------------------------------------------------------------------------------

##================================================================##
###  In longer simulations, aka computer experiments,            ###
###  you may want to                                             ###
###  1) catch all errors and warnings (and continue)             ###
###  2) store the error or warning messages                      ###
###                                                              ###
###  Here's a solution  (see R-help mailing list, Dec 9, 2010):  ###
##================================================================##

# Catch *and* save both errors and warnings, and in the case of
# a warning, also keep the computed result.
#
# @title tryCatch both warnings (with value) and errors
# @param expr an \R expression to evaluate
# @return a list with 'value' and 'warning', where
#   'value' may be an error caught.
# @author Martin Maechler;
# Copyright (C) 2010-2012  The R Core Team
#
try.CATCH <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
}



# OLD NETWORK FUNCTIONS -------------------------------------------------------------------------------------------

#
# brainButter <- function(TSmat, fs=500, band=c(lfHIp=4,hfLOp=40), Np=9){
#   # Extract frequency bands from columns of TSmat that are commonly used in Neuroscience
#   # Low Freq. High-pass (1st) and High Freq. Low-pass (2nd) FIR1 filter is applied at frequencies specified as band=c(lfHIp=...,hfLOp=...)
#   require("signal")
#
#   fHI <- butter(Np,band[[1]]*2/fs,"high")
#   fLO <- butter(Np,band[[2]]*2/fs,"low")
#
#   TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fHI,x=TS))
#   TSflt <- apply(TSflt,2,function(TS) filtfilt(f=fLO,x=TS))
#
#   #TSflt <- apply(TSflt,2,fltrIT,f=fLO)
#
#   return(TSflt)
# }
#
# brainFir1 <- function(TSmat, fs=500, band=c(lfHIp=4,hfLOp=40), Np=2/band[1]){
#   # Extract frequency bands from columns of TSmat that are commonly used in Neuroscience
#   # Low Freq. High-pass (1st) and High Freq. Low-pass (2nd) FIR1 filter is applied at frequencies specified as band=c(lfHIp=...,hfLOp=...)
#   require("signal")
#
#   if(2/band[1]>Np){print(paste("Incorrect filter order Np... using 2/",band[1]," = ",(2/band[1]),sep=""))}
#
#   fBP <- fir1(floor(Np*fs),band/(fs/2),type="pass");
#
#   TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fBP,1,x=TS))
#
#
#   #   fHI <- fir1(floor(Np*fs),band[1]/(fs/2),type="high");
#   #   fLO <- fir1(floor(Np*fs),band[2]/(fs/2),type="low");
#   #
#   #   TSflt <- apply(TSmat,2,function(TS) filtfilt(f=fHI,1,x=TS))
#   #   TSflt <- apply(TSflt,2,function(TS) filtfilt(f=fLO,1,x=TS))
#
#
#   return(TSflt)
# }
#
#
# ssi2sbi <- function(SImat,threshold){
#   # Signed Similarity matrix to "signed binary" matrix
#
#   idS   <- which(SImat<0)
#   BImat <- abs(as.matrix(SImat))
#   diag(BImat) <- 0
#   BImat[BImat <= threshold] <- 0
#   BImat[BImat >  threshold] <- 1
#   BImat[idS] <- BImat[idS]*-1
#
#   return(BImat)
# }
#
# si2bi <- function(SImat,threshold){
#   # Unsigned Similarity matrix to unsigned binary matrix
#
#   ifelse(any(SImat<0),{
#     print("Signed matrix, use: ssi2sbi()")
#     break},{
#       BImat <- as.matrix(SImat)
#       diag(BImat) <- 0
#       BImat[BImat <= threshold] <- 0
#       BImat[BImat >  threshold] <- 1})
#
#   return(BImat)
# }
#
# ssi2sth <- function(SImat,threshold){
#   # Signed Similarity matrix to "signed thresholded" matrix
#
#   idS   <- which(SImat<0)
#   THmat <- abs(as.matrix(SImat))
#   diag(THmat) <- 0
#   THmat[THmat <= threshold] <- 0
#   THmat[idS] <- THmat[idS]*-1
#
#   return(THmat)
# }
#
# si2th <- function(SImat,threshold){
#   # Similarity matrix to thresholded matrix
#
#   ifelse(any(SImat<0),{
#     print("Signed matrix, use: ssi2sth()")
#     break},{
#       THmat <- as.matrix(SImat)
#       THmat[THmat <= threshold] <- 0})
#
#   return(THmat)
# }
#
#
# plotBIN <- function(BImat){
#
#   g <- graph.adjacency(BImat, weighted=T, mode = "undirected",diag=F)
#   g <- simplify(g)
#
#   # set colors and sizes for vertices
#   V(g)$degree <- degree(g)
#
#   rev<-scaleRange(log1p(V(g)$degree))
#   rev[rev<=0.3]<-0.3
#
#   V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
#   V(g)$size        <- 25*rev
#   V(g)$frame.color <- NA
#
#   # set vertex labels and their colors and sizes
#   V(g)$label       <- V(g)$name
#   V(g)$label.color <- rgb(0, 0, 0, rev)
#   V(g)$label.cex   <- rev
#
#   # set edge width and color
#
#   E(g)$width <- 4
#   E(g)$color <- rgb(.5, .5, 0, .6)
#   set.seed(958)
#
#   #   layout1=layout.spring(g)
#   #    layout2=layout.fruchterman.reingold(g)
#   #    layout3=layout.kamada.kawai(g)
#   #   layout5 = layout.spring(g,mass=0.3,repulse=T)
#
#   #   CairoFontMatch(fontpattern="Arial")
#   #   CairoFonts(regular="Arial:style=Normal")
#
#   #   CairoPDF(pname,10,10)
#   #   plot(g, layout=layout.sphere)
#   #   dev.off()
#   #
#
#   plot(g, layout=layout.sphere)
#
#   return(g)
# }
#
# plotMAT <- function(BImat,l=NULL){
#
#   g <- graph.adjacency(BImat, weighted=T, mode = "undirected",diag=F)
#   #g <- simplify(g)
#
#   # set colors and sizes for vertices
#   V(g)$degree <- degree(g)
#
#   rev<-scaleRange(V(g)$degree)
#   rev[rev<=0.4]<-0.4
#
#   V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
#   V(g)$size        <- 20*rev
#   V(g)$frame.color <- NA
#
#   # set vertex labels and their colors and sizes
#   V(g)$label       <- V(g)$name
#   V(g)$label.color <- rgb(0, 0, 0, .8)
#   V(g)$label.cex   <- 1.1
#
#   # set edge width and color
#   #  rew<-E(g)$weight
#   #  rew[rew<=0.3]<-0.3
#   #
#   edge.central=edge.betweenness(g)
#   #
#   for (i in 1:ecount(g)) {E(g)$width[i]=0.3+sqrt((edge.central[i]))}
#
#   # E(g)$width <- 2*E(g)$weight
#   E(g)$color <- rgb(.5, .5, 0, .6)
#   set.seed(958)
#
#   if(is.null(l)){l<-layout.fruchterman.reingold(g,niter=500,area=vcount(g)^2.3,repulserad=vcount(g)^2.8)}
#
#   plot(g,layout=l)
#   return(g)
# }
#
#
# plotSIGNth <- function(sSImat){
#
#   g <- graph.adjacency(sSImat, weighted=TRUE)
#   E(g)$sign <- E(g)$weight
#   E(g)$curved <- is.mutual(g)
#   E(g)$lty <- ifelse( E(g)$sign > 0, 1, 1)
#   E(g)$arrow.size <- .2
#   E(g)$width <- 3
#   #E(g)$color <- rgb(scaleRange(abs(E(g)$weight)), 1-scaleRange(abs(E(g)$weight)), 0, 1)
#   #layout1=layout.fruchterman.reingold(g)
#
#   V(g)$label.color <- rgb(0, 0, 0, 1)
#   V(g)$label.cex <- 1.4
#   V(g)$vs   <- graph.strength(g, mode="in")
#   V(g)$vs.u <- scaleRange(graph.strength(g))
#   #V(g)$color<- ifelse( V(g)$vs > 0, rgb(V(g)$vs.u, 1-V(g)$vs.u, 0, 1), rgb(1-V(g)$vs.u, V(g)$vs.u, 0, 1))
#
#   E(g)$es.u  <- scaleRange(E(g)$weight)
#   E(g)$color <- ifelse( E(g)$sign > 0, rgb(0, 1, 0, .2), rgb(1, 0, 0, .2))
#   return(g)
# }
#
# plotSW <- function(n,k,p){
#
#   g <- watts.strogatz.game(1, n, k, p)
#
#   V(g)$degree <- degree(g)
#
#   # set colors and sizes for vertices
#   rev<-scaleRange(log1p(V(g)$degree))
#   rev[rev<=0.2]<-0.2
#   rev[rev>=0.9]<-0.9
#   V(g)$rev <- rev
#
#   V(g)$color       <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
#   V(g)$size        <- 25*V(g)$rev
#
#   # set vertex labels and their colors and sizes
#   V(g)$label       <- ""
#
#   E(g)$width <- 1
#   E(g)$color <- rgb(0.5, 0.5, 0.5, 1)
#
#   return(g)
# }
#
# plotBA <- function(n,pwr,out.dist){
#   #require("Cairo")
#
#   g <- barabasi.game(n,pwr,out.dist=out.dist,directed=F)
#   V(g)$degree <- degree(g)
#
#   # set colors and sizes for vertices
#   rev<-scaleRange(log1p(V(g)$degree))
#   rev[rev<=0.2] <- 0.2
#   rev[rev>=0.9] <- 0.9
#   V(g)$rev <- rev
#
#   V(g)$color    <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
#   V(g)$size     <- 25*V(g)$rev
#   # V(g)$frame.color <- rgb(.5, .5,  0, .4)
#
#   # set vertex labels and their colors and sizes
#   V(g)$label <- ""
#
#   E(g)$width <- 1
#   E(g)$color <- rgb(0.5, 0.5, 0.5, 1)
#
#   return(g)
# }
