#' Initialise It: Load and/or install R packages
#' @title Initialise It: Load and/or install R packages
#' @description \code{in.IT} will check if the Packages in the list argument \code{need} are installed on the system and load them. If \code{inT=TRUE} it will first install the packages if they are not present and then proceed to load them.
#'
#' @param need    A vector of package names to be loaded.
#' @param inT    Logical. If \code{TRUE} (default), packages in \code{need} wil be installed if they are not available on the system.
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family initialise packages
#' @seealso \code{\link{un.IT}}
#'
#' @examples
#' in.IT(c("reshape2", "plyr", "dplyr"))
in.IT <- function(need=NULL,inT=TRUE){
    ip <- .packages(all.available=TRUE)
    if(any((need %in% ip)==FALSE)){
        if(inT==TRUE){
            install.packages(need[!(need %in% ip)])
        } else {
            cat('Package(s):\n',paste(need[(need %in% ip)==FALSE],sep='\n'),'\nnot installed.\nUse in.IT(c("packagename1","packagename2",...),inT=TRUE)')
            need <- need[(need %in% ip)==TRUE]
        }
    }
    ok <- sapply(1:length(need),function(p) require(need[[p]],character.only=TRUE))
}

#' Wrapper for in.IT: Load I/O and data handling tools
#' @rdname in.IT
#'
#' @title Wrapper for in.IT: Load I/O and data handling tools
#'
#' @param need    Passed to in.IT: c("foreign","xlsx","plyr","doBy","reshape2","RCurl","XML","httr","dplyr")
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}, \code{\link{un.IT}}
#' @example
#' in.IO()
in.IO <- function(need = c("xlsx","plyr","doBy","reshape2","RCurl","XML","httr","dplyr")){
    in.IT(need)
}
#' Wrapper for in.IT: Parallel computing tools
#'
#' @title Wrapper for in.IT: Parallel computing tools
#'
#' @param need    Passed to in.IT: c("parallel","doParallel","foreach")
#'
#' @export
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}, \code{\link{un.IT}}
#' @example
#' in.PAR()
in.PAR <- function(need = c("parallel","doParallel","foreach")){
    in.IT(need)
}

#' Wrapper for in.IT: Tools for plotting
#'
#' @rdname in.IT
#'
#' @title Wrapper for in.IT: Tools for plotting
#'
#' @param need    Passed to in.IT: c("lattice","latticeExtra","gplots","ggplot2","grid","gridExtra","scales","beanplot","effects","RColorBrewer")
#'
#' @export
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}, \code{\link{un.IT}}
#' @example
#' in.GR()
in.GR <- function(need = c("lattice","latticeExtra","gplots","ggplot2","grid","gridExtra","scales","beanplot","effects","RColorBrewer")){
    in.IT(need)
}
#' Wrapper for in.IT: Nonlinear Time Series packages
#'
#' @rdname in.IT
#'
#' @param need    Passed to in.IT: c("fractaldim","fractalrock","RTisean","tsDyn","tseries","tseriesChaos")
#'
#' @export
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}, \code{\link{in.IT}}
#'
#' @examples
#' in.NLTS()
in.NLTS <- function(need = c("fractaldim","fractalrock","RTisean","tsDyn","tseries","tseriesChaos")){
    in.IT(need)
}

#' Wrapper for in.IT: Signal analysis packages
#'
#' @rdname in.IT
#'
#' @param need
#'
#' @export
#'
#' @family initialise packages
#' @seealso \code{\link{in.IT}}, \code{\link{un.IT}}
#'
#' @examples
#' in.SN()
in.SIGN <- function(need=c("pracma","signal","EMD","hht","matlab")){
    in.IT(need)
}
