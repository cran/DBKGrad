\name{plot.dbkGrad}
\alias{plot.dbkGrad}
\title{Plot Method for dbkGrad objects}
\description{
Plotting the \code{\link{dbkGrad}} object produces a few different plots that are of interest. The different plots, created from various \code{plottype} options, are described below:
}
\usage{
  \method{plot}{dbkGrad}(x, plottype = "obsfit", CI = T,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{a \code{\link{dbkGrad}} object}
	\item{plottype}{an optional character string. It specifies the type of plot to dispaly and must be:
    \itemize{
      \item \code{"obsfit"} to plot observed and fitted values; 
      \item \code{"fitted"} to plot fitted values; 
      \item \code{"observed"} to plot observed values; 
      \item \code{"histres"} to plot the histogram of the classical residuals;
      \item \code{"histpropres"} to plot the histogram of the proportional residuals;   
      \item \code{"exposed"} to plot the exposed to the risk of death. 
    }  
    Default value is \code{"obsfit"}.
  } 
	\item{CI}{an optional logical; if \code{TRUE}, \code{1-alpha}*100\% confidence intervals for the fitted values are displayed when \code{plottype="obsfit"} and \code{plottype="fitted"}. When the \code{alpha} argument is not provided in \code{dbkGrad()}, 95\% pointwise confidence intervals are displayed. Default value is \code{TRUE}.}
  \item{...}{additional arguments to be passed to or from methods.}
}

\value{
No values are returned from the plot function.

}
\author{
Angelo Mazza and Antonio Punzo
}

\seealso{
\code{\link{DBKGrad-package}}, \code{\link{dbkGrad}}, \code{\link{Sicily2008M}}
}
\examples{
# data loading

data("Sicily2008M")
obsqx <- Sicily2008M$qx
ex    <- Sicily2008M$ex
head(Sicily2008M)
tail(Sicily2008M)

# fixed discrete beta kernel graduation

resFX1 <- dbkGrad(obsqx=obsqx, omega=85)
plot(resFX1, plottype="observed")
plot(resFX1, plottype="fitted")
plot(resFX1, plottype="obsfit")
plot(resFX1, plottype="histres")

resFX2 <- dbkGrad(obsqx, ex=ex, omega=85, alpha=0.05)
plot(resFX2, plottype="obsfit", CI=TRUE)

resFX3 <- dbkGrad(obsqx, ex, omega=85, h=0.01, cvh=FALSE, alpha=0.05)
plot(resFX3, plottype="obsfit", CI=TRUE)
plot(resFX3, plottype="exposed")

# adaptive discrete beta kernel graduation

resEX <- dbkGrad(obsqx, omega=85, ex, bandwidth="EX", s=0.28, cvres="res", cvh=TRUE, cvs=FALSE, alpha=0.05)
plot(resEX, plottype="obsfit", CI=TRUE)

resVC <- dbkGrad(obsqx, omega=85, ex, logit=TRUE, bandwidth="VC", cvh=TRUE, cvs=TRUE, alpha=0.05)
plot(resVC, plottype="obsfit", CI=TRUE)

# summary of results

as.data.frame(resVC)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
