\name{plot.dbkGrad}
\alias{plot.dbkGrad}
\title{Plot Method for dbkGrad objects}
\description{
Plotting the \code{\link{dbkGrad}} object produces a few different plots that are of interest. The different plots, created from various \code{plottype} options, are described below:
}
\usage{
\method{plot}{dbkGrad}(x, plottype = c("obsfit", "fitted", "observed", "exposure",
    "residuals", "checksd"), plotstyle = c("mat", "level", "persp"),
    restype = c("working", "proportional", "response", "deviance", "pearson"), 
    byage = TRUE, columns, rows, CI = TRUE, CBBonf = FALSE, CBSidak = FALSE,
    logscale = TRUE, alphares = 0.05, col, \dots)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{a \code{\link{dbkGrad}} object}
	\item{plottype}{an optional character string. It specifies the type of plot to display and must be one of:
    \itemize{
      \item \code{"observed"} to plot observed values; 
      \item \code{"fitted"} to plot fitted values; 
      \item \code{"obsfit"} to plot observed and fitted values (default);
      \item \code{"exposure"} to plot the exposed to the risk of death;
      \item \code{"residuals"} to display plots related to residuals: density of residuals, residuals versus fitted values, and residuals versus the discrete variable of interest; 
      \item \code{"checksd"} to plot autocorrelogram and autodependogram (see Bagnato, Punzo, Nicolis, 2012) of residuals, only for the unidimensional case.
    }  
  }
	\item{plotstyle}{an optional character string. It specifies the style of plot; it has no effect when \code{plottype=checksd}. It must be:
    \itemize{
      \item \code{"mat"} for a matplot (default for the unidimensional case); 
      \item \code{"level"} for a levelplot (default for the bidimensional case); 
      \item \code{"persp"} for a perspective plot.
    }
  }
    \item{restype}{an optional character string.When \code{plottype=residuals} or \code{plottype=checksd}, it specifies the type of residuals displayed. It must be:
    \itemize{
      \item \code{"working"} {to use working residuals: \code{o-f} (default);} 
      \item \code{"proportional"} {to use working residuals:  \code{o/f-1};} 
      \item \code{"response"} {to use response residuals:  \code{e*(o-f)};}
      \item \code{"deviance"} {to use deviance residuals:  \code{sign(o-f) * sqrt(2*e*o*log(o/f) + 2*e*(1-o)*log((1-o)/(1-f)))};}
      \item \code{"pearson"} {to use Pearson residuals:  \code{e*(o-f)/sqrt(e*f*(1-f))},}
    }
    where \code{e}=exposures, \code{o}=observed values and \code{f}=fitted values.
  }
  \item{rows, columns}{an optional numeric vector. It specifies the rows (columns) to plot. 
  }
  \item{byage}{an optional logical. It selects the discrete variable of interest:  \code{TRUE} (default) for the variable in rows (typically age) and \code{FALSE} for the other variable (years or duration). 
  }
  \item{logscale}{an optional logical;when TRUE (default), rates are plotted on log scale.}
  \item{CI}{an optional logical; if \code{TRUE}, \code{1-alpha*100\%} confidence intervals for the fitted values are displayed in \code{plottype="obsfit"} and \code{plottype="fitted"}. When the \code{alpha} argument is not provided in \code{dbkGrad()}, 95\% pointwise confidence intervals are displayed. Default value is \code{TRUE}.}
  \item{CBBonf}{an optional logical; does the same as \code{CI} but Bonferroni correction is applied to obtain confidence bands. Default is \code{FALSE}.} 
  \item{CBSidak}{an optional logical; does the same as \code{CI} but Sidak correction is applied to obtain confidence bands. Default is \code{FALSE}.} 
  \item{alphares}{an optional scalar. 
When \code{plottype=residuals} the boundaries of the \code{(alphares)*100\%} critical region are displayed. 
Default value is 0.05.}
\item{col}{a scalar or a vector with plotting colors.}
  \item{\dots}{additional arguments to be passed to or from methods.}
}

\value{
No values are returned from the plot function.

}

\references{
Bagnato L, Punzo A, Nicolis O (2012). The autodependogram: a graphical device to investigate serial dependences. \emph{Journal of Time Series Analysis}, \bold{33}(2), 233-254.

Bagnato L, De Capitani L, Punzo A (2013a). Detecting Serial Dependencies with the Reproducibility Probability Autodependogram. \emph{Advances in Statistical Analysis}. doi:10.1007/s10182-013-0208-y.

Bagnato L, De Capitani L, Punzo A (2013b). Testing Serial Independence via Density-Based Measures of Divergence. \emph{Methodology and Computing in Applied Probability}. doi:10.1007/s11009-013-9320-4.

Mazza A, Punzo A (2014) DBKGrad: An R Package for Mortality Rates Graduation by Discrete Beta Kernel Techniques. \emph{Journal of Statistical Software, Code Snippets}, \bold{57}{2}, 1-18.
}
\author{
Angelo Mazza and Antonio Punzo
}
\seealso{
  \code{\link{DBKGrad-package}}, \code{\link{dbkGrad}}, \code{\link{ItalyM}}, \code{\link[TSA:acf]{TSA:acf}}, \code{\link[SDD:ADF]{SDD:ADF}}}
\examples{

data("ItalyM")

# unidimensional analysis

res1 <- dbkGrad(obsq=obsq, limx=c(6,71), limy=104, exposure=population, bwtypex="EX", adaptx="ab")
plot(res1, plottype="obsfit", CI=FALSE, CBBonf=TRUE)
plot(res1, plottype="residuals", restype="pearson")
plot(res1, plottype="checksd", restype="pearson")
residuals(res1, type="pearson")

# bidimensional analysis

res2 <- dbkGrad(obsq=obsq, limx=c(6,46), limy=c(60,80), exposure=population, 
transformation="logit", bwtypex="VC", bwtypey="EX", hx=0.01, hy=0.008, adaptx="ab", adapty="b")
plot(res2, plottype="obsfit")
plot(res2, plottype="obsfit", plotstyle="persp", col="black")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
