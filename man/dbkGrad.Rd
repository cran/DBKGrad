\name{dbkGrad}
\alias{dbkGrad}
\alias{print.dbkGrad}
\alias{as.data.frame.dbkGrad}
\alias{residuals.dbkGrad}
\title{
Discrete Beta Kernel Graduation of Mortality Rates
}
\description{
This function performs nonparametric graduation of mortality rates using discrete beta kernel smoothing techniques.
}
\usage{
dbkGrad(obsq, limx, limy, exposures = NULL, transformation = c("none", "log", "logit",
       "Gompertz"), bwtypex = c("FX", "VC", "EX"), bwtypey = c("FX", "VC", "EX"), 
       adaptx = c("a", "b", "ab"), adapty = c("a", "b", "ab"), hx = 0.002, hy = 0.002,
       sx = 0.2, sy = 0.2, cvres = c("propres", "res"), cvhx = FALSE, cvhy = FALSE, 
       cvsx = FALSE, cvsy = FALSE, alpha = 0.05)

\method{print}{dbkGrad}(x, \dots)

\method{as.data.frame}{dbkGrad}(x, row.names = x$limx[1]:x$limx[2], optional = FALSE, \dots)

\method{residuals}{dbkGrad}(object, type = c("working", "proportional", "response", 
      "deviance", "pearson"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsq}{
a matrix (or an object which can be coerced to a matrix using \code{as.matrix()}) of observed mortality rates. Dimnames, if provided, should be numeric; row names should be ages and column names years. 
}
  \item{limx, limy}{
optional vector of two integers; if provided, \code{limx} (\code{limy}) sets a lower and a upper row (column) limit. Only data within these intervals are graduated.
}
  \item{exposures}{
an optional matrix containing the exposed to the risk of death for each age and year. Dimensions of \code{exposures} should correspond to those of \code{obsq}.
}
  \item{transformation}{
  an optional character string; the transformation specified is applied to the observed data before graduation. Graduated data are then back-transformed. Possible values are          
  \code{"none"} for no transformation, \code{"log"}, \code{"logit"} and \code{"Gompertz"}.
  }
  \item{bwtypex, bwtypey}{
an optional character string. 
It specifies the type of bandwidth to be adopted by row (by column) and must be:
\itemize{
  \item \code{"FX"} for a fixed bandwidth (default),
  \item \code{"EX"} for an adaptive bandwidth based on \code{exposures} (see Mazza A, Punzo A, 2013a, for details);
  \item \code{"VC"} for an adaptive bandwidth based on a vector of weights derived from the variation coefficients which, in turn, depends from the \code{exposures} (see Mazza A, Punzo A, in press, for details). 
}
}
 \item{adaptx, adapty}{
an optional character string. It is the  type of adaptive bandwidth to be adopted by row (by column) and must be:
    \itemize{
      \item \code{"a"} a different bandwidth for each evaluation age x at which the rates are estimated;
      \item \code{"b"} a different bandwidth can be attributed to each age, regardless from the evaluation point;
      \item \code{"ab"} a different bandwidth can be selected for each evaluation point and for each age.
    }
  }

  \item{hx, hy}{
an optional scalar. It is the global bandwidth used for the variable on the rows (columns). 
Default value is 0.002. 
If \code{cvhx=TRUE (cvhy=TRUE)}, then the smoothing parameter is computed by means of cross-validation using this value as an initialization.
}
  \item{sx, sy}{
an optional scalar. 
It is the sensitive parameter used for the variable on the rows (columns). 
Default value is 0.2. 
If \code{cvsx=TRUE} (\code{cvsy=TRUE}), then the sensitive parameter is computed by means of cross-validation using this value as an initialization.
}
  \item{cvhx, cvhy}{
an optional logical; if \code{cvsx=TRUE (cvsy=TRUE)} then cross-validation is used to select the smoothing parameter. 
Default value is \code{TRUE}. 
Parameter \code{hx, (hy)} is the initial value used in cross-validation.
}
  \item{cvsx, cvsy}{
an optional logical; if \code{TRUE} then cross-validation is used to select the sensitive parameter. 
Default value is \code{FALSE}. 
The value of sx (sy) is used to initialize the cross-validation process.
}
  \item{cvres}{
an optional character string; if \code{cvres="propres"} (the default), then cross-validation selects the smoothing parameter and/or the sensitive parameter  by minimizing the squares of the proportional differences between observed and estimated values, while if \code{cvres="res"} then the sum of square residuals is minimized.
}

  \item{alpha}{an optional scalar. 
When the \code{exposures} argument is provided, the function returns \code{(1-alpha)*100\%} pointwise confidence intervals and pointwise confidence bands for fitted values. 
Default value is 0.05.
}	 
  \item{x}{a \code{dbkGrad} object}
  \item{row.names}{a NULL or a character vector giving the row names for the data frame. Missing values are not allowed. Default value is \code{0:x$limx, limy}.}
  \item{optional}{logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional.}
  \item{...}{additional arguments to be passed to or from methods.}

  \item{type}{ \code{"working"},\code{"proportional"},\code{"response"},\code{"deviance"}, \code{"pearson"}}
   \item{object}{a \code{dbkGrad} class object.}
}
\details{
In the cross-validation routine, minimization is performed using the Levenberg-Marquardt algorithm (More 1978) in the minpack.lm package (Elzhov, Mullen, and Bolker 2010). 
}
\value{
Returned from this function is an \code{dbkGrad} object which is a list with the following components:
  \item{fitted.values}{a matrix containing the graduated values.} 
  \item{residuals}{a matrix containing the working residuals \code{fitted.values - obsq}.}
  \item{kernels}{a matrix. 
  \code{kernels \%*\% obsq} returns the \code{fitted.values}}
  \item{cvRSS}{a scalar. 
  It is the cross-validation residual sum of squares (RSS) computed over the fitted values, using the residuals specified in \code{cvres}.}
  \item{hx (hy)}{a scalar. 
  It is the global bandwidth used for the variable on the rows (columns).}
  \item{sx (sy)}{a scalar. 
  It is the sensitive parameter used for the variable on the rows (columns). It is returned when \code{bandwidth = "EX"} or \code{bandwidth = "VC"}}
  \item{upperbound,lowerbound}{pointwise confidence interval. Returned when \code{exposures} is provided.}
  \item{bonferroniupperbound, bonferronilowerbound}{limits of the Bonferroni confidence bands. Returned when \code{exposures} is provided.}
  \item{sidakupperbound, sidaklowerbound}{limits of the Sidak confidence bands. Returned when \code{exposures} is provided.}
  \item{obsq}{a matrix containing the observed mortality rates with dimensions set by \code{limx, limy}.}
  \item{exposures}{a matrix containing the exposures with dimensions set by \code{limx, limy}.}
  \item{limx (limy)}{a vector with lower and upper row (column) limits. Only data within these interval are graduated.}
  \item{call}{an object of class \code{call}.}
}

\references{
Elzhov TV, Mullen KM, Bolker B (2010) \emph{minpack.lm:} R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK. R package version 1.1-5. URL http://CRAN.R-project.org/package=minpack.lm.

Mazza A, Punzo A (2011) Discrete Beta Kernel Graduation of Age-Specific Demographic Indicators. In S Ingrassia, R Rocci, M Vichi (eds.), \emph{New Perspectives in Statistical Modeling and Data Analysis}, \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 127-134. Springer-Verlag, Berlin-Heidelberg.

Mazza A, Punzo A (2013a) Graduation by Adaptive Discrete Beta Kernels. In A Giusti, G Ritter, M Vichi (eds.), \emph{Classification and Data Mining}, \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 77-84. Springer-Verlag, Berlin-Heidelberg.

Mazza A, Punzo A (2013b) Using the Variation Coefficient for Adaptive Discrete Beta Kernel Graduation. In P Giudici, S Ingrassia, M Vichi (eds.), \emph{Advances in Statistical Modelling for Data Analysis, Studies in Classification, Data Analysis and Knowledge Organization}, pp. 225-232, Springer International Publishing, Switzerland.

Mazza A, Punzo A (2014) DBKGrad: An R Package for Mortality Rates Graduation by Discrete Beta Kernel Techniques. \emph{Journal of Statistical Software, Code Snippets}, \bold{57}{2}, 1-18.

More J (1978) The Levenberg-Marquardt Algorithm: Implementation and Theory. In G Watson (ed.), \emph{Numerical Analysis}, volume \bold{630} of \emph{Lecture Notes in Mathematics}, pp. 104-116. Springer- Verlag, Berlin-Heidelberg.

Punzo A (2010) Discrete Beta-type Models. In H Locarek-Junge, C Weihs (eds.), \emph{Classification as a Tool for Research}, \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 253-261. Springer-Verlag, Berlin-Heidelberg.

}
\author{
Angelo Mazza and Antonio Punzo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DBKGrad-package}}, \code{\link[DBKGrad:plot.dbkGrad]{plot}}, \code{\link{ItalyM}}
}
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

