\name{dbkGrad}
\alias{dbkGrad}
\alias{print.dbkGrad}
\alias{as.data.frame.dbkGrad}
\title{
Discrete Beta Kernel Graduation of Mortality Rates
}
\description{
This function performs nonparametric graduation of mortality rates using discrete beta kernel smoothing techniques.
}
\usage{
dbkGrad(obsqx, omega, ex = NULL, logit = F, bandwidth = "FX", h = 0.002, s = 0.2, cvres = "propres", cvh = T, cvs = F, alpha = 0.05)
\method{print}{dbkGrad}(x,...)
\method{as.data.frame}{dbkGrad}(x, row.names = 0:x$omega, optional = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsqx}{
a numeric vector containing the observed mortality rates.
The first element is referred to age zero.
}
  \item{omega}{
an optional scalar; it sets the upper age limit and its default value is \code{length(obsqx)-1}.
}
  \item{ex}{
an optional numeric vector containing the exposed to the risk of death for each age.
}
  \item{logit}{
an optional logical; if \code{logit=T}, then the logit transformation is applied to the data before graduating, and then data are back-transformed to obtain the estimate of the true rates; its default value is \code{FALSE}. 
}
  \item{h}{
an optional scalar. 
It provides a value for the smoothing parameter. 
Default value is 0.002. 
If \code{cvh=T}, then the smoothing parameter is computed by means of cross-validation using this parameter as the initial value.
}
  \item{s}{
an optional scalar. 
It provides a value for the dampening factor. 
Default value is 0.2. 
If \code{cvs=T}, then the dampening factor is computed by means of cross-validation using this parameter as the initial value.
}
  \item{cvh}{
an optional logical; if \code{cvh=T} then cross-validation is used to select the smoothing parameter. 
Default value is \code{TRUE}. 
Parameter \code{h} is the initial value used in cross-validation.
}
  \item{cvs}{
an optional logical; if \code{TRUE} then cross-validation is used to select the dampening factor. 
Default value is \code{FALSE}. 
Parameter \code{s} is the initial value used in cross-validation.
}
  \item{cvres}{
an optional character string; if \code{cvres="propres"} (the default), then cross-validation selects the smoothing parameter and/or the dampening factor by minimizing the squares of the proportional differences between observed and estimated values, while if \code{cvres="res"} then the sum of square residuals is minimized.
}
 \item{bandwidth}{
an optional character string. 
It specifies the bandwidth type and must be:
\itemize{
  \item \code{"FX"} for a fixed bandwidth,
  \item \code{"EX"} for an adaptive bandwidth based on the number of the exposed to the risk of death \code{ex} (see Mazza A, Punzo A, 2013, for details);
  \item \code{"VC"} for an adaptive bandwidth based on a vector of weights derived from the variation coefficients which, in turn, depends from the exposed to the risk of death \code{ex} (see Mazza A, Punzo A, in press, for details). 
}
Default value is \code{"FX"}.
}
  \item{alpha}{an optional scalar. 
When the \code{ex} argument is provided, the function returns \code{(1-alpha)*100\%} pointwise confidence intervals for fitted values. 
Default value is 0.05.
}	 
  \item{x}{a \code{dbkGrad} object}
  \item{row.names}{a NULL or a character vector giving the row names for the data frame. Missing values are not allowed. Default value is \code{0:x$omega}.}
  \item{optional}{logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional.}
  \item{...}{additional arguments to be passed to or from methods.}
}
\details{
In the cross-validation routine, minimization is performed using the Levenberg-Marquardt algorithm (More 1978) in the minpack.lm package (Elzhov, Mullen, and Bolker 2010). 
}
\value{
Returned from this function is an \code{dbkGrad} object which is a list with the following components:
  \item{obsqx}{a numeric vector containing the observed mortality rates until age \code{omega}.
The first element is referred to age zero.}
  \item{fitted.values}{a numeric vector containing the graduated values.} 
  \item{ex}{a numeric vector containing the exposed to the risk of death for each age until \code{omega}.}
  \item{h}{a scalar. 
  It is the smoothing parameter used in discrete beta kernel smoothing.}
  \item{s}{a scalar. 
  It is the dampening parameter used in discrete beta kernel smoothing. It is returned if \code{bandwidth = "EX"} or \code{bandwidth = "VC"}}
  \item{residuals}{a numeric vector containing the residuals \code{fitted.values - obsqx}.}
  \item{prop.residuals}{a numeric vector containing the proportional residuals  \code{fitted.values / obsqx - 1}.}
  \item{kernels}{a matrix. 
  \code{kernels \%*\% obsqx} returns the \code{fitted.values}}
  \item{cvRSS}{a scalar. 
  It is the cross-validation residual sum of squares (RSS)computed over the fitted values, using the residuals specified in \code{cvres}.}
  \item{call}{an object of class \code{call}.}
  \item{lowerbound}{a numeric vector containing the lower bound of the \code{(1-alpha)*100\%} pointwise confidence intervals for fitted values. It is returned only if \code{ex} is provided bye the user.}
  \item{upperbound}{a numeric vector containing the upper bound of the \code{(1-alpha)*100\%} pointwise confidence intervals for fitted values. It is returned only if \code{ex} is provided bye the user.}
}

\references{
Elzhov TV, Mullen KM, Bolker B (2010) \emph{minpack.lm:} R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK. R package version 1.1-5. URL http://CRAN.R-project.org/package=minpack.lm.

Bagnato L, Punzo A (in press) Finite Mixtures of Unimodal Beta and Gamma Densities and the k-Bumps Algorithm. \emph{Computational Statistics}.

Mazza A, Punzo A (2011) Discrete Beta Kernel Graduation of Age-Specific Demographic Indicators. In S Ingrassia, R Rocci, M Vichi (eds.), \emph{New Perspectives in Statistical Modeling and Data Analysis}, volume \bold{42} of \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 127-134. Springer-Verlag, Berlin-Heidelberg.

Mazza A, Punzo A (2013) Graduation by Adaptive Discrete Beta Kernels. In A Giusti, G Ritter, M Vichi (eds.), \emph{Classification and Data Mining}, volume \bold{44} of \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 77-84. Springer-Verlag, Berlin-Heidelberg.

Mazza A, Punzo A (in press) Using the Variation Coefficient for Adaptive Discrete Beta Kernel Graduation. In P Giudici, S Ingrassia, M Vichi (eds.), \emph{Advances in Statistical Modelling for Data Analysis, Studies in Classification, Data Analysis and Knowledge Organization}, Springer-Verlag, Berlin-Heidelberg.

More J (1978) The Levenberg-Marquardt Algorithm: Implementation and Theory. In G Watson (ed.), \emph{Numerical Analysis}, volume \bold{630} of \emph{Lecture Notes in Mathematics}, pp. 104-116. Springer- Verlag, Berlin-Heidelberg.

Punzo A (2010) Discrete Beta-type Models. In H Locarek-Junge, C Weihs (eds.), \emph{Classification as a Tool for Research}, volume \bold{40} of \emph{Studies in Classification, Data Analysis and Knowledge Organization}, pp. 253-261. Springer-Verlag, Berlin-Heidelberg.

Punzo A, Zini A (2012) Discrete Approximations of Continuous and Mixed Measures on a Compact Interval \emph{Statistical Papers}, \bold{53}(3), 563-575.
}
\author{
Angelo Mazza and Antonio Punzo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DBKGrad-package}}, \code{\link{plot}}, \code{\link{Sicily2008M}}
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

