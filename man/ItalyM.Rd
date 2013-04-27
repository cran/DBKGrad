\name{ItalyM}
\docType{data}
\alias{ItalyM}
\alias{population}
\alias{obsq}
\title{Mortality data for the 1906-2009 male population of Italy}

\description{
This data set consists of probabilities of death and population size for the male Italian population aged from 0 to 95, years from 1906 to 2009.}

\usage{data(ItalyM)}

\format{
  \code{obsq} is a \code{[1:96, 1:104]} numeric matrix containing probabilities of death. \cr
  \code{population} is a \code{[1:96, 1:104]} numeric matrix containing the male Italian population. \cr
  In both, row names are the ages \code{(0,1,\ldots,95)} and column names are the years \code{(1905,\ldots,2009)}
}

\source{Human Mortality Database \url{http://www.mortality.org}}

\references{
Human Mortality Database (2013). University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at \url{http://www.mortality.org} 
}

\seealso{
  \code{\link{DBKGrad-package}}, \code{\link{dbkGrad}}
}

\keyword{datasets}