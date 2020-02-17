\name{pois.fe.robust}
\alias{pois.fe.robust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust standard errors of Poisson fixed effects regression
}
\description{
Compute standard errors following Wooldridge (1999) for Poisson regression with
fixed effects, and a hypothesis test of the conditional mean assumption (3.1).
}
\usage{
pois.fe.robust(outcome, xvars, group.name, data, 
               qcmle.coefs = NULL, allow.set.key = FALSE,
               index.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outcome}{
character string of the name of the dependent variable.
}
  \item{xvars}{
vector of character strings of the names of the independent variables.
}
  \item{group.name}{
character string of the name of the grouping variable.
}
  \item{data}{
data.table which contains the variables named in other arguments. See details 
for variable type requirements.
}
  \item{qcmle.coefs}{
an optional numeric vector of coefficients in the same order as \code{xvars}. 
If \code{NULL}, coefficients are estimated using \code{glmmML::glmmboot}.
}
  \item{allow.set.key}{logical. When \code{TRUE} (recommended),
  \code{data} will have its key set to
  \code{group.name}, so it may be reordered. This should reduce memory usage.
}
  \item{index.name}{
DEPRECATED (leave as NULL).
}
}
\details{
\code{data} must be a \code{data.table} containing the following:

\itemize{
\item{a column named by \code{outcome}, }{non-negative integer}
\item{columns named according to each string in \code{xvars}, }{numeric type}
\item{a column named by \code{group.name}, }{factor type}
\item{a column named by \code{index.name}, }{integer sequence increasing by one
each observation with no gaps within groups}
}

No observation in \code{data} may contain a missing value.

Setting \code{allow.set.key} to \code{TRUE} is recommended to reduce
memory usage; however, it will allow \code{data} to be modified
(sorted in-place).

\code{pois.fe.robust} also returns the p-value of the hypothesis test of the
conditional mean assumption (3.1) as described in Wooldridge (1999) section 3.3.
}
\value{
A list containing
\itemize{
\item{\code{coefficients}, }{a numeric vector of coefficients.}
\item{\code{se.robust}, }{a numeric vector of standard errors.}
\item{\code{p.value}, }{the p-value of a hypothesis test of the 
conditional mean assumption (3.1).}
}

}
\references{
Wooldridge, Jeffrey M. (1999): "Distribution-free estimation of some nonlinear
panel data models," Journal of Econometrics, 90, 77-97.
}
\author{
Evan Wright
}
\seealso{
\code{\link[glmmML:glmmboot]{glmmboot}}
}
\examples{
# ex.dt.good satisfies the conditional mean assumption
data("ex.dt.good")
pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), group.name = "id",
               index.name = "day", data = ex.dt.good)
               
# ex.dt.bad violates the conditional mean assumption
data("ex.dt.bad")
pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), group.name = "id",
               index.name = "day", data = ex.dt.bad)
}

\keyword{ regression}% use one of  RShowDoc("KEYWORDS")

