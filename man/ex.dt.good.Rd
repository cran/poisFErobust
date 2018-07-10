\name{ex.dt.good}
\alias{ex.dt.good}
\docType{data}
\title{
Poisson data satisfying conditional mean assumption
}
\description{
A data.table containing \code{id} by \code{day} observations of Poisson
random variables which satisfy the conditional mean assumption of
Wooldridge (1999).
}
\usage{data("ex.dt.good")}
\format{
  A data frame with 500 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22} \code{23} \code{24} \code{25} \code{26} \code{27} \code{28} \code{29} \code{30} \code{31} \code{32} \code{33} \code{34} \code{35} \code{36} \code{37} \code{38} \code{39} \code{40} \code{41} \code{42} \code{43} \code{44} \code{45} \code{46} \code{47} \code{48} \code{49} \code{50}}
    \item{\code{day}}{a numeric vector}
    \item{\code{fe}}{a numeric vector}
    \item{\code{x1}}{a numeric vector}
    \item{\code{x2}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
The data were simulated like 
\code{y <- rpois(1, exp(fe + x1 + x2))}
where \code{fe}, \code{x1}, and \code{x2} are standard normal random variables.
\code{fe} varies only across \code{id}.
}
\references{
Wooldridge, Jeffrey M. (1999): "Distribution-free estimation of some nonlinear
panel data models," Journal of Econometrics, 90, 77-97.
}
\examples{
data("ex.dt.good")
str(ex.dt.good)
}
\keyword{datasets}
