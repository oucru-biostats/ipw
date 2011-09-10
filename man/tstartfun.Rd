\name{tstartfun}
\Rdversion{1.1}
\alias{tstartfun}

\title{Compute Starting Time For Counting Process Notation}

\description{Function to compute starting time for intervals of follow-up, when using the counting process notation. Within each unit under observation (usually individuals), computes starting time equal to:
\itemize{
\item time of previous record when there is a previous record.
\item -1 for first record.
}}


\usage{
tstartfun(id, timevar, data)
}

\arguments{
\item{id}{numerical vector, uniquely identifying the units under observation, within which the longitudinal measurements are taken.}
\item{timevar}{numerical vector, representing follow-up time, starting at 0.}
\item{data}{dataframe containing \code{id} and \code{timevar}.}
}


\value{Numerical vector containing starting time for each record. In the same order as the records in \code{data}, to facilitate merging.}

\section{Missing values}{Currently, \code{id} and \code{timevar} should not contain missing values.}

\references{
Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{http://www.jstatsoft.org/v43/i13/}.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#data
mydata1 <- data.frame(
   patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
   time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35))

#compute starting time for each interval
mydata1$tstart <- tstartfun(patient, time.days, mydata1)

#result
mydata1

#see also ?ipwtm for example
}

\keyword{methods}
\keyword{survival}

