\name{healthdat}
\Rdversion{1.1}
\alias{healthdat}
\docType{data}
\title{IQ, Income and Health}

\description{A simulated dataset containing IQ, income and health score measurements in 1000 individuals.}

\usage{data(healthdat)}

\format{
A data frame with 1000 rows, with each row corresponding to a separate individual. The following variables are included:
\describe{
\item{\code{id}}{individual ID.}
\item{\code{iq}}{IQ score.}
\item{\code{income}}{gross monthly income (EUR).}
\item{\code{health}}{health score (0-100).}
}}

\details{
In these simulated data, IQ is a confounder for the effect of income on health.
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{healthdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#load IQ, income and health data.
data(healthdat)

#IQ is confounder for the effect of income on health score.
#Estimate inverse probability weights to correct for confounding.
#Exposure allocation model is linear regression model.
temp <- ipwpoint(
   exposure = income,
   family = "gaussian",
   numerator = ~ 1,
   denominator = ~ iq,
   data = healthdat)

#plot inverse probability weights
graphics.off()
ipwplot(weights = temp$ipw.weights, logscale = FALSE,
   main = "Stabilized weights", xlim = c(0, 8))

#Marginal structural model for the causal effect of income on health score
#corrected for confounding by IQ using inverse probability weighting.
summary(glm(health ~ income, data = healthdat, weights = temp$ipw.weights))

#Compute basic bootstrap confidence interval for income parameter.
boot.fun <- function(dat, index){
   coef(glm(
      formula = health ~ income,
      data = dat[index,],
      weights = ipwpoint(
            exposure = income,
            family = "gaussian",
            numerator = ~ 1,
            denominator = ~ iq,
            data = dat[index,])$ipw.weights))[2]
   }
bootres <- boot(healthdat, boot.fun, 499);bootres
boot.ci(bootres, type = "basic")
}

\keyword{datasets}
