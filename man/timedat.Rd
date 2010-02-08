\name{timedat}
\Rdversion{1.1}
\alias{timedat}
\docType{data}
\title{HIV: TB and Survival (Longitudinal Measurements)}

\description{Simulated dataset. Time varying CD4 measurements of 386 HIV positive individuals. Time of first active tuberculosis, time of death and individual end time of the patients are included in dataset \code{\link{basdat}}.}

\usage{data(timedat)}
\format{
	A data frame with 6291 observations on the following 3 variables.
	\describe{
		\item{\code{id}}{patient ID.}
		\item{\code{fuptime}}{follow-up time (days since HIV seroconversion).}
		\item{\code{cd4count}}{CD4 count measured at fuptime.}
	}
}
\details{
These simulated data are used together with data in \code{\link{basdat}} in a detailed causal modelling example using inverse probability weighting (IPW). See below for the example. Data were simulated using the algorithm described in Van der Wal e.a. (2009).
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Cole, S. R. & Hernán, M. A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \emph{168}(6), 656-664.

Robins, J. M., Hernán, M. A. & Brumback, B. A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \emph{11}, 550-560.

Van der Wal W. M., Prins M., Lumbreras B. & Geskus R. B. (2009). A simple G-computation algorithm to quantify the causal effect of a secondary illness on the progression of a chronic disease. \emph{Statistics in Medicine}, \emph{28}(18), 2325-2337.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{healthdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
data(basdat)
data(timedat)

#Aim: to model the causal effect of active tuberculosis (TB) on mortality.
#Longitudinal CD4 is a confounder as well as intermediate for the effect of TB.

#process original measurements
   #check for ties (not allowed)
      table(duplicated(timedat[,c("id", "fuptime")]))
   #take square root of CD4 because of skewness
      timedat$cd4.sqrt <- sqrt(timedat$cd4count)
   #add TB time to dataframe
      timedat <- merge(timedat, basdat[,c("id", "Ttb")], by = "id", all.x = TRUE)
   #compute TB status
      timedat$tb.lag <- ifelse(with(timedat, !is.na(Ttb) & fuptime > Ttb), 1, 0)
   #longitudinal CD4-model
      cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id,
      data = timedat)

#build new dataset:
#rows corresponding to TB-status switches, and individual end times
   times <- sort(unique(c(basdat$Ttb, basdat$Tend)))
   startstop <- data.frame(
      id = rep(basdat$id, each = length(times)),
      fuptime = rep(times, nrow(basdat)))
   #add baseline data to dataframe
      startstop <- merge(startstop, basdat, by = "id", all.x = TRUE)
   #limit individual follow-up using Tend
      startstop <- startstop[with(startstop, fuptime <= Tend),]
   startstop$tstart <- tstartfun(id, fuptime, startstop) #compute tstart (?tstartfun)
   #indicate TB status
      startstop$tb <- ifelse(with(startstop, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
   #indicate TB status at previous time point
      startstop$tb.lag <- ifelse(with(startstop, !is.na(Ttb) & fuptime > Ttb), 1, 0)
   #indicate death
      startstop$event <- ifelse(with(startstop, !is.na(Tdeath) & fuptime >= Tdeath),
      1, 0)
   #impute CD4, based on TB status at previous time point.
      startstop$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id = startstop$id,
         fuptime = startstop$fuptime, tb.lag = startstop$tb.lag))

#compute inverse probability weights
   temp <- ipwtm(
      exposure = tb,
      family = "survival",
      numerator = ~ 1,
      denominator = ~ cd4.sqrt,
      id = id,
      tstart = tstart,
      timevar = fuptime,
      type = "first",
      data = startstop)
   summary(temp$ipw.weights)
   ipwplot(weights = temp$ipw.weights, timevar = startstop$fuptime, binwidth = 100)

#models
   #IPW-fitted MSM, using cluster() to obtain robust standard error estimate
      summary(coxph(Surv(tstart, fuptime, event) ~ tb + cluster(id),
      data = startstop, weights = temp$ipw.weights))
   #unadjusted
      summary(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop))
   #adjusted using conditioning: part of the effect of TB is adjusted away
      summary(coxph(Surv(tstart, fuptime, event) ~ tb + cd4.sqrt, data = startstop))

#compute bootstrap CI for TB parameter (takes a few hours)
#   #robust with regard to weights unequal to 1
#   boot.fun <- function(data, index, data.tm){
#      data.samp <- data[index,]
#      data.samp$id.samp <- 1:nrow(data.samp)
#      data.tm.samp <- do.call("rbind", lapply(data.samp$id.samp, function(id.samp)
#         cbind(data.tm[data.tm$id == data.samp$id[data.samp$id.samp == id.samp],],
#         id.samp = id.samp)))
#      cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id.samp,
#         data = data.tm.samp)
#      times <- sort(unique(c(data.samp$Ttb, data.samp$Tend)))
#      startstop.samp <- data.frame(id.samp = rep(data.samp$id.samp,
#         each = length(times)), fuptime = rep(times, nrow(data.samp)))
#      startstop.samp <- merge(startstop.samp, data.samp, by = "id.samp",
#          all.x = TRUE)
#      startstop.samp <- startstop.samp[with(startstop.samp, fuptime <= Tend),]
#      startstop.samp$tstart <- tstartfun(id.samp, fuptime, startstop.samp)
#      startstop.samp$tb <- ifelse(with(startstop.samp, !is.na(Ttb) &
#          fuptime >= Ttb), 1, 0)
#      startstop.samp$tb.lag <- ifelse(with(startstop.samp, !is.na(Ttb) &
#          fuptime > Ttb), 1, 0)
#      startstop.samp$event <- ifelse(with(startstop.samp, !is.na(Tdeath) &
#          fuptime >= Tdeath), 1, 0)
#      startstop.samp$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id.samp =
#         startstop.samp$id.samp, fuptime = startstop.samp$fuptime, tb.lag =
#         startstop.samp$tb.lag))
#      return(coef(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop.samp,
#         weights = ipwtm(
#            exposure = tb,
#            family = "survival",
#            numerator = ~ 1,
#            denominator = ~ cd4.sqrt,
#            id = id.samp,
#            tstart = tstart,
#            timevar = fuptime,
#            type = "first",
#            data = startstop.samp)$ipw.weights))[1])
#      }
#   bootres <- boot(data = basdat, statistic = boot.fun, R = 999,
#      data.tm = timedat);bootres
#   boot.ci(bootres, type = "basic")
}
\keyword{datasets}
