pkgname <- "ipw"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ipw')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("basdat")
### * basdat

flush(stderr()); flush(stdout())

### Name: basdat
### Title: HIV: TB and Survival (Baseline Data)
### Aliases: basdat
### Keywords: datasets

### ** Examples

#see ?ipwtm for example



cleanEx()
nameEx("haartdat")
### * haartdat

flush(stderr()); flush(stdout())

### Name: haartdat
### Title: HAART and Survival in HIV Patients
### Aliases: haartdat
### Keywords: datasets

### ** Examples

#see ?ipwtm for example



cleanEx()
nameEx("healthdat")
### * healthdat

flush(stderr()); flush(stdout())

### Name: healthdat
### Title: IQ, Income and Health
### Aliases: healthdat
### Keywords: datasets

### ** Examples

#see ?ipwpoint for example



cleanEx()
nameEx("ipwplot")
### * ipwplot

flush(stderr()); flush(stdout())

### Name: ipwplot
### Title: Plot Inverse Probability Weights
### Aliases: ipwplot
### Keywords: hplot

### ** Examples

#see ?ipwpoint and ?ipwtm for examples



cleanEx()
nameEx("ipwpoint")
### * ipwpoint

flush(stderr()); flush(stdout())

### Name: ipwpoint
### Title: Estimate Inverse Probability Weights (Point Treatment)
### Aliases: ipwpoint
### Keywords: htest models

### ** Examples

#Simulate data with continuous confounder and outcome, binomial exposure.
#Marginal causal effect of exposure on outcome: 10.
n <- 1000
simdat <- data.frame(l = rnorm(n, 10, 5))
a.lin <- simdat$l - 10
pa <- exp(a.lin)/(1 + exp(a.lin))
simdat$a <- rbinom(n, 1, prob = pa)
simdat$y <- 10*simdat$a + 0.5*simdat$l + rnorm(n, -10, 5)
simdat[1:5,]

#Estimate ipw weights.
temp <- ipwpoint(
   exposure = a,
   family = "binomial",
   link = "logit",
   numerator = ~ 1,
   denominator = ~ l,
   data = simdat)
summary(temp$ipw.weights)

#Plot inverse probability weights
graphics.off()
ipwplot(weights = temp$ipw.weights, logscale = FALSE,
   main = "Stabilized weights", xlim = c(0, 8))

#Examine numerator and denominator models.
summary(temp$num.mod)
summary(temp$den.mod)

#Paste inverse probability weights
simdat$sw <- temp$ipw.weights

#Marginal structural model for the causal effect of a on y
#corrected for confounding by l using inverse probability weighting
#with robust standard error from the survey package.
msm <- (svyglm(y ~ a, design=svydesign(~1, weights=~sw,
   data=simdat)))
coef(msm)
confint(msm)

#Compute basic bootstrap confidence interval .
#boot.fun <- function(dat, index){
#    coef(glm(
#        formula = y ~ a,
#        data = dat[index,],
#        weights = ipwpoint(
#            exposure = a,
#            family = "gaussian",
#            numerator = ~ 1,
#            denominator = ~ l,
#            data = dat[index,])$ipw.weights))[2]
#    }
#bootres <- boot(simdat, boot.fun, 499);bootres
#boot.ci(bootres, type = "basic")




cleanEx()
nameEx("ipwtm")
### * ipwtm

flush(stderr()); flush(stdout())

### Name: ipwtm
### Title: Estimate Inverse Probability Weights (Time Varying)
### Aliases: ipwtm
### Keywords: htest models

### ** Examples

########################################################################
#EXAMPLE 1

#Load longitudinal data from HIV positive individuals.
data(haartdat)

#CD4 is confounder for the effect of initiation of HAART therapy on mortality.
#Estimate inverse probability weights to correct for confounding.
#Exposure allocation model is Cox proportional hazards model.
temp <- ipwtm(
   exposure = haartind,
   family = "survival",
   numerator = ~ sex + age,
   denominator = ~ sex + age + cd4.sqrt,
   id = patient,
   tstart = tstart,
   timevar = fuptime,
   type = "first",
   data = haartdat)

#plot inverse probability weights
graphics.off()
ipwplot(weights = temp$ipw.weights, timevar = haartdat$fuptime,
   binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability weights")

#CD4 count has an effect both on dropout and mortality, which causes informative censoring.
#Use inverse probability of censoring weighting to correct for effect of CD4 on dropout.
#Use Cox proportional hazards model for dropout.
temp2 <- ipwtm(
   exposure = dropout,
   family = "survival",
   numerator = ~ sex + age,
   denominator = ~ sex + age + cd4.sqrt,
   id = patient,
   tstart = tstart,
   timevar = fuptime,
   type = "first",
   data = haartdat)

#plot inverse probability of censoring weights
graphics.off()
ipwplot(weights = temp2$ipw.weights, timevar = haartdat$fuptime,
   binwidth = 100, ylim = c(-1.5, 1.5), main = "Stabilized inverse probability of censoring weights")

#MSM for the causal effect of initiation of HAART on mortality.
#Corrected both for confounding and informative censoring.
#With robust standard error obtained using cluster().
summary(coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient),
   data = haartdat, weights = temp$ipw.weights*temp2$ipw.weights))   

#uncorrected model
summary(coxph(Surv(tstart, fuptime, event) ~ haartind, data = haartdat))

########################################################################
#EXAMPLE 2

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
#   #taking into account the uncertainty introduced by modelling longitudinal CD4
#   #taking into account the uncertainty introduced by estimating the inverse probability weights
#   #robust with regard to weights unequal to 1
#   boot.fun <- function(data, index, data.tm){
#      data.samp <- data[index,]
#      data.samp$id.samp <- 1:nrow(data.samp)
#      data.tm.samp <- do.call("rbind", lapply(data.samp$id.samp, function(id.samp)cbind(data.tm[data.tm$id == data.samp$id[data.samp$id.samp == id.samp],], id.samp = id.samp)))
#      cd4.lme <- lme(cd4.sqrt ~ fuptime + tb.lag, random = ~ fuptime | id.samp, data = data.tm.samp)
#      times <- sort(unique(c(data.samp$Ttb, data.samp$Tend)))
#      startstop.samp <- data.frame(id.samp = rep(data.samp$id.samp, each = length(times)), fuptime = rep(times, nrow(data.samp)))
#      startstop.samp <- merge(startstop.samp, data.samp, by = "id.samp", all.x = TRUE)
#      startstop.samp <- startstop.samp[with(startstop.samp, fuptime <= Tend),]
#      startstop.samp$tstart <- tstartfun(id.samp, fuptime, startstop.samp)
#      startstop.samp$tb <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime >= Ttb), 1, 0)
#      startstop.samp$tb.lag <- ifelse(with(startstop.samp, !is.na(Ttb) & fuptime > Ttb), 1, 0)
#      startstop.samp$event <- ifelse(with(startstop.samp, !is.na(Tdeath) & fuptime >= Tdeath), 1, 0)
#      startstop.samp$cd4.sqrt <- predict(cd4.lme, newdata = data.frame(id.samp = startstop.samp$id.samp, fuptime = startstop.samp$fuptime, tb.lag = startstop.samp$tb.lag))
#      return(coef(coxph(Surv(tstart, fuptime, event) ~ tb, data = startstop.samp,
#         weights = ipwtm(
#              exposure = tb,
#              family = "survival",
#              numerator = ~ 1,
#              denominator = ~ cd4.sqrt,
#              id = id.samp,
#              tstart = tstart,
#              timevar = fuptime,
#              type = "first",
#              data = startstop.samp)$ipw.weights))[1])
#      }
#   bootres <- boot(data = basdat, statistic = boot.fun, R = 999, data.tm = timedat);bootres
#   boot.ci(bootres, type = "basic")



cleanEx()
nameEx("timedat")
### * timedat

flush(stderr()); flush(stdout())

### Name: timedat
### Title: HIV: TB and Survival (Longitudinal Measurements)
### Aliases: timedat
### Keywords: datasets

### ** Examples

#See ?ipwtm for example



cleanEx()
nameEx("tstartfun")
### * tstartfun

flush(stderr()); flush(stdout())

### Name: tstartfun
### Title: Compute Starting Time For Counting Process Notation
### Aliases: tstartfun
### Keywords: methods survival

### ** Examples

#data
mydata1 <- data.frame(
   patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
   time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35))

#compute starting time for each interval
mydata1$tstart <- tstartfun(patient, time.days, mydata1)

#result
mydata1

#see also ?ipwtm for example



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
