# About Survival Analysis
library(OIsurv) # Includes the "survival" and "KMsurv" packages used for analysis and data sets
#other packages in the market, this isn't the only one for survival analysis

# What is survival analysis? - A set of methods for analyzing data where the outcome variable is the time until the occurrence of an even of interest, e.g. death.

# Why not linear regression? - Survival times are always positive and regression canno't handle the censoring of observations, e.g. in a given study if some of the patients survive past when the data is collected, these patient observations represent a right censor.  Another cause of censoring is from patients dropping out of the study.  Unlike regression models, survival models correctly incorporate information from oth censored and uncensored observations.

# In survival analysis we can estimate two functions dependent on time:
# [1] The survival function - gives, for every time, the probability of surviving (not experiencing the event)
# [2] The hazard function - gives the potential that the event will occur, per time unit, given that an individual has surived up to the specific time.

# Functions in the survival packages apply methods to Surv objects, which are created by the Surv() function.

# Censoring
library(OIsurv)
# Here's a dataset that looks at survival times for individuals with a certain type of tumor
data(tongue)
attach(tongue)
# Let's start by looking at group 1 only:
g1 <- Surv(time[type==1],delta[type==1])
#type: only look at type 1 tongue cancer
g1
# shows us an ordered list of survival times, plus signs represent observations that are right censored.
detach(tongue)

# Here's an example of left-truncated right-censored observations:
data(psych)
p1 <- with(psych, Surv(age,age+time, death)) # note I have to use the with function here because I did not attach psych
#age + time = the age when the death is measured, either dead or still alive
p1
# Interpretation for first observation: Patient entered study at 51 years of age and survived until 52 years old.

# Estimating the Survival Function with Kaplan-Meier and Pointwise Confidence Intervals
library(OIsurv)
data(tongue)
g1 <- with(tongue, Surv(time[type==1],delta[type==1]))

# The Kaplan-Meier estimate is a nonparametric MLE of the survival function, S(t)

# Fitting a survival function like you would a regression...
# Here we use the simplest model where we look at the survival object against an intercept.
fit <- survfit(g1~1, data = tongue, conf.int = .95, conf.type = "log") # for 95% confidence interval with interval type being a log function (could be linear with "plain" or could be log(-log(t))) with "log-log"
fit
summary(fit) # survival = Kaplan Meier estimate at each time
plot(fit, main = "Kaplan-Meier estimate with 95% point-wise confidence", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 200))
#survival probability is plotted, which also is the 4th column in the summary(fit)
# shows us the survival probability for each week.  The confidence intervals are valid only pointwise; the confidence range does not capture with 95% confidence over the entire range of time values, but only the confidence range for a particular time value.

# we can also split a Kaplan Meier estimate across a specific variable, e.g.:
g2 <- with(tongue, Surv(time, delta))
#another survival subject: not only look at type1, but look at type 1 and 2.
fit2 <- survfit(g2~type, data = tongue, conf.int = .95, conf.type = "log")
summary(fit2)
plot(fit2, main = "Kaplan-Meier estimates", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 200), lty=c(1,2))
legend('topright', c("Type1","Type2"), lty=c(1,2))
#the question is: are these 2 basically the same or significantly different?

# Comparing Two Survival Curves
# Let's do a test to see if the two survival curves above are statistically different.
survdiff(Surv(time, delta)~type, data = tongue)
# We reject the null that both the survival functions are the same at the 90% confidence level; however, we fail to reject the null at the 95% level.
#p= 0.0949 , with 95% ci, they are basically the same.
#might be different result with 90% ci. 

# Simultaneous Confidence Intervals
library(OIsurv)
data(tongue)
g3 <- with(tongue, Surv(time[type==1],delta[type==1]))

# If we'd like to create confidence bands that capture the true survival function with a 95% accuracy we will need to use simultaneous confidence intervals.  This can be done with the confBands() function.
ci <- confBands(g3, confLevel = .95, confType = "log-log", type = "hall")
#confband!!!

fit3 <- survfit(g3~1, data = tongue, conf.int = .95, conf.type = "log-log")
plot(fit3, main = "Kaplan-Meier estimate with 95% point-wise confidence", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 200))
lines(ci, lty = 3, col = "red")
legend('topright', c("Survival Estimate","Pointwise Interval", "Simultaneous Interval"), lty=c(1,2,3), col = c("black", "Black", "red"))
#when it close to end, they started to expand, because there were less people left. 

# Cumulative Hazard Function
library(OIsurv)
data(tongue)
g4 <- with(tongue, Surv(time[type==1],delta[type==1]))
fit4 <- summary(survfit(g4~1, data = tongue, conf.int = .95, conf.type = "log-log"))
# The cumulative hazard function (H(t)) and the survival function S(t) are related in the following way for continuous data:
# S(t) = exp[-H(t)]

# Lets use our survival function to calculate estimates for the hazard function (potential particular event will occur):
H <- -log(fit4$surv)
H <- c(H, tail(H,1))
plot(c(fit4$time, 200), H, main = "Cumulative Hazard Functions", xlab = "Time (weeks)", ylab = "Hazard Functions", lty = 1, type = "s", ylim = range(H))
# By realizing H(t) = f(t)/S(t) H(t) can be interpreted as, "the density of events at t, divided by the probability of surviving to that duration without experiencing the event".  Essentially it's a ratio that measures how likely the event will occur in a standardized form.

# Another approximation of the cumulative hazard function is sum[(the number of individuals at risk)/(the number of events that took place after time, t)]:
H.2 <- cumsum(fit4$n.event / fit4$n.risk)
H.2 <- c(H.2, tail(H.2,1))
points(c(fit4$time, 200), H.2, lty = 2, type = "s")
legend("topleft", c("H","H.2"), lty = c(1,2))