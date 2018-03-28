library(ggfortify) # autoplot
library(dplyr)
library(OIsurv) #survival analysis tools
data(veteran)
head(veteran)

# Kaplan Analysis Again
# Analysis ignoring the treatment:
ka <- with(veteran, Surv(time, status))
ka
fit1 <- survfit(Surv(time, status) ~ 1, data=veteran)
summary(fit1)
plot(fit1, main = "Kaplan-Meier estimates", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 500))

# Analysis including the treatment:
fit2 <- survfit(Surv(time, status) ~ trt, data=veteran)
plot(fit2, main = "Kaplan-Meier estimates", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 500), lty=c(1,2))
legend('topright', c("Treat1","Treat2"), lty=c(1,2))

# Analysis looking across Age instead of treatment:
vet <- mutate(veteran, agefac = ifelse((age < 60), 1, 0))
fit3 <- survfit(Surv(time, status) ~ agefac, data=vet)
plot(fit3, main = "Kaplan-Meier estimates", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 500), lty=c(1,2))
legend('topright', c("Old","Young"), lty=c(1,2))
# Younger patients clearly have a better chance of survival.

# Cox Hazard Model
# Cox proportional hazard model that uses all covariates in the data set.
vet <- mutate(veteran, agefac = ifelse((age < 60), 1, 0), trt = factor(trt,labels=c("standard","test")), prior = factor(prior,labels=c("N0","Yes")))

cox <- coxph(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior , data = vet)
summary(cox)
fit4 <- survfit(cox)
plot(fit4, main = "Kaplan-Meier estimates", xlab = "Time (weeks)", ylab = "Survival Function", xlim = c(0, 500))
lines(fit1, col = "red")
legend('topright', c("Cox","kaplan"), col = c("black", "red"), lty=c(1,1))

# General notes:

# This model is estimated through a multiple linear estiamte of bi's and a ML estimate of the base line hazard h0(t):
# h(t,X) = h0(t)exp(Sum(biXi))

# Where the baseline hazard h0(t) is unspecified and needs to be estimated, also not it depends on t only and proportional hazard can be extracted through division.  The baseline hazard function is analogous to the intercept term in a multiple regression or logistic regression model.

# We can use the Cox hazard model in place of a logistic model where truncation exists. Notice if we rewrite the above equation we can write it as... ln(h(t,X)/h0(t)) = sum(biXi) which looks like a logit model

# Hazard ration can be directly computed as baseleine hazard will drop out: HR = h0(t)exp(Sum(biXi*)) / h0(t)exp(Sum(biXi^)) = exp(Sum(biXi)) / exp(Sum(biXi^)) = exp(Sum(bi(Xi*-Xi^)))

# Important assumption: Cox model assumes that the covariates do not vary with time.

# Because censorship is present R^2 is not appropriate, so we can use the "Concordance", which is simply the proportion of pairs of cases in which the case with the higher-risk predictor had an event before the case with the lower-risk predictor.  It ranges from 0 - 1, with 1 being perfect concordance, .5 being the expected result from random predictions, and 0 being anti-concordance.

# A positive coefficient implies the hazard is higher, and thus more likely for the subject with higher values of that variable.

# Aalen's additive regression model as a supplement to Cox
# Shows us how the effects of our covariates change over time.  We dont wan't abrupt changes in slope for these variables.
fit5 <- aareg(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior , data = vet)
autoplot(fit5)  #yay to autoplot
# Looking at karno, we can see that this variable may be an issue, due to the drastic/ obvious change in its slope.

# Other Distributions
# Cox PH is a special case of "Accelerated failure time models"
# We could also look at different distrubutions.  One example (among many distributions) is:
fit6 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist = "weibull")
summary(fit6) ##not a good model, because the p value is way way way to high. 
# This is likely not a good model as pho indicates not significant
