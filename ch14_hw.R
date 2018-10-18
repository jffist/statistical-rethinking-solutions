rm(list=ls())
library(rethinking)
library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
#source("plot_bindings.R")


## Easy ####
# 14E1 ####
# For simplicity let's consider that observed log of the population is stored in the variable logpop_observed and its error in the variable logpop_se
# T_i ~ Poisson(mu_i)
# log(mu_i) = alpha + beta*logpop_estimated_i
# logpop_observed_i ~ Normal(logpop_estimated_i, logpop_se_i)
# alpha ~ Normal(0, 10)
# beta  ~ Normal(0, 1)


# 14E2 ####
# a) Below is the solution for the simpler task - allow missed values, but don't pay attention to observations errors
# T_i ~ Poisson(mu_i)
# log(mu_i) = alpha + beta*logpop_i
# logpop_i ~ Normal(nu, sigma_logpop) # the trick is in implementation to distinguish missed values from observed values used as data to estimated nu and sigma_logpop
# alpha ~ Normal(0, 10)
# beta  ~ Normal(0, 1)
# nu ~ Normal(0,10)
# sigma_logpop ~ dcauchy(0,1)

# b) Below is my best guess of how to describe a model that encounters both for observations errors and missed predictor values.
# T_i ~ Poisson(mu_i)
# log(mu_i) = alpha + beta*logpop_estimated_i
# logpop_observed_i ~ Normal(logpop_estimated_i, logpop_se_i) #as in 14e1
# logpop_estimated_i ~ Normal(nu, sigma_logpop)  # add distribution to it as in 14e2.a
# alpha ~ Normal(0, 10)
# beta  ~ Normal(0, 1)
# nu ~ Normal(0,10)
# sigma_logpop ~ dcauchy(0,1)


## Medium ####
# 14M1 ####
# I refer to the model m14.4 from the section 14.2.2. According to this model missed neocortex values are drawn from the Gaussian distribution with parameters `nu` and `sigma_N`.
# Parameter `nu` in its turn is unique for each example and is calculated as a linear model using logpopulation as a predictor. 
# Posterior distributions of intercept and slope for the linear model of  `nu` are estimated during the model inference.

# 14M2 ####
# I think it is not correct to compare WAIC of models fitted with neocortex imputation and model on complete cases only. 
# This model use data with different number of rows, so I would expect that WAIC of the model with imputation will be higher just because it has 12 more cases in the sum.
# Let's check it. I use code from the book to fit both models.
data(milk)
d <- milk
d$neocortex.prop <- d$neocortex.perc / 100
d$logmass <- log(d$mass)

# prep data
data_list <- list(
  kcal = d$kcal.per.g,
  neocortex = d$neocortex.prop,
  logmass = d$logmass )
# fit model with imputation
m14m2.with.imputation <- map2stan(
  alist(
    kcal ~ dnorm(mu,sigma),
    mu <- a + bN*neocortex + bM*logmass,
    neocortex ~ dnorm(nu,sigma_N),
    a ~ dnorm(0,100),
    c(bN,bM) ~ dnorm(0,10),
    nu ~ dnorm(0.5,1),
    sigma_N ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  ) ,
  data=data_list , iter=1e4 , chains=2 )
show(m14m2.with.imputation)
# Log-likelihood at expected values: 22.38 
# Deviance: -44.77 
# DIC: -29.09 
# Effective number of parameters (pD): 7.84 
# 
# WAIC (SE): -29.76 (5.5)
# pWAIC: 5.98 

# select only complete cases
dcc <- d[ complete.cases(d$neocortex.prop), ]
data_list_cc <- list(
  kcal = dcc$kcal.per.g,
  neocortex = dcc$neocortex.prop,
  logmass = dcc$logmass )
# fit model on complete cases only
m14m2.cc <- map2stan(
  alist(
    kcal ~ dnorm(mu,sigma),
    mu <- a + bN*neocortex + bM*logmass,
    a ~ dnorm(0,100),
    c(bN,bM) ~ dnorm(0,10),
    sigma ~ dcauchy(0,1)
  ) ,
  data=data_list_cc , iter=1e4 , chains=2 )
show(m14m2.cc)
# Log-likelihood at expected values: 12.13 
# Deviance: -24.27 
# DIC: -16.78 
# Effective number of parameters (pD): 3.74 
# 
# WAIC (SE): -16.98 (5)
# pWAIC: 3.04 

# Unexpectedly for me, the model with imputation has lower WAIC(so is better) disregard the fact that it incorporates larger dataset. 

# 14M3 ####
data(WaffleDivorce)
d <- WaffleDivorce

dlist <- list(
  div_obs=d$Divorce,
  div_sd=d$Divorce.SE,
  mar_obs=d$Marriage,
  mar_sd=d$Marriage.SE,
  A=d$MedianAgeMarriage )
m14m3.base <- map2stan(
  alist(
    div_est ~ dnorm(mu,sigma),
    mu <- a + bA*A + bR*mar_est[i],
    div_obs ~ dnorm(div_est, div_sd),
    mar_obs ~ dnorm(mar_est, mar_sd),
    a ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bR ~ dnorm(0,10),
    sigma ~ dcauchy(0,2.5)
  ) ,
  data=dlist ,
  start=list(div_est=dlist$div_obs,mar_est=dlist$mar_obs) ,
  WAIC=FALSE , iter=5000 , warmup=1000 , chains=3 , cores=3 ,
  control=list(adapt_delta=0.95) )
precis(m14m3.base, depth=2)
m14m3.base


dlist.se2 <- list(
  div_obs=d$Divorce,
  div_sd=2*d$Divorce.SE,
  mar_obs=d$Marriage,
  mar_sd=2*d$Marriage.SE,
  A=d$MedianAgeMarriage )
m14m3.se2 <- map2stan(m14m3.base, data=dlist.se2,
  start=list(div_est=dlist.se2$div_obs, mar_est=dlist.se2$mar_obs) ,
  WAIC=FALSE , iter=1e+4, warmup=5000 , chains=3, cores=3,
  control=list(adapt_delta=0.95) )
precis(m14m3.se2, depth=2)
# iter=5000 , warmup=1000
# Model didn't converge. It looks like we couldn't rely on inferred results because an effective number of samples is small and Rhat is very large.
# iter=1e+4, warmup=5000 
# Increasing number of iterations for warmup helped to reduce Rhat to 1 for almost all parameters, but n_eff is still very low (~300).
# From these two observations I conclude that with larger deviance of observations inference becomes harder.

## Hard ###
# 14H1 ####
data("elephants")
d <- elephants

m14h1.base <- map2stan(
  alist(
    MATINGS ~ dpois(lambda),
    log(lambda) ~ a + b*AGE,
    a ~ dnorm(0,10),
    b ~ dnorm(0,1)
  ),
  data=d, iter=4000, chains=2, cores=2
)
precis(m14h1.base)
plot(m14h1.base)

m14h1.se <- map2stan(
  alist(
    MATINGS ~ dpois(lambda),
    log(lambda) ~ a + b*AGE_est[i],
    AGE ~ dnorm(AGE_est, 5),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1)
  ),
  data=d, 
  start=list(AGE_est=d$AGE),
  iter=4000, chains=2, cores=2
)
precis(m14h1.se, depth=2)
#plot(m14h1.se)


matings_sample <- sim(m14h1.se)
matings_est <- apply(matings_sample, 2, mean)

post <- extract.samples(m14h1.se)
age_estimated <- apply(post$AGE_est, 2, mean)

plot(1, 1, xlab='age', ylab='mating', xlim=c(25, 55), ylim=c(0, 10),  type='n')
points(d$AGE, d$MATINGS, pch=16, col='blue')
points(age_estimated, matings_est)

for ( i in 1:nrow(d) ){
  lines( 
    c(d$AGE[i], age_estimated[i]), 
    c(d$MATINGS[i], matings_est[i]) 
  )
}

age_seq <- seq(25, 55, by=0.5)
lambda <- sapply(age_seq, function(x) exp(post$a + x*post$b))
lambda.avg <- apply(lambda, 2, mean)
lambda.PI <- apply(lambda, 2, PI)

lines(age_seq, lambda.avg, type='l', lty=2)
shade(lambda.PI, age_seq)

## 14H2 ####
m14h2 <- map2stan(
  alist(
    MATINGS ~ dpois(lambda),
    log(lambda) ~ a + b*AGE_est[i],
    AGE ~ dnorm(AGE_est, 40),
    a ~ dnorm(0,10),
    b ~ dnorm(0,1)
  ),
  data=d, 
  start=list(AGE_est=d$AGE),
  iter=4000, chains=2, cores=2
)
precis(m14h2, depth=2)
# When standard error for AGE observation equals 40, mean estimates of its slope in the model is 0.01 (vs. 0.07 when std=5).
# AGE estimates become negative or close to zero, thus unrealistic with such a big standard error
# Just for comparison AGE range is [27, 52] with mean 36

## 14H3 ####
set.seed(100)
x <- c( rnorm(10) , NA )
y <- c( rnorm(10,x) , 100 )
d <- list(x=x, y=y)

m14h3.cc <- map2stan(alist(
      y ~ dnorm(mu, sigma),
      mu <- a + b*x,
      a ~ dnorm(0, 100),
      b ~ dnorm(0, 1),
      sigma ~ dcauchy(0,1)
  ),
  data=list(x=x[1:10], y=y[1:10])
)
precis(m14h3.cc)
post14h3.cc <- extract.samples(m14h3.cc)
dens(post14h3.cc$b)
abline(v=0, lty=2)

m14h3.full <- map2stan(alist(
    y ~ dnorm(mu, sigma),
    mu <- a + b*x,
    x ~ dnorm(0, 1),
    a ~ dnorm(0, 100),
    b ~ dnorm(0, 100),
    sigma ~ dcauchy(0,1)
  ),
  data=d, iter=3000, chains=2, cores=2
)
precis(m14h3.full)
post14h3.full <- extract.samples(m14h3.full)
dens(post14h3.full$b)
abline(v=0, lty=2)
# fun plot :)

plot(d$x, d$y, ylim=c(-5,100), pch=16, col='black')
for(i in 1:1000){
  abline(a=post14h3.full$a[i], b=post14h3.full$b[i], col=col.alpha('blue'))
  points(post14h3.full$x_impute[i], 100, pch=16, col=col.alpha('red'))
}
points(d$x, d$y, ylim=c(-5,100), pch=16, col='black')
