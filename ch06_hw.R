library(rethinking)
library(ggplot2)
library(GGally)
library(dplyr)

#source("plot_bindings.R")


### Easy ####
# 6E1
# illustration of entropy additivity for independent events
e1 <- c(0.3, 0.7) #coin 1
e2 <- c(0.6, 0.4) #coin 2
entropy <- function(x){
  x <- x[x!=0]
  -sum(x*log(x))
}
entropy(e1)
entropy(e2)
e12 <- c(0.3*0.6, 0.3*0.4, 0.7*0.6, 0.7*0.4) #joint distribution of coin1 and coin2
entropy(e12) - (entropy(e1)+entropy(e2)) #difference of joint and sum of independent is equal up to float compariosn precision

# 6E2
entropy(c(0.7, 0.3))

# 6E3
entropy(c(0.2, 0.25, 0.25, 0.3))
# as expected more choices gives higher entropy (compared to coin from 6E2)

# 6E4
entropy(c(1/3, 1/3, 1/3, 0))

### Medium ####
# 6M1
# AIC = D_train + 2*n_params 
#  assumptions:
#     * flat priors (or overwhelmed by likelihood - e.g. big sample size, small number of params)
#     * posterior is multivariative gaussian
#     * N_samples >> n_params
# DIC = mean(D) + p_d = mean(D) + (mean(D) - D_at_mean)
#     * D - sampled distribution of deviance
#     * p_d - effective number of params
#     * D_at_mean - deviance calculated at the posterior mean
#  assumptions:
#     * posterior is multivariative gaussian
# WAIC = -2(lppd - p_waic)
#     * lppd - log pointwise predictive density
#        * lppd = sum(log(Pr(y_i))) i=1:N_samples //sum through all samples
#        * Pr(y_i) - mean likelihood of observation y_i from sample. 
#                    As I understand, it's caclulated by averaging likelihood from 
#                    posterior distribution of parameters
#     * p_waic - effective number of parameters
#       * p_waic = sum(Var_ll(y_i))  i=1:N_samples //sum through all samples
#       * Var_ll(y_i) - variance of the likelihood of observation y_i. 
#                       As with Pr(y_i), Var_ll is calculated from posterior params distribution
#     * measure is pointwise
#     * application is questionable for cases where observations depends on each other as in timeseries
# Answer:
#   - From most to less general: WAIC > DIC > AIC
#   - assumptions: no > Gaussian posterior > Gaussian posterior + flat priors
#   - when the posterior predictive mean is a good representation of the posterior predictive distribution, WAIC and DIC will tend to agree, as a  consequence from the Gaussian assumptions on posterior distribution.
#   - when priors are effectively flat or overwhelmed by the amount of data, the DIC and AIC will tend to agree

# 6M2. 
#Explain the difference between model selection and model averaging. 
# What information is lost under model selection? What information is lost under model averaging?
#
# Selection - choose one best model for prediction. We lose information about differences among models. 
#   For example, if m1 and m2 were very close to each other compared on a performance metric but different in structure. 
#   Then selecting only one of them discard part of the evidence.
# Averaging - include uncertainty about 'correctness' of the models into predictions. 
#   Technically speaking, information criteria is used as a weight of certainty of the model and predictions of different models 
#   are mixed using the weights. 
#   Averaging discard information about model differences.


# 6M3
# All information criteria use sum over all samples as part of the calculation(is an additive measure of an error)
# Different number of examples (or different examples in the sample) will cause differences in the value of the criteria that influenced. Models fitted on a smaller sample will have better results (value of the criteria is lower).

# 6M4
# Effective number of parameters should reduce, as models become less 'flexible' and more rigid
# let's check
data(milk)
d <- milk[ complete.cases(milk) , ]
d$neocortex <- d$neocortex.perc / b.sd
str(d)

a.start <- mean(d$kcal.per.g) 
sigma.start <- sd(d$kcal.per.g)
m6m4.1 <- map(alist(
    kcal.per.g ~ dnorm( mu , sigma),
    mu <- Intercept + b_neocortex*neocortex,
    Intercept ~ dnorm(0,10),
    b_neocortex ~ dnorm(0,10),
    sigma ~ dcauchy(0,2)
  ), 
  data=d,
  start=list(Intercept=a.start,b_neocortex=0, sigma=sigma.start))
precis(m6m4.1)
WAIC(m6m4.1)#pWAIC = 3.874535

m6m4.2 <- map(alist(
    kcal.per.g ~ dnorm( mu , sigma),
    mu <- Intercept + b_neocortex*neocortex,
    Intercept ~ dnorm(0,10),
    b_neocortex ~ dnorm(0,0.2),#!!: changed sigma of parameter
    sigma ~ dcauchy(0,2)
  ), 
  data=d,
  start=list(Intercept=a.start,b_neocortex=0, sigma=sigma.start))
precis(m6m4.2)
WAIC(m6m4.2)#pWAIC = 2.061632 (smaller as expected)

# 6M5
# Informative priors reduce overfitting because more data is required to shift params from their initial most probable values. It means that if a pattern is not frequent in data, it will be hard to move params from the initial state, as data has a small influence.  
# Another possible explanation - regularised(more concentrated) priors are equal to adding a big set of initial data concentrated around parameters with priors.
# Excerpt from the book:
#   One way to prevent a model from getting too excited by the training sample is to give it a skeptical prior. By “skeptical,” I mean a prior that slows the rate of learning from the sample. The most common skeptical prior is a regularizing prior, which is applied to a beta-coefficient, a “slope” in a linear model. Such a prior, when tuned properly, reduces overfitting while still allowing the model to learn the regular features of a sample. If the prior is too skeptical, however, then regular features will be missed, resulting in underfitting. So the problem is really one of tuning. But as you’ll see, even mild skepticism can help a model do better, and doing better is all we can really hope for in the large world, where no model nor prior is optimal

# 6M6
# Overly informative priors result in underfitting because effective number of parameters reduces and model have less freedom to fit the data. If regularisation is too strong fitting requires more and more data, and model underfits on the small samples.
# In other words, a model becomes overconfident in parameters values and it's hard to move it from that point. It requires more evidence(data) for even tiny changes of parameters.

### Hard ####
library(rethinking)
data(Howell1)
d <- Howell1
d$age <- (d$age - mean(d$age))/sd(d$age)
set.seed( 1000 )
i <- sample(1:nrow(d),size=nrow(d)/2)
d1 <- d[ i , ]
d2 <- d[ -i , ]

formulas <- list()
sigma.rb <- 20
b.sd <- 100
a.sd <- 100


formulas[[1]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
formulas[[2]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
formulas[[3]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2 + b3*age^3,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  b3 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
formulas[[4]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2 + b3*age^3 + b4*age^4,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  b3 ~ dnorm(0, b.sd),
  b4 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
formulas[[5]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2 + b3*age^3 + b4*age^4 + b5*age^5,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  b3 ~ dnorm(0, b.sd),
  b4 ~ dnorm(0, b.sd),
  b5 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
formulas[[6]] <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2 + b3*age^3 + b4*age^4 + b5*age^5 + b6*age^6,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  b3 ~ dnorm(0, b.sd),
  b4 ~ dnorm(0, b.sd),
  b5 ~ dnorm(0, b.sd),
  b6 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)

models <- list()
for(i in 1:6){
  print(paste("fitting model ",i))
  models[[i]] <- map(formulas[[i]], data=d1)
}

### 6H1
models.cmp <- compare(models[[1]], models[[2]], models[[3]],
                      models[[4]], models[[5]], models[[6]]) 
models.cmp
plot(models.cmp)

#### 6H2
age.seq <- seq(-2, 3.5, by=0.1)
par(mfrow=c(3,3))
plot.predictions <- TRUE
d4p <- d1
d.predict <- data.frame(age=age.seq)

plot_model <- function(name, m, d.predict, d.raw, plot.predictions.pi=TRUE){
  mu <- link(m, data=d.predict)
  mu.avg <- apply(mu, 2, mean)
  mu.pi <- apply(mu, 2, PI, .97)
  
  # calculate data for points predictions
  if(plot.predictions.pi){
    predictions <- sim(m, data=d.predict)
    #predictions.avg <- apply(predictions, 2, mean)
    predictions.pi <- apply(predictions, 2, PI, .97) 
  }
  
  plot(height ~ age, d.raw, xlim=c(-2,3.5), ylim=c(0, 200), col=col.alpha('slateblue', 0.5))
  lines(d.predict$age, mu.avg, col='red')
  shade(mu.pi, d.predict$age, col = col.alpha("blue",0.15), border=1)
  mtext(paste('#params=',name))
  if(plot.predictions){
    #lines(age.seq,predictions.avg)
    shade(predictions.pi, age.seq)
  }
}
for(i in 1:6) {
  m <- models[[i]]

  plot_model(i, m, d.predict, d4p)
}

# 6H3
par(mfrow=c(2,1))
#### model 4
plot_model(4, models[[4]], d.predict, d4p)

#### ensemble
m.ensemble <- ensemble( models[[1]], models[[2]], models[[3]],
                        models[[4]], models[[5]], models[[6]] ,
                        data=d.predict )
mu.avg <- apply( m.ensemble$link , 2 , mean )
mu.pi <- apply( m.ensemble$link , 2 , PI, 0.97 )
predictions.pi <- apply(m.ensemble$sim, 2, PI, .97)

plot(height ~ age, d4p, xlim=c(-2,3.5), ylim=c(0, 200), col=col.alpha('slateblue', 0.5))
lines(age.seq, mu.avg, col='red')
shade(mu.pi, age.seq, col = col.alpha("blue",0.15), border=1)
shade(predictions.pi, age.seq)
mtext(paste('ensemble'))

## 6H4
#m1
calc_dev <- function(data, mu, sigma) {
  -2*sum(dnorm(data$height, mu, sigma,log = TRUE))
}
d4p <- d2
cf <- coef(models[[1]])
mu <- cf[['a']] + d4p$age*cf[['b1']]
tdev1 <- calc_dev(d4p, mu, cf[['sigma']])

cf <- coef(models[[2]])
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']]
tdev2 <- calc_dev(d4p, mu, cf[['sigma']])

cf <- coef(models[[3]])
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']] + d4p$age^3*cf[['b3']]
tdev3 <- calc_dev(d4p, mu, cf[['sigma']])

cf <- coef(models[[4]])
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']] + d4p$age^3*cf[['b3']] + d4p$age^4*cf[['b4']]
tdev4 <- calc_dev(d4p, mu, cf[['sigma']])

cf <- coef(models[[5]])
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']] + d4p$age^3*cf[['b3']] + d4p$age^4*cf[['b4']] + d4p$age^5*cf[['b5']]
tdev5 <- calc_dev(d4p, mu, cf[['sigma']])

cf <- coef(models[[6]])
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']] + d4p$age^3*cf[['b3']] + d4p$age^4*cf[['b4']] + d4p$age^5*cf[['b5']] + d4p$age^6*cf[['b6']]
tdev6 <- calc_dev(d4p, mu, cf[['sigma']])

# 6H5
tdev <- c(tdev1, tdev2, tdev3, tdev4, tdev5, tdev6)
d.tdev <- tdev-min(tdev)
names(d.tdev) <- paste0('m',1:6)
sort(d.tdev)
models.cmp

### 6H6
b.sd <- 5
formulas66 <- alist(
  height ~ dnorm(mu , sigma),
  mu <- a + b1*age + b2*age^2 + b3*age^3 + b4*age^4 + b5*age^5 + b6*age^6,
  a ~ dnorm(0, a.sd),
  b1 ~ dnorm(0, b.sd),
  b2 ~ dnorm(0, b.sd),
  b3 ~ dnorm(0, b.sd),
  b4 ~ dnorm(0, b.sd),
  b5 ~ dnorm(0, b.sd),
  b6 ~ dnorm(0, b.sd),
  sigma ~ dunif(0, sigma.rb)
)
m66 <- map(formulas66, data=d1)
m66
models[[4]]
coeftab(models[[4]],m66)

plot_model('6 regularised', m66, d.predict, d1)

d4p <- d2
cf <- coef(m66)
mu <- cf[['a']] + d4p$age*cf[['b1']] + d4p$age^2*cf[['b2']] + d4p$age^3*cf[['b3']] + d4p$age^4*cf[['b4']] + d4p$age^5*cf[['b5']] + d4p$age^6*cf[['b6']]
tdev66 <- calc_dev(d4p, mu, cf[['sigma']])
tdev66
tdev4
