library(rethinking)
library(rstan)
library(dplyr)
#source("plot_bindings.R")


### Easy ####

# 8E1
# 1(-) Method can work with continuous.
# 2(-) Could be any, not only Gaussian.
# 3(+) should be symmetric/ Quote from the book: 
#  The Metropolis algorithm works whenever the probability of proposing a jump to B from A is equal to the probability of proposing A from B, when the proposal distribution is symmetric

# 8E2
# It explores posterior distribution more efficiently, so there are fewer rejections compared to base Metropolis algorithm.
# Gibbs sampling leverages information about the analytic form of likelihood and conjugate priors, so it can merge proposal and reject/accept steps.
# Quote: The improvement arises from adaptive proposals in which the distribution of proposed parameter values adjusts itself intelligently, depending upon the parameter values at the moment.
# AFAIC: Gibbs sampling on each step makes a move along single coordinate of multidimensional parameters and use conditional distribution of these parameters on all others that is known because of conjugate priors and special form of likelihood. But this is only intuitive understanding.

# 8E3
# Hamiltonian Monte Carlo method couldn't handle discrete parameters. Method is based on the physical metaphor of the moving particle. This particle moves in the space of parameters with speed proportional to the gradient of likelihood, so the method should be able to differentiate the space. Discrete parameters aren't differentiable.

# 8E4
# N_effective aims to estimate the number of 'ideal' samples. Ideal samples are entirely uncorrelated. Due to way MCMC works each next sample is actually correlated with the previous one to some extent. n_eff estimate how many independent samples we should get to collect same information as presented by the current sequence of posterior samples.

# 8E5
# Rhat should be close to 1. It reflects the fact that inner variance and outer variance between chains is roughly the same, so we could expect that inference is not broken (does't depend on the chain).

### Medium ####
# 8M1  ####
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd.trim <- dd[ , c("log_gdp","rugged","cont_africa") ]

m8m1.ch <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0,2)
  ) ,
  data=dd.trim )
precis(m8m1.ch)
pairs(m8m1.ch)

m8m1.un <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ) ,
  data=dd.trim )
precis(m8m1.un)
pairs(m8m1.un)

m8m1.exp <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dexp(1)
  ) ,
  data=dd.trim )
precis(m8m1.exp)
pairs(m8m1.exp)

# plot sigma densities
sigma.ch <- extract.samples( m8m1.ch)$sigma
sigma.un <- extract.samples( m8m1.un)$sigma
sigma.ex <- extract.samples( m8m1.exp)$sigma

par(mfrow=c(1,1))
dens(sigma.ch, xlim=c(0.7, 1.2), col='red')
dens(sigma.un, add=T, col='blue')
dens(sigma.ex, add=T)

# There is no visual difference between resulting models, at least difference is indistinguishable from several runs of the same model.

# 8M2  ####
scale.seq <- c(10, 5, 2, 1, 0.1, 0.01, 0.001)
sigma.seq <- seq(0,10, by=0.1)
## cauchy ##
## priors
i=1
n=length(scale.seq)
plot(sigma.seq, dcauchy(sigma.seq, 0, scale.seq[1]), type='l', xlim=c(0,10), ylim=c(0,0.1), col=i)
for(scale in scale.seq){
  lines(sigma.seq, dcauchy(sigma.seq, 0, scale), type='l', col=i, lty=i, lwd=2)
  i=i+1
}
legend(8, 0.09, legend=scale.seq, col=1:n,  lty=1:n, lwd=2)

m10 <- map2stan( 
    alist(
      log_gdp ~ dnorm( mu , sigma ) ,
      mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
      a ~ dnorm(0,100),
      bR ~ dnorm(0,10),
      bA ~ dnorm(0,10),
      bAR ~ dnorm(0,10),
      sigma ~ dcauchy(0, 10)
    ) ,
    data=dd.trim )

m5 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0, 5)
  ) ,
  data=dd.trim )

m1 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0, 5)
  ) ,
  data=dd.trim )

mdot1 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0, .1)
  ) ,
  data=dd.trim )

mdot01 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0, .01)
  ) ,
  data=dd.trim )

sigma.samples <- list(
  '10'=extract.samples( m10)$sigma,
  '5'=extract.samples( m5)$sigma,
  '1'=extract.samples( m1)$sigma,
  '0.1'=extract.samples( mdot1)$sigma,
  '0.01'=extract.samples( mdot01)$sigma
)
i=1
n=length(sigma.samples)
dens(sigma.samples[['10']], xlim=c(0.7, 1.2), ylim=c(0,10), col=i)
for(scale.name in names(sigma.samples)){
  dens(sigma.samples[[scale.name]], col=i, add=T, lwd=2)
  i=i+1
}
legend(1.1, 8, legend=names(sigma.samples), col=1:n,  lwd=2)


scale.seq <- c(10, 5, 2, 1, 0.1, 0.01, 0.001)
sigma.seq <- seq(0,10, by=0.1)
## exp ##
## priors
i=1
n=length(scale.seq)
plot(sigma.seq, dexp(sigma.seq, scale.seq[1]), type='l', xlim=c(0,10), ylim=c(0,2), col=i)
for(scale in scale.seq){
  lines(sigma.seq, dexp(sigma.seq,  scale), type='l', col=i, lty=i, lwd=2)
  i=i+1
}
legend(8, 1.5, legend=scale.seq, col=1:n,  lty=1:n, lwd=2)

m10 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dexp(10)
  ) ,
  data=dd.trim )

m1 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dexp(1)
  ) ,
  data=dd.trim )

mdot1 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dexp(.1)
  ) ,
  data=dd.trim )

mdot01 <- map2stan( 
  alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dexp(.01)
  ) ,
  data=dd.trim )

sigma.samples <- list(
  '10'=extract.samples( m10)$sigma,
  '1'=extract.samples( m1)$sigma,
  '0.1'=extract.samples( mdot1)$sigma,
  '0.01'=extract.samples( mdot01)$sigma
)
i=1
n=length(sigma.samples)
dens(sigma.samples[['10']], xlim=c(0.7, 1.2), ylim=c(0,10), col=i)
for(scale.name in names(sigma.samples)){
  dens(sigma.samples[[scale.name]], col=i, add=T, lwd=2)
  i=i+1
}
legend(1.1, 8, legend=names(sigma.samples), col=1:n,  lwd=2)

# The answer is that choice of Cauchy priors or exponentially scaled priors does not influence the result of inference.
# The only difference that I've encountered happens for exp(10) - mean is tiny shifted towards zero in this case. Still, it can be just the effect of the random nature of sampling.
# I assume that 170 observations are enough to overcome priors.

# 8M3 ####
formula8m3 <- alist(
    log_gdp ~ dnorm( mu , sigma ) ,
    mu <- a + bR*rugged + bA*cont_africa + bAR*rugged*cont_africa ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    bAR ~ dnorm(0,10),
    sigma ~ dcauchy(0,2)
  ) 

m8m3.w1 <- map2stan(formula8m3, data=dd.trim, iter=1001, warmup=1)
precis(m8m3.w1)#awfull results, n_eff=1

m8m3.w10 <- map2stan(formula8m3, data=dd.trim, iter=1010, warmup=10)
precis(m8m3.w10)#not so awfull results, n_eff=~100..200, troubles with estimates bA & sigma
plot(m8m3.w10)

m8m3.w100 <- map2stan(formula8m3, data=dd.trim, iter=1100, warmup=100)
precis(m8m3.w100)#enough of warmup
plot(m8m3.w100)

m8m3.w500 <- map2stan(formula8m3, data=dd.trim, iter=1500, warmup=500)
precis(m8m3.w500)#"wasted" warmup
plot(m8m3.w500)

m8m3.w1k <- map2stan(formula8m3, data=dd.trim, iter=2000, warmup=1000)
precis(m8m3.w1k)#"wasted" warmup
plot(m8m3.w1k)

## Hard ####
# 8H1 ####
mp <- map2stan(
  alist(
    a ~ dnorm(0, 1),
    b ~ dcauchy(0, 1)
  ),
  data=list(y=1),
  start=list(a=0, b=0),
  iter=1e4, warmup=100, WAIC=FALSE)
# There is no model for y in the formula, so data doesn't influence inference.
# Model will be sampling from the specified prioir distributions.
# We expect to see normal distribution for a and cauchy for b.
plot(mp)
precis(mp)
samples <- extract.samples(mp)
hist(samples$a)#perfect Gausssian
hist(samples$b)#hard to examine, cause has sevral extreme values (as expected from cauchy thick tails)

# 8H2 ####
data(WaffleDivorce)
d <- WaffleDivorce

normalise <- function(x){
  (x-mean(x))/sd(x)
}

d$MedianAgeMarriage_s <- normalise(d$MedianAgeMarriage)
d$Marriage_s <- normalise(d$Marriage)

# fit model
m5.1 <- map2stan(
  alist(
    Divorce ~ dnorm( mu , sigma ) ,
    mu <- a + bA * MedianAgeMarriage_s ,
    a ~ dnorm( 10 , 10 ) ,
    bA ~ dnorm( 0 , 1 ) ,
    sigma ~ dcauchy(0 , 2)
  ), 
  data = select(d, Divorce, MedianAgeMarriage_s) )
precis(m5.1)
plot(m5.1)
pairs(m5.1)

m5.2 <- map2stan(
  alist(
    Divorce ~ dnorm( mu , sigma ) ,
    mu <- a + bR * Marriage_s ,
    a ~ dnorm( 10 , 10 ) ,
    bR ~ dnorm( 0 , 1 ) ,
    sigma ~ dcauchy( 0 , 2 )
  ) , data = select(d, Divorce, Marriage_s) )
precis(m5.2)
plot(m5.2)
pairs(m5.2)

m5.3 <- map2stan( 
             alist(
               Divorce ~ dnorm( mu , sigma ) ,
               mu <- a + bR*Marriage_s + bA*MedianAgeMarriage_s ,
               a ~ dnorm( 10 , 10 ) ,
               bR ~ dnorm( 0 , 1 ) ,
               bA ~ dnorm( 0 , 1 ) ,
               sigma ~ dcauchy( 0 , 2)
             ) ,
             data = select(d, Divorce, MedianAgeMarriage_s, Marriage_s)
             )
precis(m5.3)
plot(m5.3)
pairs(m5.3)

coeftab(m5.1, m5.2, m5.3)
compare(m5.1, m5.2, m5.3)
# Results a similar to inferred by map - m5.1 model is the best one compared by WAIC.
# bR is centered around zero in m5.3. That points to the fact that after MedianAgeMarriage_s is known, marriage rate brings a tiny portion of information in predicting divorce rate.
# Also coefficients bA and bR are highly correlated, as we can see from their plot.

# 8H3 ####
N <- 100 # number of individuals
height <- rnorm(N,10,2) # simulate total height of each individual
leg_prop <- runif(N,0.4,0.5) # leg as proportion of the height
leg_left <- leg_prop*height + # simulate left leg as proportion + error
  rnorm( N , 0 , 0.02 )
leg_right <- leg_prop*height + # simulate right leg as proportion + error
  rnorm( N , 0 , 0.02 )
# combine into data frame
d <- data.frame(height,leg_left,leg_right)

m5.8s <- map2stan(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dcauchy( 0 , 1 )
  ) ,
  data=d, chains=4,
  start=list(a=10,bl=0,br=0,sigma=1) )
plot(m5.8s)
precis(m5.8s)
pairs(m5.8s)

m5.8s2 <- map2stan( 
                    alist(
                      height ~ dnorm( mu , sigma ) ,
                      mu <- a + bl*leg_left + br*leg_right ,
                      a ~ dnorm( 10 , 100 ) ,
                      bl ~ dnorm( 2 , 10 ) ,
                      br ~ dnorm( 2 , 10 ) & T[0,] ,
                      sigma ~ dcauchy( 0 , 1 )
                    ) ,
                    data=d, chains=4,
                    start=list(a=10,bl=0,br=0,sigma=1) )
plot(m5.8s2)
precis(m5.8s2)
pairs(m5.8s2)

coeftab(m5.8s, m5.8s2)
compare(m5.8s, m5.8s2)
# second model has worse convergence and smaller n_eff for bl, br coefficients

# but sum of bl+br stay the same
par(mfrow=c(2,2))
ss5.8s <- extract.samples(m5.8s)
b1 <- ss5.8s$bl+ss5.8s$br
dens(ss5.8s$bl, col='red', ylim=c(0, .25), xlim=c(-12, 12))
abline(v=0, lty=2)
mtext('b_left model=m5.8s')
dens(ss5.8s$br, col='blue', ylim=c(0, .25), xlim=c(-12, 12))
abline(v=0, lty=2)
mtext('b_right model=m5.8s')

ss5.8s2 <- extract.samples(m5.8s2)
b2 <- ss5.8s2$bl+ss5.8s2$br
dens(ss5.8s2$bl, col='red', ylim=c(0, .25), xlim=c(-12, 12))
abline(v=0, lty=2)
mtext('b_left model=m5.8s2(truncated prior)')
dens(ss5.8s2$br, col='blue', ylim=c(0, .25), xlim=c(-12, 12))
abline(v=0, lty=2)
mtext('b_right model=m5.8s2(truncated prior)')

# plot sum of bl+br for both models
par(mfrow=c(1,1))
dens(b1, col='red', ylim=c(0,6.5), xlim=c(1.5,2.5))
dens(b2,add=T, col='blue')

# 8H4 ####
# using models from 8H3
WAIC(m5.8s)  # pWAIC = 3.2937
WAIC(m5.8s2) # pWAIC = 2.819
# Effective number of parameters for the second model is smaller. 
# Intuitively it is smaller because we restricted "freedom" of the 'br' coefficient. This parameter couldn't be negative for the second model, while the probability of having big values is still very small as for the first model. Thus overall freedom of the model declined.
# More formally, pWAIC is defined as sum of variance of the points likelihood, thus the second model has smaller variance of data likelihood(==> it's 'more restricted')

#same story with effective numbet of parameters for DIC
DIC(m5.8s) # pD=3.9
DIC(m5.8s2) # pD=3.4

# 8H5 ####
population <- c(10, 60, 20, 100, 30)
n_islands <- length(population)
n_trials <- 1e+5
positions <- rep(0, n_trials)
curr_pos <- 1
for(i in 1:n_trials){
  positions[i] <- curr_pos
  next_pos <- ifelse(runif(1)<0.5, -1, 1) + curr_pos
  if (next_pos <= 0){
    next_pos <- n_islands
  } else if (next_pos > n_islands){
    next_pos <- 1
  }
  p_ratio <- population[next_pos] / population[curr_pos]
  if ( runif(1) < p_ratio ){
    curr_pos <- next_pos
  }
}
hist(positions)
table(positions)/n_trials
population/sum(population)

# 8H6 ####
# data: number of 'water' tosses out of N total tosses
n_water = 9
n_tosses = 10
# priors: 
#    - flat prios -> just a constant  P_prior(p) = some_c
get_prior_prob <- function(p){
  #1 #flat
  dcauchy(p, 0.5, .1)
}

n_trials <- 1e+5 #number of steps of MCMC
param.samples <- rep(0, n_trials)
curr_param <- 0.01
curr_prob <- dbinom(n_water, n_tosses, prob = curr_param) * get_prior_prob(curr_param)
param_delta <- 0.01 #step for the parameter - equivalent of a jump to the next island
for(i in 1:n_trials){
  param.samples[i] <- curr_param
  next_param <- ifelse(runif(1)<0.5, -1, 1)*param_delta + curr_param #symetric move
  if (next_param <= 0){
    next_param <- 0
  } else if (next_param >= 1){
    next_param <- 1
  }
  next_prob <- dbinom(n_water, n_tosses, prob = next_param) * get_prior_prob(next_param)
  p_ratio <- next_prob / curr_prob
  if ( runif(1) < p_ratio ){
    curr_param <- next_param
    curr_prob <- next_prob
  }
}
par(mfrow=c(1,1))
hist(param.samples)
abline(v=n_water/n_tosses, col='red', lty=2, xlim=c(0,1))
plot(param.samples[1:1000], type='l')
plot(param.samples, type='l', ylim=c(0,1))
abline(h=n_water/n_tosses, col='red', lty=2)


## grid approximation
#get_prior_prob <- function(p){
#  #1 #flat
#  dcauchy(p, 0.5, .1)
#}
#param_delta <- 0.01
p.grid <- seq(0,1,by=param_delta)
prior <- get_prior_prob(p.grid)
likelihood <- sapply(p.grid, function(x) dbinom(n_water, n_tosses, prob=x))
posterior <- prior * likelihood
posterior <- posterior / sum(posterior)
par(mfrow=c(3,1))
plot(p.grid, prior, type='l', lty=2)
plot(p.grid, likelihood, type='l', lty=3, col='red')
abline(v=n_water/n_tosses)
plot(p.grid, posterior, type='l', col='blue')
abline(v=p.grid[which.max(posterior)])
print(p.grid[which.max(posterior)])
