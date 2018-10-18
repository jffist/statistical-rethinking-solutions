library(rethinking)
library(rstan)
library(dplyr)
library(ggplot2)
library(GGally)
library(MASS)
#source("plot_bindings.R")


### Easy ####

# 10E1
p <- 0.35
p/(1-p)

# 10E2
lo <- 3.2
lo/(1+lo)

# 10E3
# log ods will change in exp(1.7) times
exp(1.7)

# 10E4
# Poisson distribution models number of events per some unit of time/spatial region.
# Measurement of the outcome variable can be provided on the different scale for each observation (daily vs weekly).
# Offset is used to bring all observations on the same scale.
# * Consider an example from the chapter, where a number of manuscripts produced by the monastery is measured on the daily or weekly basis. The offset parameter is used to convert all measurements to the daily basis.  
# * Number of animals in the area is another possible example where offset is helpful. Square of the area can be treated as an offset in this case.
# * While predicting number of conversions per ad campaign, number of impressions can be treated as an offset. I think it is a Poisson regression and not Binomial, because conversion rate is extremely small usually and number of impressions is huge (1 vs 1e+6)

## Medium ####
# 10M1
# Binomial likelihood in the aggregated form contains C(n,m) multiplier. This multiplier is converted to an additional constant at the log-scale of the likelihood. For non-aggregated format, each event is modeled independently as a number of heads in the single drop of the coin.
# likelihood in the aggregated format: C(n,m)*p^m*(1-p)^n-m
# likelihood in the non-aggregated format: p^m*(1-p)^n-m

# 10M2
# It implies that the change of the predictor by 1 unit increases the lambda parameter of the Poisson distribution in exp(1.7)=5.4739 times. Thus an average number of events per interval increases in 5.47 times.

# 10M3
# 10M3. Explain why the logit link is appropriate for a binomial generalized linear model.
# Binomial likelihood is parametrised by parameter p - probability of an event
# We are interested in modelling it with linear combination of the predictors.
# But p shoudl fall in the [0, 1] range
# So instead of p we model f(p) as linear relation a + b*x, f is selected in a such way that
# p = f^-1(a+b*x) is in [0,1] scale
# f=logit=log(p/1-p), f^-1=logistic=exp(a+b*x)/(1+exp(a+b*x))
# logit function ensures required contraint.

# 10M4
# As in the previous case(10M3) we are interested in modelling lambda - parameter of the Poisson distribution with linear model. 
# Lambda should be in [0, +inf) range. An exponential function is an appropriate inverse of the link that maps (-inf, +inf)  range of linear combination of parameters into (0,+inf). Then the link function should be logarithm as the inverse of exp.

# 10M5
# Using logit link implies that a lambda parameter of the Poisson likelihood always falls in the (0,1) range.
# There is nothing wrong technically, not sure if it makes sense to set such a constraint for some real process.
# Here are simulations of such process:
prop.table(table(rpois(1e+5, .1)))
# One of possible cases to model is a process with the high number of "zero" events and with the absence of exponential growth of the predicted rate. 

# 10M6
# The binomial distribution has maximum entropy when each trial must result in one of two possible events and the expected value is constant. Poisson distribution is a special case of the Binomial one. Just because there is no other max entropy distribution for constraints of the Binomial, I assume that Poisson distribution adds some new constraint.
# From the form of the Poisson distribution it is known that its variance is equal to the expected value and both are constant. So I assume that this is the additional constraint for the Poisson, but I didn't try to prove this.

## Hard ####
# 10H1 ####
data("chimpanzees")
d <- chimpanzees
d$recipient <- NULL

m10.1q2 <- map(alist(
  pulled_left ~ dbinom( 1 , p ) ,
  logit(p) <- a[actor] + (bp + bpC*condition)*prosoc_left ,
  a[actor] ~ dnorm(0,10),
  bp ~ dnorm(0,10),
  bpC ~ dnorm(0,10)
) ,
data=d)
pairs(m10.1q2)

m10.1stan <- map2stan(m10.1q2, chains=2 , iter=2500 , warmup=500 )
pairs(m10.1stan)

plot(coeftab(m10.1q2, m10.1stan))
# The main difference lies in the posterior distribution of the intercept for the second monkey(a[2]). It is highly skewed according to the Stan estimations but is symmetric in the quadratic estimations. All other coefficients have similar estimated distribution.
# Intercept a[2] has large valuer and long tail because the second actor always pulled left lever, so inference forces intercept to be as large as possible. E.g. logistic(MAP(a[2]))=logistic(11)=0.99998


# 10H2 ####
# model with single intercept for all actors
m10.1stan.si <- map2stan(alist(
  pulled_left ~ dbinom( 1 , p ) ,
  logit(p) <- a + (bp + bpC*condition)*prosoc_left ,
  a ~ dnorm(0,10),
  bp ~ dnorm(0,10),
  bpC ~ dnorm(0,10)), 
  data=d, 
  iter=2500, 
  warmup=500 
)
pairs(m10.1stan.si)
precis(m10.1stan.si, depth=2)

# model with intercept per actor, but with no coefficient for condition variable
m10.1stan.nc <- map2stan(alist(
  pulled_left ~ dbinom( 1 , p ) ,
  logit(p) <- a[actor] + bp*prosoc_left ,
  a[actor] ~ dnorm(0,10),
  bp ~ dnorm(0,10)
  ), 
  data=d, 
  iter=2500, 
  warmup=500 
)
pairs(m10.1stan.nc)
precis(m10.1stan.nc, depth=2)

# model with intercept per actor but no coefficient for condition and for prosoc option
m10.1stan.io <- map2stan(alist(
  pulled_left ~ dbinom( 1 , p ) ,
  logit(p) <- a[actor],
  a[actor] ~ dnorm(0,10)
), 
data=d, 
iter=2500, 
warmup=500 
)
pairs(m10.1stan.io)
precis(m10.1stan.io, depth=2)

cmp <- compare(m10.1stan.io, m10.1stan.nc,  m10.1stan, m10.1stan.si)
# m10.1stan.si -  a model that has a single intercept for all actors is the worst one (has largest WAIC) and gains no score in per model weights distribution.
# m10.1stan.io - a model that consists of single intercept per actor is comparable by WAIC, but still has a difference with the best model that is greater than deviance of the difference. Thus it also gains no score in model weights distribution.
# There is a tiny difference in WAIC between models with and without `condition` variable. Model without this variable looks better in terms of estimated WAIC and gains 70% of the overall score. 
# 
plot(cmp)

plot(coeftab(m10.1stan.io, m10.1stan.nc,  m10.1stan, m10.1stan.si))
# All models estimate MAP of `bp` coefficient as non-zero

# 10H3 ####
library(MASS)
data(eagles)
d <- eagles
str(d)

# 10H3.a
d$p_is_large <- ifelse(d$P=='L', 1, 0)
d$p_is_adult <- ifelse(d$A=='A', 1, 0)
d$v_is_large <- ifelse(d$V=='L', 1, 0)

d$prob <- d$y/d$n
  
m10.3q2 <- map(alist(
    y ~ dbinom(n, p),
    logit(p) <- a + bp*p_is_large + bv*v_is_large + ba*p_is_adult,
    a ~ dnorm(0, 10),
    c(bp, bv, ba) ~ dnorm(0, 5)
  ), data=d)
precis(m10.3q2)
pairs(m10.3q2)

m10.3stan <- map2stan(m10.3q2)
precis(m10.3stan)
pairs(m10.3stan)

# While models' WAIC is comparable there is a difference in parameters estimation.
# Parameters `bp` and `bv` are negatively correlated and have long tails - both coefficients are not symmetric related to the mean.
# One of the possible causes of this phenomena is that for pairs (P=L,V=S, A={A,I}) probability of pirate to win is always equal to 1.
# It forces coefficient towards large values to shift logistic score closer to 1 (as in 10H1)

# 10H3.b
par(mfrow=c(1,1))
postcheck(m10.3stan, prob=0.89) #it plots predicted and actual probablity with their HPDI intervals
# let's do it manually
d$lbl <- with(d, paste0(P,A,V))
post_dens_probs <- function(model, d){
  prob.sample <- link(model)
  prob.mean <- apply(prob.sample, 2, mean)
  prob.pi <- apply(prob.sample, 2, PI, 0.89)
  
  dr <- d
  dr$p_pred <- prob.mean
  dr$p_pred_low <- prob.pi[1,]
  dr$p_pred_high <- prob.pi[2,]
  
  par(mfrow=c(4,2))
  for(i in 1:8){
    dens(prob.sample[,i]) #, xlim=c(0,1)
    abline(v=d$prob[i], col='red')
    abline(v=prob.mean[i], col='blue')
    abline(v=prob.pi[,i], col='blue', lty=2)
    mtext(paste('Probability for case#',i,dr$lbl[i]))
  }
  dr
}

# plot counts
post_dens_counts <- function(model, d){
  cnt.sample <- sim(model)
  cnt.mean <- apply(cnt.sample, 2, mean)
  cnt.pi <- apply(cnt.sample, 2, PI, 0.89)
  
  dr <- d
  dr$y_pred <- cnt.mean
  dr$y_pred_low <- cnt.pi[1,]
  dr$y_pred_high <- cnt.pi[2,]
  dr

  par(mfrow=c(4,2))
  for(i in 1:8){
    dens(cnt.sample[,i]) 
    abline(v=d$y[i], col='red')
    abline(v=cnt.mean[i], col='blue')
    abline(v=cnt.pi[,i], col='blue', lty=2)
    mtext(paste('Count for case#',i,dr$lbl[i]))
  }
  dr
}

dr <- post_dens_probs(m10.3stan, d)
dr <- post_dens_counts(m10.3stan, dr)
par(mfrow=c(1,1))
postcheck(m10.3stan, prob=0.89)

# case 8 has the worst predictions

par(mfrow=c(1,1))
plot(dr$y, col='blue', pch=16,  xlab="case", xaxt="n" )
points(dr$y_pred)
for(i in 1:8){
  lines( c(i, i), c(dr$y_pred_low[i], dr$y_pred_high[i]) )
}
axis(1, at=1:8, labels=dr$lbl)

# From the visual exploration we can recognise a pattern, that interaction of adult pirate with small victim always gives an advantage to the pirate (for counts plot but not for probabilities)
# Let's check this hypothesis with the model

# 10H3.c
m10.3bpa <- map2stan(alist(
  y ~ dbinom(n, p),
  logit(p) <- a + bp*p_is_large + bv*v_is_large + ba*p_is_adult + bpa*p_is_large*p_is_adult,
  a ~ dnorm(0, 10),
  c(bp, bv, ba, bpa) ~ dnorm(0, 5)
), data=d)
precis(m10.3bpa)
pairs(m10.3bpa)
compare(m10.3stan, m10.3bpa)

m10.3bpa <- map2stan(alist(
  y ~ dbinom(n, p),
  logit(p) <- a + bp*p_is_large + bv*v_is_large + ba*p_is_adult + bpa*p_is_large*p_is_adult,
  a ~ dnorm(0, 10),
  c(bp, bv, ba, bpa) ~ dnorm(0, 5)
), data=d)
precis(m10.3bpa)
pairs(m10.3bpa)
compare(m10.3stan, m10.3bpa)
postcheck(m10.3bpa, prob=0.89)



m10.3bva <- map2stan(alist(
  y ~ dbinom(n, p),
  logit(p) <- a + bp*p_is_large + bv*v_is_large + ba*p_is_adult + bva*(1-v_is_large)*p_is_adult,
  a ~ dnorm(0, 10),
  c(bp, bv, ba, bva) ~ dnorm(0, 5)
), data=d)
precis(m10.3bva)
pairs(m10.3bva)
compare(m10.3stan, m10.3bpa, m10.3bva)


m10.3bvpa <- map2stan(alist(
  y ~ dbinom(n, p),
  logit(p) <- a + bp*p_is_large + bv*v_is_large + ba*p_is_adult + bva*(1-v_is_large)*p_is_adult + bpa*p_is_large*p_is_adult,
  a ~ dnorm(0, 10),
  c(bp, bv, ba, bpa, bva) ~ dnorm(0, 5)
), data=d)
precis(m10.3bvpa)
pairs(m10.3bvpa)
cmp <- compare(m10.3stan, m10.3bpa, m10.3bva, m10.3bvpa)
plot(cmp)

par(mfrow=c(1,1))
postcheck(m10.3bvpa, prob=0.89)

dr <- post_dens_probs(m10.3bpa, d)
dr <- post_dens_counts(m10.3bpa, dr)

#bpa model - all predictors with large-x-age interaction looks like most appropriate(smallest WAIC and reasonable set of parameters)

# let's check parameters values
pairs(m10.3bpa)
precis(m10.3bpa)

# a - coefficient when all predictors are zero (SIS case)
logistic(c(-0.78, -2.45, 0.77))
# bp - change in odds if pirate is large
exp(6.6)
# bv - change in odds if victim is large
exp(-5.31)
# I think there is a issue with interpreting this coefficitents, as bp and bv are correlated, also there is an interaction
plot(coeftab(m10.3stan, m10.3bpa))
# from this comparison we can see, that adding interaction term increases bp and ba coefficients as we have negative bpa coefficent. Also as ba and bpa are negatively correlated they influence each other.

pd <- position_dodge(0.1) # move points .05 to the left and right
# check influence of the pirate is large varibale on the probability (join by to other variables)
dr$lablel_av <- paste0('v.size=',dr$V,' ','p.age=',dr$A)
ggplot(dr, aes(x=lablel_av, y=p_pred, color=P, group=P)) + 
  geom_point(size=3, position=pd) +
  geom_line( position=pd) + 
  geom_errorbar(aes(ymin=p_pred_low, ymax=p_pred_high), width=.1, position=pd)

# check influence of the pirate is adult varibale on the probability (join by to other variables)
dr$lablel_pv <- paste0('p.size=',dr$P,' ','v.size=',dr$V)
ggplot(dr, aes(x=lablel_pv, y=p_pred, color=A, group=A)) + 
  geom_point(size=3, position=pd) + 
  geom_line( position=pd) +
  geom_errorbar(aes(ymin=p_pred_low, ymax=p_pred_high), width=.1, position=pd)

# check influence of the victim is large varibale on the probability (join by to other variables)
dr$lablel_pa <- paste0('p.size=',dr$P,' ','p.age=',dr$A)
ggplot(dr, aes(x=lablel_pa, y=p_pred, color=V, group=V)) + 
  geom_point(size=3, position=pd) + 
  geom_line( position=pd) + 
  geom_errorbar(aes(ymin=p_pred_low, ymax=p_pred_high), width=.1, position=pd)


#####
dr.base <- post_dens_probs(m10.3stan, d)
dr.base$model = 'base'
dr.bpa <- post_dens_probs(m10.3bpa, d)
dr.bpa$model = 'bpa'

dr2 <- bind_rows(dr.base, dr.bpa)

dr <- dr2
dr$lablel_av <- paste0('v.size=',dr$V,' ','p.age=',dr$A)
ggplot(dr, aes(x=lablel_av, y=p_pred, color=paste(P,model), group=paste(P,model))) + 
  geom_point(size=3, position=pd) +
  geom_line(position=pd) + 
  geom_errorbar(aes(ymin=p_pred_low, ymax=p_pred_high), width=.1, position=pd) +
  geom_point(aes(x=lablel_av, y=prob), size=3, color='black')
# This image shows that bpa models better reflect fact, that move from pirate.age=A to pirate.age=I for small victims has bigger impact on probability due to pirate size-age interaction
# Still, the simplest interpreation is that being large and adult pirate gives additional bonus in the fights for food :)

## 10H4 ####
data(salamanders)
d <- salamanders
nrow(d)
str(d)
pairs(d)
ggpairs(d)

par(mfrow=c(1,1))
# 10H4.a
summary(d$PCTCOVER)
sd(d$PCTCOVER)
d$PCTCOVER_cntr <- (d$PCTCOVER - mean(d$PCTCOVER))/sd(d$PCTCOVER)
m10h4 <- map(alist(
      SALAMAN ~ dpois(lambda),
      log(lambda) <- a + bc*PCTCOVER,
      a ~ dnorm(0,10),
      bc ~ dnorm(0, 5)
  ), 
  data=d
)
pairs(m10h4)
precis(m10h4)
exp(0.03)
# change of the coverage on 1% changes avg number of salamanders in 1.03 times

# from model with centered pctcover
# The change in one standart deviance if cover(~35%) increases avg number of salamander in 3 times

m10h4.s <- map2stan(m10h4)
pairs(m10h4.s)
precis(m10h4.s)
# stan and quadratic results are comparable in shape and values, let's proceed with quadratic as it's faster

#postcheck(m10h4)
d.predict <- data.frame(PCTCOVER=1:100)
lambda.sample <- link(m10h4, n = 1e+4, data=d.predict)
lambda.avg <- apply(lambda.sample, 2, mean)
lambda.PI <- apply(lambda.sample, 2, PI, .89)

counts.sample <- sim(m10h4, data=d.predict)
counts.PI <- apply(counts.sample, 2, PI, .89)

plot(d$PCTCOVER, d$SALAMAN, col=col.alpha('blue', 0.5), pch=16, xlim=c(0,100))
lines(d.predict$PCTCOVER, lambda.avg, lwd=2)
shade(lambda.PI, d.predict$PCTCOVER, col=col.alpha('blue', 0.5))
shade(counts.PI, d.predict$PCTCOVER,  col=col.alpha('green', 0.5))

# 10H4b
d$FORESTAGE_cnt <- (d$FORESTAGE - mean(d$FORESTAGE))/sd(d$FORESTAGE)
m10h4.f <- map(alist(
  SALAMAN ~ dpois(lambda),
  log(lambda) <- a  + bf*FORESTAGE_cnt,
  a ~ dnorm(0,10),
  bf ~ dnorm(0, 2)
), 
data=d
)
precis(m10h4.f)
compare(m10h4, m10h4.f)
pairs(m10h4.f)

d.predict <- data.frame(FORESTAGE_cnt=seq(-1,3,by=0.1))
lambda.sample <- link(m10h4.f, n = 1e+4, data=d.predict)
lambda.avg <- apply(lambda.sample, 2, mean)
lambda.PI <- apply(lambda.sample, 2, PI, .89)

counts.sample <- sim(m10h4.f, data=d.predict)
counts.PI <- apply(counts.sample, 2, PI, .89)

plot(d$FORESTAGE_cnt, d$SALAMAN, col=col.alpha('blue', 0.5), pch=16, xlim=c(-1,3))
lines(d.predict$FORESTAGE_cnt, lambda.avg, lwd=2)
shade(lambda.PI, d.predict$FORESTAGE_cnt, col=col.alpha('blue', 0.5))
shade(counts.PI, d.predict$FORESTAGE_cnt,  col=col.alpha('green', 0.5))

m10h4.cf <- map(alist(
  SALAMAN ~ dpois(lambda),
  log(lambda) <- a  + bf*FORESTAGE_cnt + bc*PCTCOVER,
  a ~ dnorm(0,10),
  c(bc, bf) ~ dnorm(0, 2)
), 
data=d
)
precis(m10h4.cf)
(cmp <- compare(m10h4, m10h4.f, m10h4.cf))
plot(cmp)
pairs(m10h4.cf)

# Adding forestage variable to the model doesn't help anyhow (One can expect that age of the trees doesn't influence life of salamanders, as they live on the ground)

