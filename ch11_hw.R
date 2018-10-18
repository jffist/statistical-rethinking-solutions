library(rethinking)
library(rstan)
library(dplyr)
library(ggplot2)
#source("plot_bindings.R")

## Easy ####
# 11E1 #####
# There is a natural ordering of levels within ordered category variable. For example rating of the movie 4 is greater(better) than rating 3. But the difference between levels is not equal and usually unknown/subjective. Increasing a rating from 1 to 2, in general, is very different from moving it from 4 to 5. 
# On the contrary, levels of the unordered categorical variable are not comparable. Examples of unordered categorical variable include gender (Male/Female) or geographical region (EMEA, NAM).

# 11E2 #####
# Ordered logistic regression employs 'cumulative logit link' function. It's a logit function calculated on the cumulative probability function of levels. For each level, this function returns a sum of probabilities of all levels less than or equal to a given one (P(y<=k))

# 11E3 #####
# It will underestimate true value of "lambda" for Poisson regression or parameter "p" for Binomial regression, as it would try to explain large number of zeros with small values of lambda.
# It can be illustrated using an example from the chapter:
### generate data
prob_drink <- 0.2 # 20% of days
rate_work <- 1 # average 1 manuscript per day
# sample one year of production
N <- 365
# simulate days monks drink
drink <- rbinom( N , 1 , prob_drink )
# simulate manuscripts completed
y <- (1-drink)*rpois( N , rate_work )

# fit zero-inflated Poisson 
m11e3 <- map( 
  alist(
    y ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm(0,1),
    al ~ dnorm(0,10)
  ) ,
  data=list(y=y) )
precis(m11e3)
# al==-0.05
exp(-0.05) #0.95 - close to true value 1

# fit ordinary Poisson
m11e3.p <- map( 
  alist(
    y ~ dpois(lambda ),
    log(lambda) <- al,
    al ~ dnorm(0,10)
  ) ,
  data=list(y=y) )
precis(m11e3.p)
# al=-0.24
exp(-0.24) #0.78 - understimated(lower) then real 1 and infered by zero-inflated model

compare(m11e3, m11e3.p)# m11e3 is slightly better 

# 11E4 #####
# One of possible artificial examples of under-dispersed variable is a Binomial varibale with all observations smaller than a threshold substituted with the threshold itself. 
# Also Binomial variable would be underdispersed if it exposes autocorrelation. By construction, Binomial is a sum of independent Bernoulli variables. If these variables are correlated, then the final sum will be underdispersed and no longer Binomial.

n <- 10
p <- 0.8
# plain binomial
y <- rbinom(1e+4, n, p)
print(sprintf("theoretical mean: %f, sample mean %f", n*p, mean(y)))
print(sprintf("theoretical sd: %f, sample sd %f", sqrt(n*p*(1-p)), sd(y)))
simplehist(y)

# thresholded binomial
yc <- ifelse(y<6, 6, y)
print(sprintf("theoretical mean: %f, sample mean %f", n*p, mean(yc)))
print(sprintf("theoretical sd: %f, sample sd %f", sqrt(n*p*(1-p)), sd(yc)))

# Medium ####
# 11M1 #####
r_cnt <- c(12, 36, 7, 41)
p_cum <- cumsum(r_cnt) /  sum(r_cnt)
log_cum_odds <- log(p_cum/(1-p_cum)) 
print(log_cum_odds) #last value is Inf and shoould be ignored

# 11M2 ####
plot(1:4, p_cum)
lines(1:4, p_cum)
prev <- 0
prev2 <- 0
for(i in 1:4){
  lines(c(i,i),c(0,p_cum[i]), lwd=4)
  lines(c(i+0.03,i+0.03), c(prev, p_cum[i]), lwd=4, col='blue')
  #if(i>1){
  #  lines(c(i-1+0.03, i+0.03), c(prev2,  prev))
  #}
  prev2 <- prev
  prev <- p_cum[i]
}

# 11M3 #####
# just change Poisson likelihood to Binomial 
# pz - probability of zero
# n - number of trials
# p - probability of success in trial
# Pr(0|pz,n,p) = pz + (1-pz)*(1-p)^n
# Pr(k|pz,n,p) = (1-pz)*P_binom(k,n,p) = (1-pz)*C(n,k)*p^k*(1-p)^n-k

# Hard ####
# 11H1 #####
data(Hurricanes)
d <- Hurricanes
str(d)
ggplot(d, aes(as.factor(female), femininity))+geom_boxplot()
ggplot(d, aes(as.factor(female), deaths))+geom_violin()
ggplot(d, aes(femininity, deaths))+geom_point() + geom_smooth()
# Visual exploration shows that there is a bunch of "outliers" among female cases, but as with UCBadmit data there could be a confounding variable

m11h1.int <- map(alist(
  deaths ~ dpois(lambda),
  log(lambda) ~ a ,
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h1.int)

m11h1.fem <- map(alist(
  deaths ~ dpois(lambda),
  log(lambda) ~ a + b_fem*femininity,
  b_fem ~ dnorm(0, 10),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h1.fem) #b_fem = 0.07
exp(0.07) # 1.07
coef <- extract.samples(m11h1.base)
dens(coef$b_fem) #coefficient is strongly positive

(cmp <- compare(m11h1.int, m11h1.fem))
plot(cmp)
plot(coeftab(m11h1.int, m11h1.fem))
# According to WAIC comparison model with femininity is better and take all score, but the dispersion of the difference and SE of WAIC itself is huge. 
# The difference is less than standard error of the difference. 

postcheck(m11h1.fem, window = 100)
abline(h=40, col='red')
abline(h=10, col='blue')
abline(h=mean(d$deaths), lty=2)
# Visual exploration shows that model significantly underestimates deaths for hurricanes with a number of deaths greater than 40 and overestimates counts for cases with a number of death less than 10.
# Model is good at predicting number of deaths for hurricanes with target variable close to the average across the sample.

# Because model contains only a single variable, we can draw counter factual plot to illustrate dependency
d.predict <- data.frame(femininity=seq(1,11,0.1))
lambda.sample <- link(m11h1.fem, d.predict)
lambda.avg <- apply(lambda.sample, 2, mean )
lambda.pi <- apply(lambda.sample, 2, PI )

# predict actual counts
count.sample <- sim(m11h1.fem, data = d.predict)
count.avg <- apply(count.sample, 2, mean )
count.pi <- apply(count.sample, 2, PI )

#plot
plot(d$femininity, d$deaths, xlim=c(0,12), col='blue', pch=16)
lines(d.predict$femininity, lambda.avg)
shade(lambda.pi, d.predict$femininity) #gives very narrow shade

lines(d.predict$femininity, count.avg, col='red')
shade(count.pi, d.predict$femininity) #shade of counts predictions

# Summary: Intuitively there is some hidden variable that better explains deaths.  
# Visually relation induced by the model looks suspicious for me because it looks like being caused by several outliers.

# 11H2 #####
m11h1.fem.gamma <- map(alist(
  deaths ~ dgampois(mu, theta),
  log(mu) ~ a + b_fem*femininity,
  b_fem ~ dnorm(0, 10),
  a ~ dnorm(0, 10),
  theta ~ dexp(1)
),
data=d
)
precis(m11h1.fem.gamma)
# now b_fem is centered around zero, so looks like no correlation with deaths

postcheck(m11h1.fem.gamma, window = 100)
# Overdispersed model, as expected, has a wide interval for predictions of counts. 
# These intervals cover actual counts for almost all cases, except 8 cases with the number of deaths greater than 50. 
# Also all predictions are nearly the same and are close to the average of the sample.

## yet anothe plot of femininity vs deaths
lambda.sample <- link(m11h1.fem.gamma, d.predict)
lambda.avg <- apply(lambda.sample, 2, mean )
lambda.pi <- apply(lambda.sample, 2, PI )

# predict actual counts
count.sample <- sim(m11h1.fem, data = d.predict)
count.avg <- apply(count.sample, 2, mean )
count.pi <- apply(count.sample, 2, PI )

#plot
plot(d$femininity, d$deaths, xlim=c(0,12), col='blue', pch=16)
lines(d.predict$femininity, lambda.avg)
shade(lambda.pi, d.predict$femininity) #gives very narrow shade

lines(d.predict$femininity, count.avg, col='red')
shade(count.pi, d.predict$femininity) #shade of counts predictions

# 11H3 #####
normalise <- function(x){
  (x-mean(x))/sd(x)
}

d$damage_norm_c <- normalise(d$damage_norm)
d$femininity_c <- normalise(d$femininity)
d$min_pressure_c <- normalise(d$min_pressure)



m11h3 <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*damage_norm_c + b_mp*min_pressure_c,
  c(b_fem,b_dam,b_mp) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3)
postcheck(m11h3, window = 100)

m11h3.fxd <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*damage_norm_c + b_mp*min_pressure_c +
                b_fem_dam*femininity_c*damage_norm_c,
  c(b_fem,b_dam,b_mp,b_fem_dam) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3.fxd)
postcheck(m11h3.fxd, window = 100)
pairs(m11h3.fxd)


m11h3.fxmp <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*damage_norm_c + b_mp*min_pressure_c +
    b_fem_mp*femininity_c*min_pressure_c,
  c(b_fem,b_dam,b_mp,b_fem_mp) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3.fxmp)
postcheck(m11h3.fxmp, window = 100)
pairs(m11h3.fxmp)

m11h3.fxd.fxmp <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*damage_norm_c + b_mp*min_pressure_c +
    b_fem_mp*femininity_c*min_pressure_c + b_fem_dam*femininity_c*damage_norm_c,
  c(b_fem,b_dam,b_mp,b_fem_dam,b_fem_mp) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3.fxd.fxmp)
postcheck(m11h3.fxd.fxmp, window = 100)
pairs(m11h3.fxd.fxmp)


m11h3.no.fem <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a +  b_dam*damage_norm_c + b_mp*min_pressure_c,
  c(b_dam,b_mp) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)

m11h3.no.mp <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*damage_norm_c,
  c(b_fem,b_dam) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)

m11h3.fxd.only <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem_dam*femininity_c*damage_norm_c,
  c(b_fem_dam) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)

cmp <- compare(m11h3, m11h3.fxd, m11h3.fxmp, m11h3.fxd.fxmp, m11h3.no.fem, m11h3.no.mp, m11h3.fxd.only)
cmp
plot(cmp)

# Best model according to the WAIC comparison changes across runs. So it is hard to distinguish which one of models m11h3.fxd or m11h3 is better.
# Thus adding interaction of damage vs feminity doesn't improve inference a lot. 
# Removing femininity variable from the model increases WAIC, but it is still within single deviance of the WAIC difference between models.
# It's interesting that removing min_pressure_c variable produces results comparable to the best model, but with huge deviance of difference.

# adhoc visualisation of poisson models built on male/female subsets
d %>% ggplot(aes(damage_norm_c, deaths, group=female, color=as.factor(female))) + geom_point(size=3) + geom_smooth(method = 'glm', method.args=list(family=poisson)) 

# visualisation 
predict_lambda_counts <- function(model, data){
  lambda.sample <- link(model, data = data)
  lambda.avg <- apply(lambda.sample, 2, mean )
  lambda.pi <- apply(lambda.sample, 2, PI )
  
  count.sample <- sim(model, data = data)
  count.avg <- apply(count.sample, 2, mean )
  count.pi <- apply(count.sample, 2, PI )
  
  list(
    l_avg=lambda.avg,
    l_pi=lambda.pi,
    cnt_avg=count.avg,
    cnt_pi=count.pi
  )
}

plot_lambda_cnt <- function(x, pred, color_name) {
  lines(x, pred$l_avg, col=color_name)
  shade(pred$l_pi, x) 
  
  lines(x, pred$cnt_avg, col=color_name, lty=2)
  shade(pred$cnt_pi, x) #shade of counts predictions
}

## plot for avg min_pressure_c
model <- m11h3.no.mp#m11h3
damage_seq <- seq(-1,5.5,0.1)
d.predict.male <- data.frame(
  femininity_c=-1.3,
  damage_norm_c=damage_seq,
  min_pressure_c=0
)
d.predict.female <- data.frame(
  femininity_c=1,
  damage_norm_c=damage_seq,
  min_pressure_c=0
)

p.male <- predict_lambda_counts(model, d.predict.male)
p.female <- predict_lambda_counts(model, d.predict.female)


idx.male <- d$female!=1
idx.female <- d$female==1
plot(d$damage_norm_c[idx.male], d$deaths[idx.male], xlim=range(damage_seq), pch=16, col='blue', ylim=range(d$deaths))
points(d$damage_norm_c[idx.female], d$deaths[idx.female], pch=16, col='red')
plot_lambda_cnt(d.predict.male$damage_norm_c, p.male, 'blue')
plot_lambda_cnt(d.predict.male$damage_norm_c, p.female, 'red')

# >>triptych plot for min_pressure_c ####

model <- m11h3 #m11h3.no.mp#m11h3
damage_seq <- seq(-1,5.5,0.1)
qvar_name <- "min_pressure_c"
par(mfrow=c(1,3))
left_q <- 0
for(right_q in c(0.33, 0.66, 1)){
  # calculate qunatiles range
  qvar_vec <- d[[qvar_name]]
  qq <- quantile(qvar_vec, probs=c(left_q, right_q))
  lval = qq[1]
  rval = qq[2]
  # filter data subset
  if(right_q!=1){
    d.raw <- d[(qvar_vec >= lval) & (qvar_vec < rval),]
  } else {
    d.raw <- d[(qvar_vec >= lval) & (qvar_vec <= rval),]
  }
  # calc avg of the variable
  qvar_avg = mean(d.raw[[qvar_name]])
  # create data for prediction
  d.predict.male <- data.frame(
    femininity_c=-1.3,
    damage_norm_c=damage_seq,
    min_pressure_c=qvar_avg
  )
  d.predict.female <- data.frame(
    femininity_c=0.66,
    damage_norm_c=damage_seq,
    min_pressure_c=qvar_avg
  )
  ## predict
  p.male <- predict_lambda_counts(model, d.predict.male)
  p.female <- predict_lambda_counts(model, d.predict.female)
  ## plot
  idx.male <- d.raw$female!=1
  idx.female <- d.raw$female==1
  plot(d.raw$damage_norm_c[idx.male], d.raw$deaths[idx.male], xlim=range(damage_seq), pch=16, col='blue', ylim=range(d$deaths))
  points(d.raw$damage_norm_c[idx.female], d.raw$deaths[idx.female], pch=16, col='red')
  plot_lambda_cnt(d.predict.male$damage_norm_c, p.male, 'blue')
  plot_lambda_cnt(d.predict.female$damage_norm_c, p.female, 'red')
  mtext(sprintf("%s=%5.4f (%3.2f, %3.2f)", qvar_name, qvar_avg, left_q, right_q))
  ## end
  left_q <- right_q
}
# These plots show that for low min_pressure values difference between male and female lines is bigger.  
# Still, there are more female dots with higher death values.
# For high values of min_pressure difference is tiny.

# 11H4 ####
d$log_damage_norm_c <- normalise(log(d$damage_norm))
m11h3.logd <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a + b_fem*femininity_c + b_dam*log_damage_norm_c + b_mp*min_pressure_c,
  c(b_fem,b_dam,b_mp) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3.logd)
par(mfrow=c(1,1))
postcheck(m11h3.logd, window = 100)

m11h3.logd.fxd <- map(alist(
  deaths ~ dpois(mu),
  log(mu) ~ a  + b_dam*log_damage_norm_c  + b_fem*femininity_c 
    b_fem_dam*femininity_c*log_damage_norm_c,
  c(b_fem,b_dam,b_fem_dam) ~ dnorm(0, 2),
  a ~ dnorm(0, 10)
),
data=d
)
precis(m11h3.logd.fxd)
postcheck(m11h3.logd.fxd, window = 100)
# interesting observation - b_fem becomes zero, but b_fem_dam - not

(cmp <- compare(m11h3, m11h3.logd, m11h3.logd.fxd, m11h3.fxd))
plot(cmp)
plot(coeftab(m11h3, m11h3.logd, m11h3.logd.fxd, m11h3.fxd))

# m11h3.logd.fxd - is the winner
model <- m11h3.logd.fxd
log_damage_seq <- seq(-3.1, 2, 0.1)
d.predict.male <- data.frame(
  femininity_c=-1.3,
  log_damage_norm_c=log_damage_seq,
  min_pressure_c=0
)
d.predict.female <- data.frame(
  femininity_c=1,
  log_damage_norm_c=log_damage_seq,
  min_pressure_c=0
)

p.male <- predict_lambda_counts(model, d.predict.male)
p.female <- predict_lambda_counts(model, d.predict.female)


idx.male <- d$female!=1
idx.female <- d$female==1
plot(d$log_damage_norm_c[idx.male], d$deaths[idx.male], xlim=range(log_damage_seq), pch=16, col='blue', ylim=range(d$deaths))
points(d$damage_norm_c[idx.female], d$deaths[idx.female], pch=16, col='red')
plot_lambda_cnt(log_damage_seq, p.male, 'blue')
plot_lambda_cnt(log_damage_seq, p.female, 'red')

# With log scale of damage there is no need for min_pressure variable.
# Model becomes more accurate in predictions.
# According to coefficients it's only interaction of femininicity and log(damage_norm) that really matters.


d$pred_m11h3 <- apply(sim(m11h3), 2, mean)
d$pred_m11h4 <- apply(sim(m11h3.logd.fxd), 2, mean)

calc_r2 <- function(d, pred_col){
  avg=mean(d$deaths)
  1-sum((d$deaths - d[pred_col])**2)/sum((d$deaths - avg)**2)
}
calc_r2(d,'pred_m11h3')
calc_r2(d,'pred_m11h4')
# R-squared of log-damage model is much better, but this is R2 on the train sample, that could mislead because of the overfitting

# 11H5 ####
data("Trolley")
d <- Trolley
str(d)

# quick look through the visualisation
d %>% filter(contact==1) %>% ggplot( aes(x=as.factor(response), group=as.factor(male), fill=as.factor(male))) + 
  geom_bar(aes(y=..prop..), position = "dodge") + 
  ggtitle("Probabiliy of response per gender for questions that involves contact")
# Plot illustrates that women tend to have lower proportion of answers rated as 6 or 7 

# We are interested in checking how gender influences decision.
# We expect that women have more responses that will qualify story as immoral when contact is included.
# I suggest to take the best model from the chapter and add gender variable to it for testing the hypothesis.

# best model from the chapter
m11h5.base <- map( 
  alist(
    response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ) ,
    phi <- bA*action + bI*intention + bC*contact + bAI*action*intention + bCI*contact*intention,
    c(bA,bI,bC,bAI,bCI) ~ dnorm(0,10),
    c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
  ) ,
  data=d ,
  start=list(a1=-1.9,a2=-1.2,a3=-0.7,a4=0.2,a5=0.9,a6=1.8) )
precis(m11h5.base)

# model that contains gender and gender vs contact interaction
m11h5.f <- map( 
  alist(
    response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ) ,
    phi <- bA*action + bI*intention + bC*contact + 
      bAI*action*intention + bCI*contact*intention + bF*(1-male) +
      bFC*(1-male)*contact,
    c(bA,bI,bC,bAI,bCI,bFC,bF) ~ dnorm(0,10),
    c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
  ) ,
  data=d ,
  start=list(a1=-1.9,a2=-1.2,a3=-0.7,a4=0.2,a5=0.9,a6=1.8) 
)
precis(m11h5.f)

# model that contains only gender vs contact interaction
m11h5.fc <- map( 
  alist(
    response ~ dordlogit( phi , c(a1,a2,a3,a4,a5,a6) ) ,
    phi <- bA*action + bI*intention + bC*contact + 
           bAI*action*intention + bCI*contact*intention +
           bFC*(1-male)*contact,
    c(bA,bI,bC,bAI,bCI,bFC) ~ dnorm(0,10),
    c(a1,a2,a3,a4,a5,a6) ~ dnorm(0,10)
  ) ,
  data=d ,
  start=list(a1=-1.9,a2=-1.2,a3=-0.7,a4=0.2,a5=0.9,a6=1.8) 
)
precis(m11h5.fc)

(cmp <- compare(m11h5.base, m11h5.fc, m11h5.f))
plot(cmp)
plot(coeftab(m11h5.base, m11h5.fc, m11h5.f))
# According to WAIC comparison model with gender and gender-contact is significantly better and takes all weights.

# Sign of coefficient near gender-contact interaction term changes across models, but it's correlated with the coefficient for gender variable, so let's examine changes in posterior distribution of predictions instead of coefficients

## >>>11h5 visualisation of distrubution shift due to gender variable #####

post <- extract.samples( m11h5.f )
str(post)

plot(1, 1, type="n", xlab="gender==female", ylab="probability", xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )
abline(h=c(0,1), lty=2, col='blue')

kI <- 0 # value of intention
kA <- 0 # value for action 
kC <- 1 # value for contact
kF <- c(0,1)# values of gender==female flag to calculate over
for ( s in 1:100 ) {
  p <- post[s,]
  ak <- as.numeric(p[1:6])
  phi <- p$bA*kA + p$bI*kI + p$bC*kC + p$bAI*kA*kI + p$bCI*kC*kI + p$bF*kF + p$bFC*kF*kC
  pk <- pordlogit( 1:6 , a=ak , phi=phi )
  for ( i in 1:6 )
    lines( kF , pk[,i] , col=col.alpha(rangi2,0.1) )
}
mtext( concat( "action=",kA,", contact=",kC,", intention=",kI ) )

# calculate phi for all posterior samples for male
kF <- 0
post$phi_m <- with(post, bA*kA + bI*kI + bC*kC + bAI*kA*kI + bCI*kC*kI + bF*kF + bFC*kF*kC)
# calculate phi for all posterior samples for female
kF <- 1
post$phi_f <- with(post, bA*kA + bI*kI + bC*kC + bAI*kA*kI + bCI*kC*kI + bF*kF + bFC*kF*kC)

n <- nrow(post)
a.mx = data.matrix(select(post, a1:a6))
# calculate mean cumulative and per response probabilities from all posterior sample
# male
p_m <- sapply(1:n, function(idx) pordlogit(1:6 , a=a.mx[idx,] , phi=post$phi_m[idx] )  )
p_m <- apply(p_m, 1, mean)
p_m <- c(p_m,1)
cp_m <- p_m
p_m <- p_m-c(0,p_m)[1:7]
# female
p_f <- sapply(1:n, function(idx) pordlogit(1:6 , a=a.mx[idx,] , phi=post$phi_f[idx] )  )
p_f <- apply(p_f, 1, mean)
p_f <- c(p_f,1)
cp_f <- p_f
p_f <- p_f-c(0,p_f)[1:7]
# add values to the plot
text(1, cp_f-0.04, sprintf("%3.2f",p_f))
text(0, cp_m-0.04, sprintf("%3.2f",p_m))

# From the plot we see that females tend to have a larger proportion of responses equal to 1, 2 or 3 and smaller proportion of responses equal to 5 or 6. 
# So it looks that data support the hypothesis.

# 11H6 ####
data("Fish")
d <- Fish
str(d)

m11h6.base <- map(alist(
  fish_caught ~ dzipois(p, lambda),
  logit(p) <-  ap,
  log(lambda) <- log(hours) + al,
  ap ~ dnorm(0,10),
  al ~ dnorm(0,10)
), data=d)
precis(m11h6.base)
postcheck(m11h6.base, window=250 )
logistic(-0.75) # probability of not fishing is 0.32
exp(-0.14) # avg number of caught fish per hour when fishing is 0.87 vs. mean(d$fish_caught/d$hours)=1.1

m11h6 <- map(alist(
  fish_caught ~ dzipois(p, lambda),
  logit(p) <-  ap + bp_c*camper + bp_p*persons + bp_nchd*child,
  log(lambda) <- log(hours) + al + bl_lb*livebait + bl_c*camper + bl_p*persons + bl_nchd*child,
  ap ~ dnorm(0,10),
  al ~ dnorm(0,10),
  c(bp_c, bp_p, bp_nchd) ~ dnorm(0,2),
  c(bl_lb, bl_c, bl_p, bl_nchd) ~ dnorm(0,2)
), data=d)
precis(m11h6)
postcheck(m11h6, window=250 )
pairs(m11h6)


m11h6.no.child <- map(alist(
  fish_caught ~ dzipois(p, lambda),
  logit(p) <-  ap + bp_c*camper + bp_p*persons,
  log(lambda) <- log(hours) + al + bl_lb*livebait + bl_c*camper + bl_p*persons,
  ap ~ dnorm(0,10),
  al ~ dnorm(0,10),
  c(bp_c, bp_p, bp_nchd) ~ dnorm(0,2),
  c(bl_lb, bl_c, bl_p, bl_nchd) ~ dnorm(0,2)
), data=d)
precis(m11h6.no.child)

(cmp <- compare(m11h6.base, m11h6, m11h6.no.child))
plot(cmp)
# adding child to the model has almost no influence on WAIC

# compare to stan model
m11h6.s <- map2stan(alist(
  fish_caught ~ dzipois(p, lambda),
  logit(p) <-  ap + bp_c*camper + bp_p*persons + bp_nchd*child,
  log(lambda) <- log(hours) + al + bl_lb*livebait + bl_c*camper + bl_p*persons + bl_nchd*child,
  ap ~ dnorm(0,10),
  al ~ dnorm(0,10),
  c(bp_c, bp_p, bp_nchd) ~ dnorm(0,2),
  c(bl_lb, bl_c, bl_p, bl_nchd) ~ dnorm(0,2)
), data=d)
precis(m11h6.s)
#
plot(coeftab(m11h6, m11h6.s))
# stan model gives the same results, so can proceed with map

# 
# m11h6.x <- map(alist(
#   fish_caught ~ dzipois(p, lambda),
#   logit(p) <-  ap + bp_c*camper,
#   log(lambda) <- log(hours) + al + bl_lb*livebait  + bl_p*persons,
#   ap ~ dnorm(0,10),
#   al ~ dnorm(0,10),
#   c(bp_c, bp_p, bp_nchd) ~ dnorm(0,2),
#   c(bl_lb, bl_c, bl_p, bl_nchd) ~ dnorm(0,2)
# ), data=d)
# precis(m11h6.x)
# (cmp <- compare(m11h6.base, m11h6, m11h6.x))
# plot(cmp)

# predict
model <- m11h6.no.child
d.predict <- data.frame(
  hours=1,
  livebait=c(0,1,0,1),
  camper=c(0,0,1,1),
  persons=c(1,1,1,1),
  child=c(0,0,0,0)
)
lambda.sample <- link(model, data = d.predict)
lambda.avg <- apply(lambda.sample$lambda, 2, mean )
lambda.pi <- apply(lambda.sample$lambda, 2, PI )

p.avg <- apply(lambda.sample$p, 2, mean )
p.pi <- apply(lambda.sample$p, 2, PI )

count.sample <- sim(model, data = d.predict)
count.avg <- apply(count.sample, 2, mean )
count.pi <- apply(count.sample, 2, PI )

d.predict$lambda <- lambda.avg
d.predict$p <- p.avg
d.predict$cnt <- count.avg

# Group with camper has more chances to start fishing, but lower number of caught fishes per hour, that's counterfactual. I suspect that there is some correlation inside the model. 
# As expected, a group that uses livebait has larger expected number of caught fishes per hour.
