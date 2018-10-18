rm(list=ls())
library(rethinking)
library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
source("plot_bindings.R")

## Easy ####
# 12E1 ##
# First one, as it is more 'regularising'(narrow because of lower sigma), thus I expect it forces estimates to be shifted to the mean more than the second option. 
# It results in stronger shrinkage.

# 12E2 ##
# y_i ~ Binomial(1,p_i)
# logit(p_i) = a[i] + beta*x_i
# a[i] ~ Normal(a_l2, sigma_l2)
# beta ~ Normal(0,1)
# a_l2 ~ Normal(0, 1)
# sigma_l2 ~ HalfCauchy(0,1)

# 12E3 ## (same as 12E2)
# ...
# a[i] ~ Normal(a_l2, sigma_l2)
# a_l2 ~ Normal(0, 1)
# sigma_l2 ~ HalfCauchy(0,1)
# ...

# 12E4 ##
# y_i ~ Poisson(lambda_i)
# log(lambda_i) = a[i] + beta*x_i
# a[i] ~ Normal(a_l2, sigma_l2)
# beta ~ Normal(0,1)
# a_l2 ~ Normal(0, 1)
# sigma_l2 ~ HalfCauchy(0,1)

# 12E5 ##
# y_i ~ Poisson(lambda_i)
# log(lambda_i) = a_l1 + a1[i] + a2[i] + beta*x_i
# a_l1 ~ Normal(0, 10)
# beta ~ Normal(0,1)
# a1[i] ~ Normal(0, sigma1_l2)
# a2[i] ~ Normal(0, sigma2_l2)
# sigma1_l2 ~ HalfCauchy(0,1)
# sigma2_l2 ~ HalfCauchy(0,1)

## Medium ####
# 12M1, 12M2 ####
data("reedfrogs")
d <- reedfrogs
str(d)

d$has_pred<- as.integer(d$pred=='pred')
d$is_big<- as.integer(d$size=='big')
d$tank <- 1:nrow(d)

m12m1.base <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] ,
    a_tank[tank] ~ dnorm( a , sigma ),
    a ~ dnorm(0, 1) ,
    sigma ~ dcauchy(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.base, depth=2)

m12m1.size <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] + beta_s*is_big,
    a_tank[tank] ~ dnorm( a , sigma ),
    a ~ dnorm(0, 10) ,
    sigma ~ dcauchy(0,1),
    beta_s ~ dnorm(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.size, depth=2)
# First unexpected observation - adding size variable makes inference less stable(n_eff for some tanks significantly reduced from 4000 to ~500)
cmp <- compare(m12m1.base, m12m1.size)
plot(cmp)
plot(coeftab(m12m1.base, m12m1.size))
# accoring to WAIC models are identical

m12m1.pred <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] + beta_p*has_pred,
    a_tank[tank] ~ dnorm( a , sigma ),
    a ~ dnorm(0,1) ,
    sigma ~ dcauchy(0,1),
    beta_p ~ dnorm(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.pred, depth=2)

cmp <- compare(m12m1.base, m12m1.size, m12m1.pred)
plot(cmp)
plot(coeftab(m12m1.base, m12m1.size, m12m1.pred))
coeftab(m12m1.base, m12m1.size, m12m1.pred)
# Interesting observation that signs for some a_tank coefficient have changed.
# Predator variable explained part of the variance in the observations, thus sigma reduced.

m12m1.sp <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] +  beta_s*is_big + beta_p*has_pred,
    a_tank[tank] ~ dnorm( a , sigma ),
    a ~ dnorm(0,1) ,
    sigma ~ dcauchy(0,1),
    beta_p ~ dnorm(0,1),
    beta_s ~ dnorm(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.sp, depth=2)
(cmp <- compare(m12m1.base, m12m1.size, m12m1.pred, m12m1.sp))
plot(cmp)
plot(coeftab( m12m1.pred, m12m1.sp)) # size gives almost no influence
coeftab(m12m1.base, m12m1.size, m12m1.pred, m12m1.sp)


m12m1.sp.i2 <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] +  beta_s*is_big + beta_p*has_pred + b_sp*is_big*has_pred,
    a_tank[tank] ~ dnorm( a , sigma ),
    a ~ dnorm(0,1) ,
    sigma ~ dcauchy(0,1),
    c(beta_p, beta_s, b_sp) ~ dnorm(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.sp.i2, depth=2)
(cmp <- compare(m12m1.base, m12m1.size, m12m1.pred, m12m1.sp,  m12m1.sp.i2))
plot(cmp)

coeftab(m12m1.pred, m12m1.sp.i2)
plot(coeftab(m12m1.pred, m12m1.base))

# Models m12m1.pred, m12m1.sp and  m12m1.sp.i2 are very similar according to WAIC comparison.
# The most significant improvement in WAIC occurred after adding has_pred variable. Beside decreasing WAIC it reduced sigma estimates (dispersion across tanks).

postcheck(m12m1.base, window=60)
postcheck(m12m1.pred, window=60)


predict_p_counts <- function(model, data, model_name){
  p.sample <- link(model, data = data)
  p.avg <- apply(p.sample, 2, mean )
  p.pi <- apply(p.sample, 2, PI )
  
  count.sample <- sim(model, data = data)
  count.avg <- apply(count.sample, 2, mean )
  count.pi <- apply(count.sample, 2, PI )
  
  pred <- list(
    p_avg=p.avg,
    p_pi=p.pi,
    cnt_avg=count.avg,
    cnt_pi=count.pi
  )
  d.res <- data.frame(data)
  d.res$model <- model_name
  d.res$pred_p <- pred$p_avg
  d.res$pred_p_low <- pred$p_pi[1,]
  d.res$pred_p_high <- pred$p_pi[2,]
  d.res
}
d1 <- predict_p_counts(m12m1.base, d, 'base')
d2 <- predict_p_counts(m12m1.pred, d, '+predator')
d.res <- bind_rows(d1, d2)

pd <- position_dodge(0.2)
d.res %>% ggplot(aes(x=tank,color=model)) + 
  geom_point(aes(y=propsurv, color='raw'), size=2) + 
  geom_point(aes(y=pred_p), size=2, position=pd) + 
  geom_errorbar(aes(ymin=pred_p_low, ymax=pred_p_high), position=pd, width=1) 

# According to visual exploration width of credible intervals is smaller for the model with the predator variable, and prediction is better(again only visually)

# Plots illustrates that for most tanks m12m1.pred has smaller credible interval than the base model
pd <- position_dodge(0.2)
d.res %>% mutate(p_pred_witdth=pred_p_high-pred_p_low) %>% 
  ggplot(aes(x=tank, y=p_pred_witdth, fill=model, color=model, group=model)) + 
  #geom_bar(stat='identity', position = 'dodge') + 
  geom_line() + 
  geom_point(size=2) +
  ggtitle('width of credible interval per model')

# 12M3 ####
# reuse from prev exercises
#m12m1.base <- map2stan( 
#  alist(
#    surv ~ dbinom(density , p),
#    logit(p) <- a_tank[tank] ,
#    a_tank[tank] ~ dnorm( a , sigma ),
#    a ~ dnorm(0, 1) ,
#    sigma ~ dcauchy(0,1)
#  ), 
#  data=d , iter=2000 , chains=2, cores=2 
#)

m12m1.base.cauchy <- map2stan( 
  alist(
    surv ~ dbinom(density , p),
    logit(p) <- a_tank[tank] ,
    a_tank[tank] ~ dcauchy( a , sigma ),
    a ~ dnorm(0, 1) ,
    sigma ~ dcauchy(0,1)
  ), 
  data=d , iter=2000 , chains=2, cores=2 
)
precis(m12m1.base.cauchy, depth=2)
(cmp <- compare(m12m1.base, m12m1.base.cauchy))
#WAIC is almost the same

plot(coeftab(m12m1.base, m12m1.base.cauchy))
# Model with Cauchy sigma has several intercepts per tank with much greater variance than the model with the Gaussian sigma.
# Also Cauchy model has bigger estimate of mean `a` and smaller `sigma`(yet we couldn't compare sigmas as they come from different distributions)

# plot a distributions
s1 <- extract.samples(m12m1.base)
s2 <- extract.samples(m12m1.base.cauchy)

a_tank_normal <- apply(s1$a_tank, 2, mean)
a_tank_cauchy <- apply(s2$a_tank, 2, mean)

plot(a_tank_normal, a_tank_cauchy, xlab="Gaussian prior" , ylab="Cauchy prior", xlim=c(-2,10), ylim=c(-2,10) )
abline(a=0, b=1, lty=2)
# check examples with very big a_tanks for cauchy priors
d[a_tank_cauchy>4, ]
# All such examples have 100% survival rate, and because of this a_tank can go to infinity unless restricted by priors.
# Big values for a_tank become more probable in the posterior sample with Cauchy priors, because this distribution has wider tales that Gaussian.

# 12M4 ####
data("chimpanzees")
d <- chimpanzees
d$block_id <- d$block # name 'block' is reserved by Stan

m12m4.base <- map2stan(
  alist(
    pulled_left ~ dbinom( 1 , p ),
    logit(p) <- a + a_actor[actor] + a_block[block_id] +
      (bp + bpc*condition)*prosoc_left,
    a_actor[actor] ~ dnorm( 0 , sigma_actor ),
    a_block[block_id] ~ dnorm( 0 , sigma_block ),
    c(a,bp,bpc) ~ dnorm(0,10),
    sigma_actor ~ dcauchy(0,1),
    sigma_block ~ dcauchy(0,1)
  ) ,
  data=d, warmup=1000 , iter=2000 , chains=2 , cores=2 )
precis(m12m4.base, depth=2)
#pairs(m12m4.base)
plot(m12m4.base)

m12m4 <- map2stan(
  alist(
    pulled_left ~ dbinom( 1 , p),
    logit(p) <- a_actor[actor] + a_block[block_id] +
      (bp + bpc*condition)*prosoc_left,
    a_actor[actor] ~ dnorm( alpha , sigma_actor ),
    a_block[block_id] ~ dnorm( gamma , sigma_block ),
    c(alpha, gamma, bp,bpc) ~ dnorm(0,10),
    sigma_actor ~ dcauchy(0,1),
    sigma_block ~ dcauchy(0,1)
  ) ,
  data=d, warmup=1000 , iter=2000 , chains=2 , cores=2 )
precis(m12m4, depth=2)
# The first observation is that sampling of the second model takes enormous time compared to the first one  - this is a bad sign :)
# There were 141 divergent transitions after warmup - also looks like a bad sign.
# Rhat for second model ~1.01 for all coefficients - sign of bad convergence.
#pairs(m12m4)
plot(m12m4)

compare(m12m4.base, m12m4)
# WAIC measure is roughly the same

## 12H1 ####
data("bangladesh")
d <- bangladesh

str(d)
#df <- d %>% mutate(district=as.factor(district), urban=as.factor(urban))
#df %>% 
#    group_by(district,urban) %>% 
#    summarise(avg_use=mean(use.contraception)) %>% 
#    ggplot(aes(district, avg_use, color=urban)) + geom_point()
#df %>% ggplot(aes(as.factor(use.contraception), living.children, fill=urban, color=urban)) + geom_jitter()
d$district_id <- as.integer(as.factor(d$district))
d$use_contraception <- d$use.contraception
d$age_centered <- d$age.centered
d <- select(d,-use.contraception, -age.centered, -district)

m12h1.d.fixed <- map2stan( 
  alist(
    use_contraception ~ dbinom(1, p),
    logit(p) <- a_district[district_id],
    a_district[district_id] ~ dnorm(0, 5)
  ), 
  data=d
)
precis(m12h1.d.fixed, depth=2)

m12h1.d.pooled <- map2stan( 
  alist(
    use_contraception ~ dbinom(1, p),
    logit(p) <- a_district[district_id] ,
    a_district[district_id] ~ dnorm(a, sigma),
    a ~ dnorm(0,5) ,
    sigma ~ dcauchy(0,1)
  ), 
  data=d
)
precis(m12h1.d.pooled, depth=2)

compare(m12h1.d.fixed, m12h1.d.pooled)
# Pooled model looks superior according to WAIC comparison

cacl_district_rate_estimate <- function(model, model_name, district_data){
  sf <- extract.samples(model)
  sf$a_district <- logistic(sf$a_district)
  
  d1 <- data.frame(district_data)
  d1$rate <- apply(sf$a_district, 2, mean)
  pi <- apply(sf$a_district, 2, PI)
  d1$rate_low <- pi[1,]
  d1$rate_high <- pi[2,]
  d1$model <- model_name
  d1
}

d.res <- d %>% 
  group_by(district_id) %>% 
  summarise(
    cnt=n(), 
    ttl_use_c=sum(use_contraception),
    rate=ttl_use_c/cnt
  ) %>% 
  as.data.frame() %>%
  mutate(
    d_label = reorder(as.factor(paste0(district_id,'/n=',cnt)), cnt)
  ) %>% 
  arrange(district_id)
d1 <- cacl_district_rate_estimate(m12h1.d.fixed, 'fixed', d.res)
d2 <- cacl_district_rate_estimate(m12h1.d.pooled, 'pooled', d.res)

d.res.models <- bind_rows(d1,d2)
d.res <- arrange(d.res, cnt) %>% mutate(cnt_sorted_id=1:nrow(d.res))

pd <- position_dodge(0.3)
d.res.models %>% 
  ggplot(aes(d_label, rate, color=model, fill=model)) + 
  geom_point(data=d.res, size=3, aes(fill='raw', color='raw')) +
  geom_point(size=2, position = pd) + 
  geom_errorbar(aes(ymin=rate_low, ymax=rate_high), position = pd) + 
  geom_hline(yintercept = mean(d$use_contraception), linetype='longdash')  + 
  geom_hline(yintercept = mean(d.res$rate), linetype='dotdash' ) + 
  theme(axis.text.x = element_text(angle = 90))
  ggtitle('Predicted vs. empirical rate of women who use contraception per district')
# This plot shows per district rate of contraception use. Districts are sorted by the number of observations from the smallest to the largest.
# We can see that  differences between pooled, fixed and empirical estimates are tiny for the right part of the plot, where number of observations per district is ~100.
# But for teh first 10 districts, that have less than 15 observations, pooled estimates are shrinked towards the sample mean.
  
# 12H2 ####
data("Trolley")
d <- Trolley
str(d)
d$person_id <- as.integer(d$id)

m12h2.base <- map2stan( 
    alist(
      response ~ dordlogit(phi, cutpoints),
      phi <- bA*action + bI*intention + bC*contact + bAI*action*intention + bCI*contact*intention,
      c(bA, bI, bC, bAI, bCI) ~ dnorm(0,10),
      cutpoints ~ dnorm(0,10)
    ) ,
    data=d ,
    start=list(cutpoints=c(-1.9, -1.2, -0.7, 0.2, 0.9, 1.8)),
    types=list(cutpoints="ordered") 
)
precis(m12h2.base, depth=2 )
set.seed(123)
(d.pred <- sample_n(d, 10))
sm <- sim(m12h2.base,data=d.pred)
simplehist(sm[,5])

# fixed effect model
# !!! it took ~2hrs to run it
m12h2.id.fe <- map2stan(
  alist(
   response ~ dordlogit(phi, cutpoints),
     phi <- a_person[person_id] + bA*action + bI*intention + bC*contact + bAI*action*intention + bCI*contact*intention,
    c(bA, bI, bC, bAI, bCI) ~ dnorm(0,10),
    cutpoints ~ dnorm(0,10),
    a_person[person_id] ~ dnorm(0, 10)
  ) ,
  data=d ,
  start=list(cutpoints=c(-1.9, -1.2, -0.7, 0.2, 0.9, 1.8)),
  types=list(cutpoints="ordered")
)
precis(m12h2.id.fe, depth=2) #enormous Rhat(>1.3), n_eff<10, we couldn't trust this model
post <- extract.samples(m12h2.id.fe)
str(post)
hist(post$a_person[,1])
compare(m12h2.base, m12h2.id.fe) #still it's much better according to WAIC

# !!! it took me ~45min to run
m12h2.id.ve <- map2stan(
  alist(
    response ~ dordlogit(phi, cutpoints),
    phi <- a_person[person_id] + bA*action + bI*intention + bC*contact + bAI*action*intention + bCI*contact*intention,
    c(bA, bI, bC, bAI, bCI) ~ dnorm(0,10),
    cutpoints ~ dnorm(0,10),
    a_person[person_id] ~ dnorm(0, sigma_person),
    sigma_person ~ dcauchy(0, 2)
  ) ,
  data=d ,
  start=list(cutpoints=c(-1.9, -1.2, -0.7, 0.2, 0.9, 1.8)),
  types=list(cutpoints="ordered")
)
precis(m12h2.id.ve, depth=2)
compare(m12h2.base, m12h2.id.ve) #m12h2.id.ve is significantly better

# I don't know how to compare uncertainty in predictions among seven ordered categories for 331 individuals.
# So let's plot a_person coefficients. I assume if these coefficients are far from zero it would mean that they have a significant influence on the result.

post <- extract.samples(m12h2.id.ve)
str(post)
a_pi <- apply(post$a_person, 2, PI)
a_avg <- apply(post$a_person, 2, mean)
d.res <- data.frame(person_id=1:length(unique(d$person_id)), a=a_avg, a_low=a_pi[1,], a_high=a_pi[2,])
head(d.res)
d.res$plabel <- reorder(levels(d$id), d.res$a)


d.res %>% ggplot(aes(plabel, a)) + 
        geom_point(size=2) + 
        geom_errorbar(aes(ymin=a_low, ymax=a_high)) + 
        geom_hline(yintercept = 0, linetype='dashed') + 
        theme(axis.text.x = element_text(angle = 90))
# This plots illustrate that clustering per person captures a lot of variances.
# It's interesting to check individuals with a_person approaching extreme values -4 or 8 (p==0 or 1 on logit scale)
head(arrange(d.res, a))
d %>% filter(id=='96;550') #a_person ~= -5, almost all responses are equal to 1
head(arrange(d.res, desc(a)))
d %>% filter(id=='97;644') #a_person ~= 6,  all responses are equal to 7
d %>% filter(id=='98;214') #a_person ~= 6,  all responses are equal to 7

# Check visualisation of predictions for a single person and single case later in the file (around line 522)


# 12H3 #####
# !!!samping of this model takes ~2.5hrs
d$story_id <- as.integer(d$story)
m12h2.id.story.ve <- map2stan(
  alist(
    response ~ dordlogit(phi, cutpoints),
    phi <- a_story[story_id] + a_person[person_id] + bA*action + bI*intention + bC*contact + bAI*action*intention + bCI*contact*intention,
    c(bA, bI, bC, bAI, bCI) ~ dnorm(0,10),
    cutpoints ~ dnorm(0, 10),
    a_person[person_id] ~ dnorm(0, sigma_person),
    a_story[story_id] ~ dnorm(0, sigma_story),
    sigma_person ~ dcauchy(0, 2),
    sigma_story ~ dcauchy(0, 2)
  ) ,
  data=d ,
  start=list(cutpoints=c(-1.9, -1.2, -0.7, 0.2, 0.9, 1.8)),
  types=list(cutpoints="ordered")
)
precis(m12h2.id.story.ve, depth=2)
(cmp <- compare(m12h2.id.ve, m12h2.id.story.ve)) 
plot(cmp)
# model with cluster per story looks better

post <- extract.samples(m12h2.id.story.ve)
str(post)
a_pi <- apply(post$a_person, 2, PI)
a_avg <- apply(post$a_person, 2, mean)
d.res2 <- data.frame(d.res) %>% 
  arrange(person_id) %>%
  mutate(
          a=a_avg, 
          a_low=a_pi[1,],
          a_high=a_pi[2,],
          model='person_id+story'
 )
head(d.res2)
d.res$model='person_id'
d.m.res <- bind_rows(d.res, d.res2)

pd <- position_dodge(0.2)
# plot first 50 a_person
d.m.res %>% arrange(person_id) %>% head(100) %>% ggplot(aes(plabel, a, color=model, fill=model)) + 
  geom_point(size=2, position = pd, pch=21, colour='black') + 
  geom_errorbar(aes(ymin=a_low, ymax=a_high), position = pd) + 
  geom_hline(yintercept = 0, linetype='dashed') + 
  theme(axis.text.x = element_text(angle = 90,  vjust=0.5))
# Quite interesting that values for a_person intercepts stayed almost the same.
# Only for extreme cases, they shifted back from the mean, as, I assume,  variance was partially explained by the story parameters.

# let's check story intercepts
s_pi <- apply(post$a_story, 2, PI)
s_avg <- apply(post$a_story, 2, mean)
d.s.res <- data.frame(story_id=1:length(unique(d$story_id)), a=s_avg, a_low=s_pi[1,], a_high=s_pi[2,])
head(d.s.res)
d.s.res$slabel <- reorder(levels(d$story), d.s.res$a)
d.s.res
d.s.res %>% ggplot(aes(slabel, a)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=a_low, ymax=a_high)) + 
  geom_hline(yintercept = 0, linetype='dashed') + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Avg value and 89 percentile interval of varying intercept per story")
# There are stories that tends to shift decision towards 1 (pon) or towards 7 (aqu) significantly.


## Plot probability per response level for single case and unique individual ####
model_name <- 'm12h2.id.story.ve' #m12h2.id.story.ve #m12h2.base 
model <-  m12h2.id.story.ve       #m12h2.id.story.ve #m12h2.base 

pid <- '96;474' #'96;474' #'96;512'
case_name <- 'cfaqu' #cfaqu, #ikpon
(d.pred <- d %>% filter(id==pid, case==case_name)) 
post <- extract.samples(model)

plot(1,1,type='n',xlim=c(1,7), ylim=c(0,1.0), xlab="response", ylab="probability")
n_samples <- 1000
mx <- matrix(0, nrow=n_samples, ncol=7)
for(i in 1:n_samples) {
  ak <- post$cutpoints[i,]
  phi <- 0
  if('a_story' %in% names(post)){
    phi <- phi + post$a_story[i, d.pred$story_id]
  }
  if('a_person' %in% names(post)){
    phi <- phi + post$a_person[i, d.pred$person_id]
  }
  phi <- phi +
         post$bA[i]*d.pred$action + 
         post$bI[i]*d.pred$intention + 
         post$bC[i]*d.pred$contact + 
         post$bAI[i]*d.pred$action*d.pred$intention + 
         post$bCI[i]*d.pred$contact*d.pred$intention
  pk <- pordlogit( 1:6 , a=ak , phi=phi )
  levels_prob <- c(pk,1) - c(0,pk)
  if(i%%5==0){
    lines(levels_prob, col=col.alpha(rangi2,0.1))
  }
  mx[i,] <- levels_prob
}
lvl_avg <- apply(mx, 2, mean)
lvl_pi <- apply(mx, 2, PI)
points(1:7, lvl_avg)
for(i in 1:7){
  lines(c(i,i), lvl_pi[,i])
}
abline(v=d.pred$response+0.01, col='red', lty=2)
mtext(sprintf("Probability per level for id=%s case=%s | model:%s", pid, case_name, model_name))
