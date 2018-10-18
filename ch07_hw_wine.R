rm(list=ls())
library(rethinking)
library(ggplot2)
library(GGally)
library(dplyr)
#source("plot_bindings.R")

# part 1 - I've started exploration by asking questions and digging into the dataset in the following order
# 1.1 Do avg scores of red and white wines differ?
# 1.2 Do avg scores of American and non-American wines differ for each of 2 subgroups (red/white)
# 1.3 Do avg scores of American and non-American judges differ for each of 4 subgroups (red-american, red-nonamerican, white-american, white-nonamerican)

# This is kind of drill down exploration that uses dimension hierarchy color->region->"judge origin"
# On each level of analysis t.test is used to compare groups averages
# This is a purely frequentist approach. Not sure if I need to use some kind of Bonfernoni correction for multiple tests
# Result of such analysis - avg scores for red wines from america and non-american regions are different

# Part2 - attempt to reproduce same results using Bayesian approach
# I've construcred models that include all varibales that influence subgroup choice
# Make prediction for a subgroups
# Compare distrubutions of mean's predicted for a subgroup
# Depending on choosed model(best one or ensemble) P(avg red wines diff<0|model) ~ 0.88..0.97


data("Wines2012")
d <- Wines2012

### PART 1 ####
# some manual intutive explorations
d%>% group_by(flight) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(d,flight=='red')$score, filter(d,flight!='red')$score)
# no difference in avg scores between red and white wines 

d%>% group_by(wine.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(d,wine.amer==0)$score, filter(d,wine.amer==1)$score)
# no diff between amer/not amer wine scores


# red wines - amer/non-amer score diff
d%>% group_by(wine.amer,flight) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
amer.red <- filter(d, wine.amer==1, flight=='red')$score
nonamer.red <- filter(d, wine.amer==0, flight=='red')$score

hist(amer.red, col=col.alpha('red'))
hist(nonamer.red, col=col.alpha('blue'), add=T)
t.test(amer.red, nonamer.red)
# t.test tells us there is a significant difference in scores

# white wines - amer/non-amer
amer.white <- filter(d, wine.amer==1, flight!='red')$score
nonamer.white <- filter(d, wine.amer==0, flight!='red')$score
hist(amer.white, col=col.alpha('red'))
hist(nonamer.white, col=col.alpha('blue'), add=T)
t.test(amer.white, nonamer.white)
# no diff

# scores among amer/non-amer judges
d%>% group_by(judge.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(d, judge.amer==0)$score, filter(d, judge.amer==1)$score)
# no diff between amer/not amer judge scores


red.wines <- filter(d, flight=='red')
amer.red <- filter(red.wines, wine.amer==1)
nonamer.red <- filter(red.wines, wine.amer==0)

nrow(amer.red)
amer.red%>% group_by(judge.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(amer.red, judge.amer==0)$score, filter(amer.red, judge.amer==1)$score)
#while avg scores are different, there is not enough evidence to prove significance

nonamer.red%>% group_by(judge.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(nonamer.red, judge.amer==0)$score, filter(nonamer.red, judge.amer==1)$score)


white.wines <- filter(d, flight!='red')
amer.wht <- filter(white.wines, wine.amer==1)
nonamer.wht <- filter(white.wines, wine.amer==0)

nrow(amer.wht)
amer.wht%>% group_by(judge.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))

nonamer.wht%>% group_by(judge.amer) %>% summarise(mean(wine.amer), mean(score), n(), max(score), min(score))
t.test(filter(nonamer.wht, judge.amer==0)$score, filter(nonamer.wht, judge.amer==1)$score)
# there is significant diff in scores of non-amer white wines by amer/non-amer judges

### PART 2 ####

##### attempt to do same stuff with Bayesian regression ####
# find diff on color and region (red amer vs red non amer)
d$is_red <- as.integer(d$flight=='red')

m.wf <- map(alist(
  score ~ dnorm(mu, sigma),
  mu <- a + bW*wine.amer + bF*is_red,
  a ~ dnorm(14, 20),
  bW ~ dnorm(0, 1),
  bF ~ dnorm(0, 1),
  sigma ~ dunif(0, 100)
), 
data=d)
precis(m.wf)

m.wf_i <- map(alist(
  score ~ dnorm(mu, sigma),
  mu <- a + bW*wine.amer + bF*is_red + bWF*wine.amer*is_red,
  a ~ dnorm(14, 20),
  bW ~ dnorm(0, 1),
  bF ~ dnorm(0, 1),
  bWF ~ dnorm(0, 1),
  sigma ~ dunif(0, 100)
), 
data=d)
precis(m.wf_i)
compare(m.wf, m.wf_i)

# compare red wines amer vs non-amer using interaction model (wf_i)
samples <- extract.samples(m.wf_i)
str(samples)
mu.red.amer <- with(samples, a+bW+bF+bWF)
mu.red.nonamer <- with(samples, a+bF)
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)

diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)

# compare white wines amer vs non-amer
mu.wht.amer <- with(samples, a+bW)
mu.wht.nonamer <- with(samples, a)
dens(mu.wht.amer, col='red', xlim=c(12,16.5))
dens(mu.wht.nonamer, col='blue', add=T)

diff.wht.region <- mu.wht.amer - mu.wht.nonamer
dens(diff.wht.region)
str(diff.wht.region)
mean(diff.wht.region<0)

# compare red wines amer vs non-amer using plain model (wf)
samples <- extract.samples(m.wf)
str(samples)
mu.red.amer <- with(samples, a+bW+bF)
mu.red.nonamer <- with(samples, a+bF)
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)

diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)

# compare red wines amer vs non-amer using link function and interaction model (wf_i)
d.predict <- data.frame(wine.amer=c(1,0), is_red=c(1,1))
mu <- link(m.wf_i, data=d.predict, n = 1e+4)
str(mu)
mu.red.amer <- mu[,1]
mu.red.nonamer <- mu[,2]
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)
diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)

# compare red wines amer vs non-amer using ensemble model (wf_i+wf)
d.predict <- data.frame(wine.amer=c(1,0), is_red=c(1,1))
m.ens <- ensemble(m.wf, m.wf_i, data = d.predict, n = 1e+4)
mu <- m.ens$link #link(m.wf_i, data=d.predict, n = 1e+4)
str(mu)
mu.red.amer <- mu[,1]
mu.red.nonamer <- mu[,2]
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)
diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)


# add judges to the model
m.jwf_iwf <- map(alist(
  score ~ dnorm(mu, sigma),
  mu <- a + bW*wine.amer + bF*is_red + bWF*wine.amer*is_red + bJ*judge.amer,
  a ~ dnorm(14, 20),
  bW ~ dnorm(0, 1),
  bF ~ dnorm(0, 1),
  bWF ~ dnorm(0, 1),
  bJ ~ dnorm(0, 1),
  sigma ~ dunif(0, 100)
), 
data=d)
precis(m.jwf_iwf)
compare(m.wf, m.wf_i, m.jwf_iwf)

# terms were commented after obsevation that regularised priors set their MAP values to zero
m.jwf_i3 <- map(alist(
  score ~ dnorm(mu, sigma),
  mu <- a + bJ*judge.amer + bW*wine.amer +
    #bF*is_red +  
    #bJW*judge.amer*wine.amer + 
    bJF*judge.amer*is_red + bWF*wine.amer*is_red, #+
    #bJWF*(judge.amer)*(wine.amer)*(is_red),
  a ~ dnorm(14, 20),
  bJ ~ dnorm(0, .2),
  bW ~ dnorm(0, .2),
  #bF ~ dnorm(0, .2),
  #bJW ~ dnorm(0, .2),
  bJF ~ dnorm(0, .2),
  bWF ~ dnorm(0, .2),
  #bJWF ~ dnorm(0, .2),
  sigma ~ dunif(0, 100)
), 
data=d)
compare(m.wf, m.wf_i, m.jwf_iwf, m.jwf_i3)
precis(m.jwf_i3)


# compare red wines amer vs non-amer using link function and best model (jwf_i3)
d.predict <- data.frame(wine.amer=c(1,1,0,0), is_red=c(1,1,1,1), judge.amer=c(1,0,1,0))
mu <- link(m.jwf_i3, data=d.predict, n = 1e+4)
str(mu)
mu.red.amer <- c(mu[,1], mu[,2])
mu.red.nonamer <- c(mu[,3], mu[,4])
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)
diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)

# compare red wines amer vs non-amer using ensemble model (m.wf, m.wf_i, m.jwf_iwf, m.jwf_i3)
d.predict <- data.frame(wine.amer=c(1,1,0,0), is_red=c(1,1,1,1), judge.amer=c(1,0,1,0))
m.ens <- ensemble(m.wf, m.wf_i, m.jwf_iwf, m.jwf_i3, data = d.predict, n = 1e+4)
mu <- m.ens$link #link(m.wf_i, data=d.predict, n = 1e+4)
str(mu)
mu.red.amer <- c(mu[,1], mu[,2])
mu.red.nonamer <- c(mu[,3], mu[,4])
dens(mu.red.amer, col='red', xlim=c(12,16.5))
dens(mu.red.nonamer, col='blue', add=T)
diff.red.region <- mu.red.amer - mu.red.nonamer
dens(diff.red.region)
str(diff.red.region)
mean(diff.red.region<0)

