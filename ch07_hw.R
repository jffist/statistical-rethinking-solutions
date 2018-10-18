library(rethinking)
library(ggplot2)
library(GGally)
library(dplyr)
#source("plot_bindings.R")

### Easy ####
# 7E1
# (1) Bread dough rises because of yeast.
# -- and amount of flour/eggs/"soure cream"/water (any other ingredients of the dough)
#
# (2) Education leads to higher income.
# -- but also <parents income> does matter, or just country. In different countries slope between education and income is be different

# (3) Gasoline makes a car go.
# in this example it's not obvious what is the predicted variable - let's assume it's speed
# -- volume of the engine or number of cylinders, or proxy variable like a model of the car
#

# 7E2+7E3. 
# Which of the following explanations invokes an interaction?
# For each of the explanations in 7E2, write a linear model that expresses the stated relationship
#  (1) Caramelizing onions requires cooking over low heat and making sure the onions do not dry out.
#      caramelizing ~ heat + dryness + heat*dryness
#      One can expect that for the higher temperature there are more chances for the onion to dry out, so I added interaction of those terms.
#  (2) A car will go faster when it has more cylinders or when it has a better fuel injector.
#      Both models - with and without interaction look reasonable for me.
#      max_speed ~ n_cylinders + fuel_injector
#      max_speed ~ n_cylinders + fuel_injector + n_cylinders*fuel_injector
#  (3) Most people acquire their political beliefs from their parents, unless they get them instead from their friends.
#      political_beliefs ~ parents_beliefs*(1-has_political_friends) + has_political_friends*friends_beliefs

#  (4) Intelligent animal species tend to be either highly social or have manipulative appendages (hands, tentacles, etc.).
#      Without any domain knowledge, I assume, that there is no interaction between social level and presence of manipulative appendages as evidence of intelligence.
#      intelligence ~ social_level + has_manipulalative_appendages 

### Medium ####
# 7M1
# Temperature is the main factor, shade and water depends on it. 
# or there is an interaction between T and W, T and S, T and W and S

# 7M2
# bloom ~ is_temp_cold * (a + bW*w + bS*s + bWS*w*s)

# 7M3
# The simplest dataset contains two variables - "raven population size" and "wolf population size"  measured across different areas. 
# I expect that these variables are correlated. Also, I assume, this relation is non-linear. As with the increase of wolves population but constant amount of food on the area, wolves tend to eat all food not leaving anything to ravens. 
# So hypothetical relation look like  
#  raven_pop_size ~ beta_wolf*wolf_pop_size +beta_wolf2*raven_pop_size^2 with beta_wolf2<0
# Extended data set would consist of 4 variables
#  1) raven population 
#  2) wolves population 
#  3) amount of potential food (herbivorous population in are) 
#  4) area size
# Then I would expect model like 
#  raven_pop_size ~ a1 + (beta_food + beta_wolf*wolf_pop_size)*food, that explicitely states that slope on food depends from wolves population size


### Hard ####
# 7H1 + 7H2 ####
data(tulips)
d <- tulips
str(d)
pairs(d)
d$shade.c <- d$shade - mean(d$shade)
d$water.c <- d$water - mean(d$water)
d$bed.idx <- coerce_index(d$bed)
m7.8 <- map(
  alist(
    blooms ~ dnorm( mu , sigma ) ,
    mu <- a + bW*water.c + bS*shade.c ,
    a ~ dnorm( 130 , 100 ) ,
    bW ~ dnorm( 0 , 100 ) ,
    bS ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 100 )
  ) ,
  data=d ,
  start=list(a=mean(d$blooms),bW=0,bS=0,sigma=sd(d$blooms)) )
m7.9 <- map(
  alist(
    blooms ~ dnorm( mu , sigma ) ,
    mu <- a + bW*water.c + bS*shade.c + bWS*water.c*shade.c ,
    a ~ dnorm( 130 , 100 ) ,
    bW ~ dnorm( 0 , 100 ) ,
    bS ~ dnorm( 0 , 100 ) ,
    bWS ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 100 )
  ) ,
  data=d ,
  start=list(a=mean(d$blooms),bW=0,bS=0,bWS=0,sigma=sd(d$blooms)) )

coeftab(m7.8,m7.9)

m7.10 <- map(
  alist(
    blooms ~ dnorm( mu , sigma ) ,
    mu <- a[bed.idx] + bW*water.c + bS*shade.c + bWS*water.c*shade.c ,
    a[bed.idx] ~ dnorm( 130 , 100 ) ,
    bW ~ dnorm( 0 , 100 ) ,
    bS ~ dnorm( 0 , 100 ) ,
    bWS ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 100 )
  ) ,
  data=d ,
  start=list(a=c(mean(d$blooms),mean(d$blooms),mean(d$blooms)),bW=0,bS=0,bWS=0,sigma=sd(d$blooms)),
  method="Nelder-Mead" ,
  control=list(maxit=1e4))

coeftab(m7.8,m7.9, m7.10)
compare(m7.8,m7.9, m7.10)

samples <- extract.samples(m7.10)
str(samples)

hist(samples$a[,1], xlim=c(0,200), ylim=c(0, 1500), col=col.alpha('blue',0.3))
hist(samples$a[,2], add=T, col=col.alpha('red',0.3))
hist(samples$a[,3], add=T, col=col.alpha('yellow',0.3))

dens(samples$a[,1], xlim=c(0,200),  col='blue')
dens(samples$a[,2], add=T, col='red')
dens(samples$a[,3], add=T, col='black')

dens(samples$a[,1]-samples$a[,2], col='blue')
dens(samples$a[,1]-samples$a[,3], add=T, col='red')
dens(samples$a[,2]-samples$a[,3], add=T, col='black')

d$is_first_bed <- as.integer(d$bed.idx!=1)+1
m7.11 <- map(
  alist(
    blooms ~ dnorm( mu , sigma ) ,
    mu <- a[is_first_bed] + bW*water.c + bS*shade.c + bWS*water.c*shade.c ,
    a[is_first_bed] ~ dnorm( 130 , 100 ) ,
    bW ~ dnorm( 0 , 100 ) ,
    bS ~ dnorm( 0 , 100 ) ,
    bWS ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 100 )
  ) ,
  data=d ,
  start=list(a=c(mean(d$blooms),mean(d$blooms)),bW=0,bS=0,bWS=0,sigma=sd(d$blooms)),
  method="Nelder-Mead" ,
  control=list(maxit=1e4))
# Several executions of compare function below gave slightly different resutls for models' weights(!!!). 
# Weight of the top model varies from 0.85 to 0.96
compare(m7.8, m7.9, m7.11)
coeftab(m7.8, m7.9, m7.10, m7.11)

# add interaction of water and shade with bed
m7.12 <- map(
  alist(
    blooms ~ dnorm( mu , sigma ) ,
    mu <- a[bed.idx] + bW[bed.idx]*water.c + bS[bed.idx]*shade.c + bWS*water.c*shade.c ,
    a[bed.idx] ~ dnorm( 130 , 100 ) ,
    bW[bed.idx] ~ dnorm( 0 , 100 ) ,
    bS[bed.idx] ~ dnorm( 0 , 100 ) ,
    bWS ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 100 )
  ) ,
  data=d ,
  start=list(a=c(mean(d$blooms),mean(d$blooms),mean(d$blooms)),bW=0,bS=0,bWS=0,sigma=sd(d$blooms)),
  method="Nelder-Mead" ,
  control=list(maxit=1e4))
compare(m7.8,m7.9, m7.11, m7.12)
precis(m7.12, depth=2)


# 7H3 ####
data("rugged")
d <- rugged
d$log_gdp <- log( d$rgdppc_2000 )
d <- d[complete.cases(d$rgdppc_2000), ]
nrow(d)

d$rugged.c <- d$rugged - mean(d$rugged) 
slope.sigma <- 10 #(!) try to make several runs varying this value(e.g. 0.1 correspondts to regularized priors and difference between full and ns model becomes smaller)
m7h3.full <- map(alist(
          log_gdp ~ dnorm(mu, sigma),
          mu <- a + bR*rugged.c + bA*cont_africa + rugged.c*cont_africa*bAR,
          a ~ dnorm(8, 100),
          bR ~ dnorm(0, slope.sigma),
          bA ~ dnorm(0, slope.sigma),
          bAR ~ dnorm(0, slope.sigma),
          sigma ~ dunif(0, 10)
     ),
     data=d) 
precis(m7h3.full)

## plot predictions for Africa and not Africa
#### helper functions ####
predict_mu <- function(model, d.predict){
  mu <- link(model, data=d.predict)
  mu.mean <- apply(mu, 2, mean)
  mu.pi <- apply(mu, 2, PI)
  list(mean=mu.mean, pi=mu.pi)  
}
plot_model_mu <- function(d.raw, d.predict, mu.pred, title){
  plot(log_gdp ~ rugged.c, data=d.raw, col='blue', ylim=c(5,12), xlim=c(-2, 6))
  lines(d.predict$rugged.c, mu.pred$mean, col='red')
  shade(mu.pred$pi, d.predict$rugged.c)
  mtext(title)
}
d.predict.af <- data.frame(rugged.c=seq(-2,6,by=0.1), cont_africa=1)
d.predict.naf <- data.frame(rugged.c=seq(-2,6,by=0.1), cont_africa=0)

mu7h3.full.af <- predict_mu(m7h3.full, d.predict.af)
mu7h3.full.naf <- predict_mu(m7h3.full, d.predict.naf)

par(mfrow=c(1,2))
plot_model_mu(d[d$cont_africa==1,], d.predict.af,  mu7h3.full.af, 'Africa, full')
plot_model_mu(d[d$cont_africa==0,], d.predict.naf,  mu7h3.full.naf, 'not Africa, full')

# fit model on data without Seychelles
d.ns <- d[d$country!='Seychelles',]
#slope.sigma <- 10
m7h3.ns <- map(alist(
  log_gdp ~ dnorm(mu, sigma),
  mu <- a + bR*rugged.c + bA*cont_africa + rugged.c*cont_africa*bAR,
  a ~ dnorm(8, 100),
  bR ~ dnorm(0, slope.sigma),
  bA ~ dnorm(0, slope.sigma),
  bAR ~ dnorm(0, slope.sigma),
  sigma ~ dunif(0, 10)
),
data=d.ns) 
precis(m7h3.ns)
mu7h3.ns.af <- predict_mu(m7h3.ns, d.predict.af)
mu7h3.ns.naf <- predict_mu(m7h3.ns, d.predict.naf)

par(mfrow=c(1,2))
plot_model_mu(d.ns[d.ns$cont_africa==1,], d.predict.af,  mu7h3.ns.af, 'Africa, no Seychelles')
plot_model_mu(d.ns[d.ns$cont_africa==0,], d.predict.naf,  mu7h3.ns.naf, 'not Africa, no Seychelles')

plot(coeftab(m7h3.full, m7h3.ns))
#compare(m7h3.full, m7h3.ns)

par(mfrow=c(1,1))
#plot Africa's only predictions
plot_model_mu(d[d$cont_africa==1,], d.predict.af,  mu7h3.full.af, 'Africa full (red+grey) vs Africa withou Seychelles(yellow+green)')
lines(d.predict.af$rugged.c, mu7h3.ns.af$mean, col='yellow')
shade(mu7h3.ns.af$pi, d.predict.af$rugged.c, col=col.alpha('green'))

# compare distributions of slope for Africa countries in two models
m7h3.full.sample <- extract.samples(m7h3.full)
m7h3.ns.sample <- extract.samples(m7h3.ns)

m7h3.full.sample <- mutate(m7h3.full.sample, af_slope=bR+bAR)
m7h3.ns.sample <- mutate(m7h3.ns.sample, af_slope=bR+bAR)

summary(m7h3.full.sample$af_slope)
mean(m7h3.full.sample$af_slope>0)

summary(m7h3.ns.sample$af_slope)
mean(m7h3.ns.sample$af_slope>0)


summary(m7h3.full.sample$bR)
summary(m7h3.ns.sample$bR)

dens(m7h3.full.sample$af_slope, xlim=c(-0.3, 0.5),   col='blue')
dens(m7h3.ns.sample$af_slope, add=T, col='red')
dens(m7h3.full.sample$af_slope-m7h3.ns.sample$af_slope, add=T, col='black')#not sure it's a correct procedure

# summary for 7H3 (a+b)
# Without Seychelles data point slope for African countries reduces almost twice but is still positive compared to non-African countries that have a negative slope.
# The use of strong regularised priors with sigma = 0.1 decreases all slopes and make differences between full and NoSeychelles(NS) models almost neglectable.
# Even for NS model there are 0.77 probability that slope for African countries is greater than zero.

# 7H3.c
slope.sigma <- 10
m7h3.ns.r <- map(alist(
      log_gdp ~ dnorm(mu, sigma),
      mu <- a + bR*rugged.c,
      a ~ dnorm(8, 100),
      bR ~ dnorm(0, slope.sigma),
      sigma ~ dunif(0, 10)
  ),
  data=d.ns) 

m7h3.ns.ra <- map(alist(
    log_gdp ~ dnorm(mu, sigma),
    mu <- a + bR*rugged.c + bA*cont_africa,
    a ~ dnorm(8, 100),
    bR ~ dnorm(0, slope.sigma),
    bA ~ dnorm(0, slope.sigma),
    sigma ~ dunif(0, 10)
  ),
  data=d.ns, 
  start=list(a=c(mean(d.ns$log_gdp)),bR=0,bA=0,sigma=sd(d.ns$log_gdp))
  ) 

m7h3.ns <- map(alist(
  log_gdp ~ dnorm(mu, sigma),
  mu <- a + bR*rugged.c + bA*cont_africa + rugged.c*cont_africa*bAR,
  a ~ dnorm(8, 100),
  bR ~ dnorm(0, slope.sigma),
  bA ~ dnorm(0, slope.sigma),
  bAR ~ dnorm(0, slope.sigma),
  sigma ~ dunif(0, 10)
),
data=d.ns) 

compare(m7h3.ns.r, m7h3.ns.ra, m7h3.ns)

mu7h3.ns.af <- predict_mu(m7h3.ns, d.predict.af)
mu7h3.ns.naf <- predict_mu(m7h3.ns, d.predict.naf)

par(mfrow=c(3,2))
d.predict.af <- data.frame(rugged.c=seq(-2,6,by=0.1), cont_africa=1)
d.predict.naf <- data.frame(rugged.c=seq(-2,6,by=0.1), cont_africa=0)
for(model in list(m7h3.ns.r, m7h3.ns.ra, m7h3.ns)){
  mu.af <- predict_mu(model, d.predict.af)
  plot_model_mu(d.ns[d.ns$cont_africa==1,], d.predict.af,  mu.af, 'Africa, no Seychelles')
  
  m.naf <- predict_mu(model, d.predict.naf)
  plot_model_mu(d.ns[d.ns$cont_africa==0,], d.predict.naf,  m.naf, 'not Africa, no Seychelles')
}
par(mfrow=c(1,1))
coeftab(m7h3.ns.r, m7h3.ns.ra, m7h3.ns)
compare(m7h3.ns.r, m7h3.ns.ra, m7h3.ns)

### averaged across models
mu.ensemble <- ensemble(m7h3.ns.r, m7h3.ns.ra, m7h3.ns, data = d.predict.af)
mu.mean <- apply(X = mu.ensemble$link, MARGIN = 2, FUN = mean)
mu.PI <- apply(X = mu.ensemble$link, MARGIN = 2, FUN = PI)
plot_model_mu(d.ns[d.ns$cont_africa==1,], d.predict.af,  list(mean=mu.mean, pi=mu.PI), 'Africa, no Seychelles')

#### 7H4 ####
data(nettle)
d <- tbl_df(nettle)
d <- d %>% mutate(
        lang.per.cap = num.lang / k.pop,
        log_pop = log(k.pop),
        log_area = log(area),
        log_lang.per.cap = log(lang.per.cap),
        pop_density = k.pop/area,
        
        log_area_c = log_area-mean(log_area),
        mgc_c = mean.growing.season - mean(mean.growing.season)
)
d
ggpairs(select(d,-country))
d <- as.data.frame(d)


dens(log(d$lang.per.cap))
hist(log(d$lang.per.cap),breaks = 20)

hist(d$mgc_c)
hist(d$log_area_c)

slope.sigma <- 10
m7h4.1 <- map(alist(
        log_lang.per.cap ~ dnorm(mu, sigma),
        mu <- a + b.ms*mgc_c + b.area*log_area_c,
        a ~ dnorm(-5, 100),
        b.ms ~ dnorm(0, slope.sigma),
        b.area ~ dnorm(0, slope.sigma),
        sigma~dunif(0,10)
    ), 
    data=d)
precis(m7h4.1)

# I'm just curious how priors influence coefficients
slope.sigma <- .1
m7h4.1rp <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.ms*mgc_c + b.area*log_area_c,
  a ~ dnorm(-5, 100),
  b.ms ~ dnorm(0, slope.sigma),
  b.area ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d)
precis(m7h4.1rp)
compare(m7h4.1, m7h4.1rp)
coeftab(m7h4.1, m7h4.1rp)

hist(d$sd.growing.season)
d$sdgs_c <- d$sd.growing.season - mean(d$sd.growing.season)

m7h4.2 <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.sds*sdgs_c + b.area*log_area_c,
  a ~ dnorm(-5, 100),
  b.sds ~ dnorm(0, slope.sigma),
  b.area ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d)
precis(m7h4.2)
compare(m7h4.1, m7h4.2)

m7h4.3 <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.ms*mgc_c + b.sds*sdgs_c,# + b.area*log_area_c, 
  a ~ dnorm(-5, 100),
  b.ms ~ dnorm(0, slope.sigma),
  b.sds ~ dnorm(0, slope.sigma),
#  b.area ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d)
precis(m7h4.3)
compare(m7h4.1, m7h4.2, m7h4.3)


m7h4.3i <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.ms*mgc_c + b.sds*sdgs_c + b.ms.sds*mgc_c*sdgs_c, 
  a ~ dnorm(-5, 100),
  b.ms ~ dnorm(0, slope.sigma),
  b.sds ~ dnorm(0, slope.sigma),
  b.ms.sds ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d)
precis(m7h4.3i)
compare(m7h4.1, m7h4.2, m7h4.3, m7h4.3i)

##
plot_model_mu_mgc <- function(d.raw, d.predict, mu.pred, title){
  plot(log_lang.per.cap ~ mgc_c, data=d.raw, col='blue', ylim=c(-10,0), xlim=c(-7.1, 5))
  lines(d.predict$mgc_c, mu.pred$mean, col='red')
  shade(mu.pred$pi, d.predict$mgc_c)
  mtext(title)
}
plot_model_mu_sds <- function(d.raw, d.predict, mu.pred, title){
  plot(log_lang.per.cap ~ sdgs_c, data=d.raw, col='blue', ylim=c(-10,0), xlim=c(-2, 4.1))
  lines(d.predict$sdgs_c, mu.pred$mean, col='red')
  shade(mu.pred$pi, d.predict$sdgs_c)
  mtext(title)
}

## tryptych for ~mgc_c
summary(d$sdgs_c)
par(mfrow=c(1,4))
mgc_c.seq <- seq(min(d$mgc_c),max(d$mgc_c),by=0.1)
model <- m7h4.3i
model.name <- ' m7h4.3i'
left_q <- 0
for(right_q in c(0.25, 0.5, 0.75, 1)){
  qq <- quantile(d$sdgs_c, probs=c(left_q, right_q))
  lval = qq[1]
  rval = qq[2]
  if(right_q!=1){
    d.raw <- d[(d$sdgs_c >= lval) & (d$sdgs_c < rval),]
  } else {
    d.raw <- d[(d$sdgs_c >= lval) & (d$sdgs_c <= rval),]
  }
  sdg_c = mean(d.raw$sdgs_c)
  d.predict=data.frame(mgc_c=mgc_c.seq, sdgs_c=sdg_c, log_area_c=0)
  mu <- predict_mu(model, d.predict)
  plot_model_mu_mgc(d.raw, d.predict, mu, paste0(c(model.name,'sdgs_c=',sdg_c),collapse = ' '))
  left_q <- right_q
}




sdgs_c.seq <- seq(min(d$sdgs_c),max(d$sdgs_c),by=0.1)
model <- m7h4.2
model.name <- 'm7h4.2'
for(mgc_c in c(-7, 0, 4.9)){
  d.predict=data.frame(mgc_c=mgc_c, sdgs_c=sdgs_c.seq, log_area_c=0)
  mu <- predict_mu(model, d.predict)
  plot_model_mu_sds(d, d.predict, mu, paste0(c(model.name,'mgc_c=',mgc_c),collapse = ' '))
}



###### add population density to models and compare them with already fitted #####
m7h4.4d <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.ms*mgc_c + b.dens*log(pop_density),
  a ~ dnorm(-5, 100),
  b.ms ~ dnorm(0, slope.sigma),
  b.dens ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d, method="Nelder-Mead" ,
control=list(maxit=1e4))
precis(m7h4.4d)

m7h4.3i.d <- map(alist(
  log_lang.per.cap ~ dnorm(mu, sigma),
  mu <- a + b.ms*mgc_c + b.sds*sdgs_c + b.ms.sds*mgc_c*sdgs_c+ b.dens*log(pop_density),
  a ~ dnorm(-5, 100),
  b.ms ~ dnorm(0, slope.sigma),
  b.sds ~ dnorm(0, slope.sigma),
  b.ms.sds ~ dnorm(0, slope.sigma),
  b.dens ~ dnorm(0, slope.sigma),
  sigma~dunif(0,10)
), 
data=d, 
method="Nelder-Mead" ,
control=list(maxit=1e5))
precis(m7h4.3i.d)


compare(m7h4.1, m7h4.4d, m7h4.3, m7h4.3i, m7h4.3i.d)
plot(coeftab(m7h4.1, m7h4.4d, m7h4.3, m7h4.3i, m7h4.3i.d))
