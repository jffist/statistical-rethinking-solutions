library(rethinking)
library(MASS) 
library(ggplot2)
library(GGally)
library(dplyr)
library(htmltab)
library(stringr)
library(tidyr)
#source("plot_bindings.R") 

## Easy ####

# 5E3 (related)
# I am curious about generating dataset with similar properties as described in task
# I'm looking for dependent variable z that 
# z ~ a*x+b*y -> a>>0 and b>>0 (significantly greater then zero)
# but
# z ~ a*x --> a ~ 0 (a is close to zero)
# z ~ b*y --> b ~ 0
N <- 1000
rho <- 0.8
set.seed(11)
x <- rnorm(N, 10)
y <- rnorm(N, -rho*x , sqrt(1-rho^2) )
d <- data.frame(x=x, y=y, z=rnorm(1000, x + y, 0.5))
ggpairs(d)
summary(lm(z~y, d))
summary(lm(z~x, d))
summary(lm(z~x+y, d))

m5e3x <- map(
  alist(
    z ~ dnorm(mu, sigma),
    mu <- a + bx*x,
    a ~ dnorm(0,10),
    bx ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data = d
) 
summary(m5e3x)
plot(z~x,d, col = col.alpha("blue",0.5))
x.seq <- seq(7, 14, by=0.1)
mu.5e3x <- link(m5e3x, data=list(x=x.seq), n=1000)
mu.5e3x.mean <- apply(mu.5e3x, 2, mean)
mu.5e3x.pi <- apply(mu.5e3x, 2, PI)
lines(x.seq, mu.5e3x.mean, col='red')
shade(mu.5e3x.pi, x.seq, col = col.alpha("red",0.15))

m5e3y <- map(
  alist(
    z ~ dnorm(mu, sigma),
    mu <- a + by*y,
    a ~ dnorm(0,10),
    by ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data = d
) 
summary(m5e3y)
plot(z~y,d, col = col.alpha("blue",0.5))
y.seq <- seq(-14, -4, by=0.1)
mu.5e3y <- link(m5e3y, data=list(y=y.seq), n=1000)
mu.5e3y.mean <- apply(mu.5e3y, 2, mean)
mu.5e3y.pi <- apply(mu.5e3y, 2, PI)
lines(y.seq, mu.5e3y.mean, col='red')
shade(mu.5e3y.pi, y.seq, col = col.alpha("red",0.15))

m5e3xy <- map(
  alist(
    z ~ dnorm(mu, sigma),
    mu <- a + bx*x + by*y,
    a ~ dnorm(0,10),
    bx ~ dnorm(0, 1),
    by ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data = d
) 
summary(m5e3xy)

# add to the previous plot of z~y mean of relation between z and y in the new model(5e3xy) 
mu.5e3xy.y <- link(m5e3xy, data=data.frame(y=y.seq, x=mean(d$x)), n=1000)
mu.5e3xy.y.mean <- apply(mu.5e3xy.y, 2, mean)
mu.5e3xy.y.pi <- apply(mu.5e3xy.y, 2, PI)
lines(y.seq, mu.5e3xy.y.mean, col='green')
shade(mu.5e3xy.y.pi, y.seq, col = col.alpha("green",0.5))

# 5E4
# lets explore inferential equivalency of m1 and m3 using milk data
data(milk) 
d <- milk
d$catvar <- as.factor(d$clade)
levels(d$catvar) <- c('A','B','C','D')

(d$catvar.A <- as.integer(d$catvar=='A'))
(d$catvar.B <- as.integer(d$catvar=='B'))
(d$catvar.C <- as.integer(d$catvar=='C'))
(d$catvar.D <- as.integer(d$catvar=='D'))

ggplot(d, aes(y=kcal.per.g, x=catvar, fill=catvar)) + geom_boxplot()


m5e4.m1 <- map(
  alist(
    kcal.per.g ~ dnorm( mu , sigma ) ,
    mu <- a + b.a*catvar.A + b.b*catvar.B + b.d*catvar.D, #a is mean val for level C
    a ~ dnorm( 0.6 , 10 ) ,
    b.a ~ dnorm( 0 , 1 ) ,
    b.b ~ dnorm( 0 , 1 ) ,
    b.d ~ dnorm( 0 , 1 ) ,
    sigma ~ dunif( 0 , 10 )
  ) ,
  data=d )
precis(m5e4.m1)

post.5e4.m1 <- extract.samples(m5e4.m1, n=1e3)
head(post.5e4.m1)

m5e4.m3 <- map(
  alist(
    kcal.per.g ~ dnorm( mu , sigma ) ,
    mu <- a + b.b*catvar.B + b.c*catvar.C + b.d*catvar.D, #a is mean val for level A
    a ~ dnorm( 0.6 , 10 ) ,
    b.b ~ dnorm( 0 , 1 ) ,
    b.c ~ dnorm( 0 , 1 ) ,
    b.d ~ dnorm( 0 , 1 ) ,
    sigma ~ dunif( 0 , 10 )
  ) ,
  data=d )
precis(m5e4.m3)

post.5e4.m3 <- extract.samples(m5e4.m3, n=1e3)
head(post.5e4.m3)

# transform m1 into m3 and compare histograms of params
post.m3from1 <- with(post.5e4.m1,
  data.frame(
  a = a + b.a,
  b.b = b.b - b.a,
  b.c = -b.a,
  b.d = b.d - b.a,
  sigma=sigma,
  source='3from1'
  )
)
post.5e4.m3$source='m3'

dd <- bind_rows(post.m3from1 , post.5e4.m3)

ggplot(dd, aes(x=a, fill=source)) + 
  geom_density(alpha=0.5)
#  geom_histogram(alpha=0.5, position="identity", bins=20) 

ggplot(dd, aes(x=b.b, fill=source)) + 
  geom_density(alpha=0.5)
#  geom_histogram(alpha=0.5, position="identity", bins=20) 

ggplot(dd, aes(x=sigma, fill=source)) + 
  geom_density(alpha=0.5)

# let's check what happens if we add all indicators and intercept to the model
# intercept now is correspondent to unexistent examples with no categories
m5e4.m2 <- map(
  alist(
    kcal.per.g ~ dnorm( mu , sigma ) ,
    mu <- a + b.a*catvar.A + b.b*catvar.B + b.c*catvar.C + b.d*catvar.D, 
    a ~ dnorm( 0 , 10 ) ,
    b.a ~ dnorm( 0 , 1 ) ,
    b.b ~ dnorm( 0 , 1 ) ,
    b.c ~ dnorm( 0 , 1 ) ,
    b.d ~ dnorm( 0 , 1 ) ,
    sigma ~ dunif( 0 , 10 )
  ) ,
  data=d)
precis(m5e4.m2)

post.5e4.m2 <- extract.samples(m5e4.m2, n=1e3)
head(post.5e4.m2)

# still possible to convert to m3
post.m3from2 <- with(post.5e4.m2,
                     data.frame(
                       a = a + b.a,
                       b.b = b.b - b.a,
                       b.c = b.c - b.a,
                       b.d = b.d - b.a,
                       sigma=sigma,
                       source='3from2'
                     ))
dd <- bind_rows(post.m3from2 , post.5e4.m3)

ggplot(dd, aes(x=a, fill=source)) + 
  geom_density(alpha=0.5)
#  geom_histogram(alpha=0.5, position="identity", bins=20) 

ggplot(dd, aes(x=b.b, fill=source)) + 
  geom_density(alpha=0.5)
#  geom_histogram(alpha=0.5, position="identity", bins=20) 

ggplot(dd, aes(x=sigma, fill=source)) + 
  geom_density(alpha=0.5)

#### Medium ####
# 5M1
# example give in the book looks like
# x -> x_sp (x is used to generate x_sp as correlated)
# x -> y (x is used to generate y as correlated)
# then
# x_sp is correlated to y
# but in presense of x, dependency between y and x_sp is almost degraded
# Another example that can be adopted to solve the task is masked treatment
# (t)reatment -> (c)onsequences -> (r)esult
# t ~ c
# c ~ r
# ==> t ~ r, but modelling r ~ t+c should make beta_t close to zero (degrade order of dependency)
N = 1000
r1 = 0.7
r2 = 0.9
set.seed(112)
x <- rnorm(N, 5)
y <- rnorm(N, r1*x , sqrt(1-r1^2) )
cor(x,y)
z <- rnorm(N, r2*y , sqrt(1-r2^2) )
cor(x,z)
cor(y,z)
d <- data.frame(x=x, y=y, z=z)
pairs(d)
confint(lm(z~x, d))
confint(lm(z~y, d))
confint(lm(z~y+x, d))

# 5M2
# obvious solution is same as in book
# x - any random
# y - correlated with x
# z = x-y
N <- 1000
rho <- 0.8
set.seed(11)
x <- rnorm(N, 1)
y <- rnorm(N, rho*x , sqrt(1-rho^2) )
z <- rnorm(N, x - y, 0.5)
d <- data.frame(x=x, y=y, z=z)
ggpairs(d)
summary(lm(z~y, d))
summary(lm(z~x, d))
summary(lm(z~x+y, d))

# 5M4
data(WaffleDivorce)
d <- WaffleDivorce

# data about LDS per state took from https://www.worldatlas.com/articles/mormon-population-by-state.html
q <- htmltab('https://www.worldatlas.com/articles/mormon-population-by-state.html',
             which="//div[@id='artReg-table']/table")
stopifnot("Percentage of Mormon Residents" %in% names(q))
stopifnot(nrow(q) == 51)
names(q)[2] <- 'Location'

lcd <- q %>% select(Location, 'Percentage of Mormon Residents') %>% 
        rename(ratio_str='Percentage of Mormon Residents') %>%
        mutate(lcd_ratio=as.double(str_replace(ratio_str, "[%]", "")))
d2 <- left_join(d, select(lcd, Location, lcd_ratio), by='Location')

normalise <- function(x){
  (x-mean(x))/sd(x)
}
d2 <- d2 %>% mutate(
      marriage.s=normalise(Marriage),
      age.marriage.s=normalise(MedianAgeMarriage),
      lcd.ratio.s=normalise(lcd_ratio)
)
m1 <- lm(Divorce ~ marriage.s + age.marriage.s, data=d2)
summary(m1)
confint(m1)
m2 <- lm(Divorce ~ marriage.s + age.marriage.s + lcd.ratio.s, data=d2)         
summary(m2)
confint(m2)

predict_with_conf <- function(model, data, name){
  mu <- link( model )
  # summarize samples across cases
  mu.mean <- apply( mu , 2 , mean )
  mu.PI <- apply( mu , 2 , PI )
  data.frame(divorce=data$Divorce, 
             prediction=mu.mean,
             pred.low = mu.PI[1,],
             pred.high = mu.PI[2,],
             model=name
  )
}
res <- bind_rows(
  predict_with_conf(m1, d2, 'm1'),
  predict_with_conf(m2, d2, 'm2.w.lcd')
)

pd <- position_dodge(0.1)
ggplot(res, aes(x=divorce, y=prediction, color=model, fill=model)) +
  geom_point(alpha=0.7, size=2, shape=21, color='black', position=pd) + 
  geom_abline(slope=1,intercept = 0) + 
  geom_errorbar(aes(ymin=pred.low, ymax=pred.high), width=.5, position=pd) +
  ylim(4, 16) + xlim(4,16) + xlab('Observed divorce')

### Hard ####
data(foxes)
d <- foxes
str(d)
unique(d$group)
gd <- d %>% group_by(group) %>% summarise(
  avgfood=first(avgfood),
  groupsize=first(groupsize),
  area=first(area),
  weight.min = min(weight),
  weight.avg = mean(weight),
  weight.max = max(weight)
)
ggpairs(gd)

ggpairs(d)

### 5H1
ggplot(d, aes(area, weight)) + geom_point(size=2) + geom_smooth(method='lm')
ggplot(d, aes(groupsize, weight)) + geom_point(size=2) + geom_smooth(method='lm')

### area
m5h1.a <- map(alist(
      weight ~ dnorm( mu , sigma ),
      mu <- Intercept + b_area*area,
      Intercept ~ dnorm(0,10),
      b_area ~ dnorm(0,10),
      sigma ~ dunif(0,10)
   ),
   data=d
)
summary(m5h1.a)

area.seq <- seq(1,6,0.1)
mu5h1.a <- link(m5h1.a, data=list(area=area.seq))
mu5h1.a.mean <- apply(mu5h1.a, 2, mean)
mu5h1.a.PI <- apply(mu5h1.a, 2, PI, .89)

plot(weight~area, d, col='blue')
lines(area.seq, mu5h1.a.mean, col='red')
shade(mu5h1.a.PI, area.seq)

### groupsize
m5h1.g <- map(alist(
  weight ~ dnorm( mu , sigma ),
  mu <- Intercept + b_gs*groupsize,
  Intercept ~ dnorm(0,10),
  b_gs ~ dnorm(0,10),
  sigma ~ dunif(0,10)
),
data=d
)
summary(m5h1.g)

gs.seq <- seq(1,9,0.1)
mu5h1.g <- link(m5h1.g, data=list(groupsize=gs.seq))
mu5h1.g.mean <- apply(mu5h1.g, 2, mean)
mu5h1.g.PI <- apply(mu5h1.g, 2, PI, .89)

plot(weight~groupsize, d, col='blue')
lines(gs.seq, mu5h1.g.mean, col='red')
shade(mu5h1.g.PI, gs.seq)

# 5H2
m5h2 <- map(alist(
  weight ~ dnorm( mu , sigma ),
  mu <- Intercept + b_gs*groupsize + b_area*area,
  Intercept ~ dnorm(0,10),
  b_gs ~ dnorm(0,10),
  b_area ~ dnorm(0,10),
  sigma ~ dunif(0,10)
),
data=d
)
summary(m5h2)
# 2 area
mu5h2.a <- link(m5h2, data=data.frame(groupsize=mean(d$groupsize), area=area.seq))
mu5h2.a.mean <- apply(mu5h2.a, 2, mean)
mu5h2.a.PI <- apply(mu5h2.a, 2, PI, .89)

plot(weight~area, d, col='blue')
lines(area.seq, mu5h1.a.mean, col='red')
shade(mu5h1.a.PI, area.seq,  col = col.alpha("goldenrod",0.5))
lines(area.seq, mu5h2.a.mean, col='cyan3')
shade(mu5h2.a.PI, area.seq, col = col.alpha("cyan3",0.5))

# 2 groupsize
mu5h2.g <- link(m5h2, data=data.frame(area=mean(d$area), groupsize=gs.seq))
mu5h2.g.mean <- apply(mu5h2.g, 2, mean)
mu5h2.g.PI <- apply(mu5h2.g, 2, PI, .89)

plot(weight~groupsize, d, col='blue')
lines(gs.seq, mu5h1.g.mean, col='red')
shade(mu5h1.g.PI, gs.seq, col = col.alpha("goldenrod",0.5))
lines(gs.seq, mu5h2.g.mean, col='cyan3')
shade(mu5h2.g.PI, gs.seq, col = col.alpha("cyan3",0.5))

# 5H3
# quick check
d <- mutate(d, 
            avgfood.s = normalise(avgfood),
            area.s = normalise(area),
            groupsize.s=normalise(groupsize))

m.fd.gs <- lm(weight ~ avgfood.s + groupsize.s, d)
summary(m.fd.gs)
precis(m.fd.gs, digits=3)

m.fd.gs.a <- lm(weight ~ avgfood.s + groupsize.s + area.s, d)
summary(m.fd.gs.a)
precis(m.fd.gs.a, digits=3)

m.gs.a <- lm(weight ~  groupsize.s + area.s, d)
summary(m.gs.a)
precis(m.gs.a, digits=3)


# group size from area and food
m.f.a <- lm(groupsize ~ avgfood.s + area.s, d)
summary(m.f.a)
precis(m.f.a, digits=3)

gd <- mutate(gd, 
             avgfood.s = normalise(avgfood),
             area.s = normalise(area),
             groupsize.s=normalise(groupsize))
m.f.a2 <- lm(groupsize ~ avgfood.s + area.s, gd)
summary(m.f.a2)
precis(m.f.a2, digits=3)
