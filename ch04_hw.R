rm(list=ls())
library(rethinking)
library(ggplot2)
library(MASS)

## Medium ####
# 4M1
trials <- 1e4
mu <- rnorm(trials, 0, 10)
sigma <- runif(trials, 0, 10)
height <- rnorm(trials, mu, sigma)
dens(height)

# 4M2
formula4.2 <- alist(
  y ~ dnorm(mu, sigma),
  mu ~ dnorm(0, 10),
  sigma ~ dunif(0, 10)
)

## Hard ####
# 4H1
### load data
data("Howell1")
d <- Howell1
### check age and height for people with weight > 30
d2 <- d[d$weight>=30,]
hist(d2$age)
summary(d2)
ggplot(d2, aes(weight, height)) + geom_point()
### to simplify life lets use linear regression as in chapter 4 for individuals with age>18
d2 <- d[d$age>=18,]
summary(d2)
ggplot(d2, aes(weight, height)) + geom_point()
### fit map model 
m4h1 <- map(
  alist(
     height ~ dnorm(mu, sigma),
     mu <- a + b*weight,
     a ~ dnorm(155,50),
     b ~ dnorm(0, 10),
     sigma ~ dunif(0, 50)
  ),
  data=d2
)
summary(m4h1)
### create dataframe with records of interest
dnew <- data.frame(id=1:5, weight=c(46.95, 43.72, 64.78, 32.59, 54.63))
sim.dnew <- sim(m4h1, data = dnew, n=1e4)
str(sim.dnew)
dnew$height <- apply(sim.dnew, 2, mean)
hpdi.new <- apply(sim.dnew, 2, HPDI)
dnew$hpdi.low <- hpdi.new[1,]
dnew$hpdi.high <- hpdi.new[2,]
dnew

# draw from posterior for the first individual
post4h1 <- extract.samples(m4h1)
sim4h1.p1 <- rnorm(n = 1e4, mean = post4h1$a + post4h1$b*dnew$weight[1], sd = post4h1$sigma)
(p1.height <- mean(sim4h1.p1))
(p1.hpdi <- HPDI(sim4h1.p1))

# 4H2
d2 <- d[d$age<18,]
stopifnot(nrow(d2)==192)
plot(height~weight, d2)

## 4h2.a
m4h2 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(100,50),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data=d2
)
summary(m4h2)
# for every 10 units increase in weight on average we expect 10*coef(m4h2)['b'] increase in height
10*coef(m4h2)['b']

weight.seq <- seq(round(min(d2$weight)-2), round(max(d2$weight)+2), 1)
mu4h2 <- link(m4h2, data=list(weight=weight.seq), n=1e4)
str(mu4h2)
mu4h2.mean <- apply(mu4h2, 2, mean)
mu4h2.hpdi <- apply(mu4h2, 2, HPDI)

sim4h2 <- sim(m4h2, data=list(weight=weight.seq), n=1e4)
points4h2.hpdi <- apply(sim4h2, 2, HPDI)

plot(height~weight, d2, col=col.alpha("red",0.9))
lines( weight.seq , mu4h2.mean )
shade(mu4h2.hpdi, weight.seq)
shade(points4h2.hpdi, weight.seq, col = col.alpha("green",0.15))

# same stuff calculated manually
post4h2.params <- data.frame( mvrnorm(n=1e4, mu=coef(m4h2), Sigma=vcov(m4h2)) )
head(post4h2.params)
mu4h2.mean.man <- sapply(weight.seq, function(w) mean(post4h2.params$a + post4h2.params$b*w))
mu4h2.hpdi.man <- sapply(weight.seq, function(w) HPDI(post4h2.params$a + post4h2.params$b*w))

height.link <- function(w){
  # generates sample for the given height
  rnorm(n = nrow(post4h2.params), 
        mean = post4h2.params$a + post4h2.params$b*w,
        sd = post4h2.params$sigma)
} 
sim4h2.man <- sapply(weight.seq, height.link)
points4h2.hpdi.man <- apply(sim4h2.man, 2, HPDI)

plot(height~weight, d2, col=col.alpha("red",0.9))
lines( weight.seq , mu4h2.mean.man )
shade(mu4h2.hpdi.man, weight.seq)
shade(points4h2.hpdi.man, weight.seq, col = col.alpha("blue",0.15))

# 4H3
m4h3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*log(weight),
    a ~ dnorm(178,100),
    b ~ dnorm(0, 100),
    sigma ~ dunif(0, 50)
  ),
  data=d
)
summary(m4h3)

weight.seq <- seq(1, 70, 1)
mu4h3 <- link(m4h3, data=list(weight=weight.seq), n=1e4)
str(mu4h3)
mu4h3.mean <- apply(mu4h3, 2, mean, 0.97)
mu4h3.hpdi <- apply(mu4h3, 2, HPDI, 0.97)

sim4h3 <- sim(m4h3, data=list(weight=weight.seq), n=1e4)
points4h3.hpdi <- apply(sim4h3, 2, HPDI, 0.97)

plot(height~weight, d, col=col.alpha("red",0.9))
lines( weight.seq , mu4h3.mean )
shade(mu4h3.hpdi, weight.seq)
shade(points4h3.hpdi, weight.seq, col = col.alpha("blue",0.15))
