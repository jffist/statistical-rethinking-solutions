rm(list=ls())
library(rethinking)

## Easy ####
p_grid <- seq( from=0 , to=1 , length.out=1000 ) 
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

hist(samples)


# 3e1
mean(samples<0.2) #5e-04

mean(samples>0.8) #0.1117

mean(samples>0.2 & samples<0.8) #0.8878

# 3e4 3e5
quantile(samples, c(0.2, 0.8))

# 3e6
HPDI( samples , prob=0.66)

# 3e7
PI( samples , prob=0.66)

## Medium ####
# 3m1
p_grid <- seq( from=0 , to=1 , length.out=1000 ) 
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 8 , size=15 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
plot(p_grid, posterior, type='l')
p_grid[ which.max(posterior) ] #0.5335335

# 3m2
set.seed(123)
samples <- sample(p_grid, prob = posterior, replace = TRUE, size = 1e+4)
hist(samples)
HPDI(samples, prob = .9) #0.3543544 0.7437437
mean(samples) #0.5295707
median(samples) #0.5305305

# 3m3
n = 15
dumdata <- rbinom(10000, size=n, prob=samples)
simplehist(dumdata)
mean(dumdata==8) #0.1445

# 3m4
likelihood_6of9 <- dbinom( 6 , size=9 , prob=p_grid )
prior_6of9 <- posterior
(p_6of9 <- sum(likelihood_6of9*prior_6of9)) #0.1763898

# alt check with generative approach
set.seed(22)
dumdata_6of9 <- rbinom(10000, size=9, prob=samples)
simplehist(dumdata_6of9)
mean(dumdata_6of9==6) #0.1735

# 3m5
#### 5.1
p_grid <- seq( from=0 , to=1 , length.out=1000 ) 
prior <- ifelse(p_grid<0.5, 0, 0.5)
likelihood <- dbinom( 8 , size=15 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
plot(p_grid, posterior, type='l')
p_grid[ which.max(posterior) ] # 0.5335335

#### 5.2
set.seed(123)
samples <- sample(p_grid, prob = posterior, replace = TRUE, size = 1e+4)
hist(samples)
HPDI(samples, prob = .9) #0.5005005 0.7137137
mean(samples) #0.6070759
median(samples) #0.5945946

#### 5.3
n = 15
dumdata <- rbinom(10000, size=n, prob=samples)
simplehist(dumdata)
mean(dumdata==8) #0.1576
table(dumdata)/1e+4

#### 5.4
likelihood_6of9 <- dbinom( 6 , size=9 , prob=p_grid )
prior_6of9 <- posterior
(p_6of9 <- sum(likelihood_6of9*prior_6of9)) #0.2323071

# alt check with generative approach
set.seed(22)
dumdata_6of9 <- rbinom(10000, size=9, prob=samples)
simplehist(dumdata_6of9)
mean(dumdata_6of9==6) #0.2319

## Hard ####
data(homeworkch3)
birth1
birth2
length(birth1)

# 3H1
n_boys = sum(c(birth1, birth2))
n_ttl = length(birth1) + length(birth2)
n_pgrid = 1000
p_grid = seq(0, 1, length.out = n_pgrid)
prior = rep(1, n_pgrid)
likelihood = dbinom(n_boys, size=n_ttl, prob = p_grid)
posterior = likelihood * prior
posterior = posterior / sum(posterior)
plot(p_grid, posterior, type='l')
p_grid[which.max(posterior)] #0.5545546

# 3H2
n_ptrials = 1e4
p_samples = sample(p_grid,size = n_ptrials, prob = posterior, replace = TRUE)
(hpi_50 = HPDI(samples, .5))
(hpi_89 = HPDI(samples, .89))
(hpi_97 = HPDI(samples, .97))
for(w in c(.5, .89, .97)){
  hpi = HPDI(samples, w)
  print(sprintf("HPDI %d%% [%f, %f]",w*100, hpi[1], hpi[2]))
}
mean(p_samples)
median(p_samples)

# 3H3
n_btrials = 1e4 #birth observations
b_sample = rbinom(n_btrials, size=n_ttl, prob=p_samples)
simplehist(b_sample)
mean(b_sample) #110.9
median(b_sample) #111
dens(b_sample)
abline(v = n_boys, col = "red")

# 3H4
n_boys_b1 = sum(birth1)
n_ttl_b1 = length(birth1)
n_btrials = 1e4 #birth observations
b_sample = rbinom(n_btrials, size=n_ttl_b1, prob=p_samples)
simplehist(b_sample)
abline(v = n_boys_b1, col = "red")
mean(b_sample) #55.39
median(b_sample) #55
dens(b_sample)
abline(v = n_boys_b1, col = "red")
# model overestimates number of boys for the first child

# 3H5
n_ttl_g1 = sum(birth1==0)
n_ttl_g1b2 = sum(birth2[birth1==0])
print(sprintf("there were %d boys born after first girl. There were ttl %d cases",n_ttl_g1b2,n_ttl_g1))

n_btrials = 1e4 #birth observations
b_sample = rbinom(n_btrials, size=n_ttl_g1, prob=p_samples)
mean(b_sample) #27
median(b_sample) #27
#dens(b_sample)
simplehist(b_sample)
abline(v = n_ttl_g1b2, col = "red")
# model underestimates number of boys for the second child after the first girl

# conclusion - gender of the second child is not independent from the first one

