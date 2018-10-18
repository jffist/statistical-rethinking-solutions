rm(list=ls())
library(rethinking)
library(rstan)
library(MASS)
library(ellipse)
library(dplyr)
library(ggplot2)
library(tidyr)
#source("plot_bindings.R")

## Easy ####
# 13E1 ####

# y_i ~ Normal(mu[i], sigma)
# mu[i] = a_g[group[i]] + beta_g[group[i]]*x[i]
# c(a_g, beta_g)[group[i]] ~ MVNormal( [a,b], S)
# a ~ Normal(0, 10)
# b ~ Normal(0, 1)
# S = diagonal_sigma_matrix * R * diagonal_sigma_matrix
# sigma_matrix ~ HalfCauchy(0,2) #sigma_matrix is the diagonal matrix of size N_groups, so it has N_groups independent parameters, each distributed as HalfCauchy 
# R ~ LKJcorr(2) #LKJcorr produces priors matrix of the size N_groups x N_groups
# sigma ~ HalfCauchy(0,2)

## 13E2 ####
# case 1)
# Consider a hypothetical example of the average gas mileage of the ride. It depends on the car model (as a proxy of engine's volume) and the speed you travel. 
# The faster you ride, the more gas you need. Also, car with larger engine volume uses more gas and needs more gas to accelerate per each next mile/hour.
# Thus intercept(avg gas mileage for the car) and slope(increase per each next mile/hour)  are positively correlated. 
#
# gas_mileage ~ a[model] + b[model]*speed
#
# case 2)
# Another possible example is a dependency between pollution in the region(county) and distance to the industrial area.
#  pollution ~ a[county] + b[county]*is_industrial_area
# If a county has a lot of factories than average pollution is large. But also the closer you are to the actual factory, the bigger is pollution. 
# So a[county] and b[county] are positively correlated. Still, I'm not sure that such dependency holds on practice.

## 13E3 ####
# I assume it could happen when slopes and intercepts are highly correlated across groups. The unpooled model treats correlation between slope and intercept independently for each group, while pooled relies on common distribution that can reduce the number of parameters.
# From another point of view, effective number of parameters is by construction an average variance of the likelihood across all cases in the sample, so in the pooled model this variability reduces due to shrinkage to the mean of slopes estimations.
# There is an illustration of this phenomena in the exercise 13M2 

## Medium ####
## 13M1 ####
#' Functio to generate cafe data
#'
#' @param N - number of caffes
#' @param N_visits  - number of visits per caffe
#' @param a - mean of the distribution of average wait time across all cafe
#' @param b  - mean of the distribution of the average difference between morning and afternoon wait times (slope in the model)
#' @param sigma_a - deviance of the distribution of avg. wait times
#' @param sigma_b - deviance of the distribution of of the average difference between morning and afternoon wait times
#' @param rho - correlations between parameters a and b
#'
#' @return list with keys
#'       data - dataframe with N*N_visits rows and 3 columns: cafe(int, identifier),  afternoon(boolean encoded as {0,1}), wait(real)
#'       a_cafe - vector of true a
#'       b_cafe - vector of true b
#'       Mu
#'       Sigma
generate_caffe_data <- function(N_cafes, N_visits, a, b, sigma_a, sigma_b, rho) {
  Mu <- c( a , b )
 
  # Build covariance matrix from factorised representation
  sigmas <- c(sigma_a,sigma_b) # standard deviations
  Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
  Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas) # covariance matrix
  
  set.seed(5) # used to replicate examples
  vary_effects <- mvrnorm( N_cafes , Mu , Sigma )
  a_cafe <- vary_effects[,1]
  b_cafe <- vary_effects[,2]
  
  afternoon <- rep(0:1, N_cafes*N_visits/2)
  cafe_id <- rep( 1:N_cafes , each=N_visits )
  mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
  sigma <- 0.5 # std dev within cafes
  wait <- rnorm( N_visits*N_cafes , mu , sigma )
  d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )
  list(
    data=d,
    a_cafe=a_cafe,
    b_cafe=b_cafe,
    Mu=Mu,
    Sigma=Sigma
  )
}
plot_true_cafe_data <- function(gen.cafe.data) {
  plot( gen.cafe.data$a_cafe , gen.cafe.data$b_cafe , col='red', pch=16,
        xlab="intercepts (a_cafe)" , ylab="slopes (b_cafe)", xlim=c(0,6), ylim=c(-2.5,0) )
  for ( l in c(0.1, 0.3, 0.5, 0.8, 0.99) ) {
    lines(ellipse(gen.cafe.data$Sigma, centre=gen.cafe.data$Mu, level=l), col=col.alpha("red", 0.7))
  }
}
plot_posterior_rho <- function(model, true_rho) {
  post <- extract.samples(model)
  rho_est <- mean( post$Rho[,1,2] )
  dens(post$Rho[,1,2])
  abline(v=true_rho, col='blue', lty=2)
  abline(v=rho_est, col='red', lty=2)
}
## generate data as in the chapter, but use 50 cafes instead of 20 to have better rho estimates
data.rho7 <- generate_caffe_data(50, 10, 3.5, -1, 1, 0.5, -0.7)
plot_true_cafe_data(data.rho7)

m13m1.rho7 <- map2stan(
  alist(
    wait ~ dnorm( mu , sigma ),
    mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
    c(a_cafe, b_cafe)[cafe] ~ dmvnorm2(c(a,b), sigma_cafe, Rho),
    a ~ dnorm(0,10),
    b ~ dnorm(0,10),
    sigma_cafe ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2),
    Rho ~ dlkjcorr(2)
  ) ,
  data=data.rho7$data,
  iter=5000 , warmup=2000 , chains=2, cores = 2)
plot_posterior_rho(m13m1.rho7, -0.7)

# generate data with correlation set to zero
data.rhoZero <- generate_caffe_data(50, 10, 3.5, -1, 1, 0.5, 0.0)
plot_true_cafe_data(data.rhoZero)
m13m1.rhoZero <- map2stan(m13m1.rho7, data=data.rhoZero$data, iter=5000 , warmup=2000 , chains=2, cores = 2)
plot_posterior_rho(m13m1.rhoZero, 0.0)
# Answer: posterior distribution of the correlation coefficient is centred around zero
#   It means, that even if we set up priors to encounter for correlation, but there is no correlation in the data, 
#    then the model successfully infer the absence of a relationship.
#   In other words, don't be afraid to use multivariate normal as a prior for intercepts and slopes. 
#   It doesn't hurt the inference if there is no correlation but it helps to find correlation if there is any according to data.

## 13M2 ####
m13m2.rho7 <- map2stan(
  alist(
      wait ~ dnorm(mu, sigma),
      mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
      a_cafe[cafe] ~ dnorm(a, sigma_a),
      b_cafe[cafe] ~ dnorm(b, sigma_b),
      a ~ dnorm(0,10),
      b ~ dnorm(0,10),
      sigma_a ~ dcauchy(0, 1),
      sigma_b ~ dcauchy(0, 1),
      sigma ~ dcauchy(0, 1)
  ),
  data=data.rho7$data,
  iter=5000 , warmup=2000 , chains=2, cores = 2
)
(cmp <- compare(m13m1.rho7, m13m2.rho7))
plot(cmp)
# First model (with multi-variate Gaussian priors(MVN)) has smaller WAIC and smaller number of effective parameters, thus it is better.
# Smaller number of effective parameters can be explained by the fact that MVN takes into account correlation between slope and intercept 
# and therefore requires fewer parameters to describe data. 
# The second model(m13m2.rho7) treats intercept and slope independently and uses more parameters.

extract_mean_posterior_slope_intercept <- function(model){
  post <- extract.samples(model)
  a <- apply( post$a_cafe , 2 , mean )
  b <- apply( post$b_cafe , 2 , mean )
  list(a=a, b=b)
}
post13m1 <- extract_mean_posterior_slope_intercept(m13m1.rho7)
post13m2 <- extract_mean_posterior_slope_intercept(m13m2.rho7)

# plot mean intercept(a) and slope(b) estimated by MVN model
plot( post13m1$a , post13m1$b,  xlab="intercept" , ylab="slope" ,
      pch=16 , col=rangi2, 
      xlim=c(0,6), ylim=c(-2.5,0)
)
points(post13m2$a , post13m2$b,  pch=16 , col='black') #estimated by pooled model with no correlation
for ( i in 1:length(post13m1$a) ) { 
  lines( c(post13m1$a[i], post13m2$a[i]) , c(post13m1$b[i], post13m2$b[i]) )
}

# plot contours of joint posterior distribution of slope and intercept
post <- extract.samples(m13m1.rho7)
Mu_est <- c( mean(post$a) , mean(post$b) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma_cafe[,1] )
sb_est <- mean( post$sigma_cafe[,2] )
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
# draw contours
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) {
  lines(ellipse(Sigma_est,centre=Mu_est,level=l), col=col.alpha("black",0.2))
}
# It can be seen from the plot that black points with small intercept(between 1.5 and 2.5) are shifted up 
#  in the direction of bigger slope(-0.5..0) and, on the contrary, points with large intercept(between 5 and 6) 
#  are shifted towards smaller slope (-1.5).
# This is the effect of shrinkage to the mean in the model that takes into account correlation among slope and intercept.

## 13M3 ####
data(UCBadmit)
d <- UCBadmit
d$male <- ifelse( d$applicant.gender=="male" , 1 , 0 )
d$dept_id <- coerce_index( d$dept )
d$acc_rate <- d$admit/d$applications

m13m3 <- map2stan(
  alist(
    admit ~ dbinom( applications , p ),
    logit(p) <- a_dept[dept_id] + bm_dept[dept_id]*male,
    c(a_dept, bm_dept)[dept_id] ~ dmvnorm2( c(a,bm) , sigma_dept , Rho ),
    a ~ dnorm(0,10),
    bm ~ dnorm(0,1),
    sigma_dept ~ dcauchy(0,2),
    Rho ~ dlkjcorr(2)
  ) ,
  data=d , warmup=1000 , iter=5000 , chains=4 , cores=3 )
precis(m13m3, depth=2)
plot( precis(m13m3,pars=c("a_dept","bm_dept"),depth=2) )

m13m3.nc <- map2stan(
  alist(
    admit ~ dbinom( applications , p ),
    logit(p) <- a_dept[dept_id] + bm_dept[dept_id]*male,
    c(a_dept, bm_dept)[dept_id] ~ dmvnormNC(sigma_dept , Rho ),
    #a ~ dnorm(0,10),
    sigma_dept ~ dcauchy(0,2),
    Rho ~ dlkjcorr(2)
  ) ,
  data=d , warmup=1000 , iter=5000 , chains=4 , cores=3 )
precis(m13m3.nc, depth=2)
plot( precis(m13m3.nc,pars=c("a_dept","bm_dept"),depth=2) )
(cmp <- compare(m13m3, m13m3.nc))
plot(cmp)
# models are almost identical in terms of the WAIC

plot(coeftab(m13m3, m13m3.nc))

# extract n_eff values for each model
neff_c <- precis(m13m3, 2)@output 
neff_c['param'] <- rownames(neff_c)

neff_nc <- precis(m13m3.nc, 2)@output
neff_nc['param'] <- rownames(neff_nc)

params_df <- inner_join(dplyr::select(neff_c, Mean, param, n_eff), 
                        dplyr::select(neff_nc, Mean, param, n_eff), 
                        by=c('param'),  suffix=c('_centred','_noncentred'))
# plot distributions
boxplot( list( 'm13m3'=params_df$n_eff_centred , 'm13m3.non_centered'=params_df$n_eff_noncentred ) ,
         ylab="effective samples" , xlab="model" )
# on average, the non-centred model has larger effective samples values for a_dept and bm_dept
# but smaller values for deviance estimates (sigma_dept and Rho[1,2])

t.test(params_df$n_eff_centred, params_df$n_eff_noncentred, paired = TRUE)
# but p-value of the difference is 0.2257 so from the perspective of significanse testing it is not possible to consider models as different

## 13M4 ####
data(islandsDistMatrix)
data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
d$society <- 1:10 # index observations
d$has_high_contact <- ifelse( d$contact=="high" , 1 , 0 )

# model with the Gaussian process from the chapter
m13m4.gp <- map2stan(
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + g[society] + bp*logpop,
    g[society] ~ GPL2( Dmat , etasq , rhosq , 0.01 ),
    a ~ dnorm(0,10),
    bp ~ dnorm(0,1),
    etasq ~ dcauchy(0,1),
    rhosq ~ dcauchy(0,1)
  ),
  data=list(
    total_tools=d$total_tools,
    logpop=d$logpop,
    society=d$society,
    Dmat=islandsDistMatrix),
  warmup=2000, iter=1e4 , chains=3, cores=3)
#plot(m13m4.gp)
precis(m13m4.gp, depth=2)

# one of the top model from chapter 10 with the fixed slope and intercept
m13m4.fixed <- map2stan(
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + bp*logpop + bc*has_high_contact,
    a ~ dnorm(0, 10),
    c(bp, bc) ~ dnorm(0, 1)
  ), 
  data=d, warmup=500, iter=3000 , chains=3, cores=3 
)
precis(m13m4.fixed, depth=2)

# model from the chapter 12 with varying intercept
m13m4.varintercept <- map2stan(
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + a_society[society] + bp*logpop,
    a ~ dnorm(0, 10),
    bp ~ dnorm(0, 1),
    a_society[society] ~ dnorm(0, sigma_society),
    sigma_society ~ dcauchy(0, 1)
  ),
  data=d, iter=4000, chains=3, cores=3
)
precis(m13m4.varintercept, depth=2)
(cmp <- compare(m13m4.fixed, m13m4.varintercept, m13m4.gp))

m13m4.var.intercept.slope <- map2stan(
  alist(
    total_tools ~ dpois(lambda),
    log(lambda) <- a + a_society[society] + bp_society[society]*logpop,
    a ~ dnorm(0, 10),
    c(a_society, bp_society)[society] ~ dmvnormNC(sigma_society , Rho),
    sigma_society ~ dcauchy(0, 2),
    Rho ~ dlkjcorr(2)
  ),
  data=d, warmup=2000, iter=1e4, chains=3, cores=3
)
precis(m13m4.var.intercept.slope, depth=2)
(cmp <- compare(m13m4.fixed, m13m4.varintercept, m13m4.var.intercept.slope, m13m4.gp))
# model                     WAIC pWAIC dWAIC weight    SE   dSE
# m13m4.gp                  67.2   3.9   0.0   0.63  2.13    NA
# m13m4.var.intercept.slope 69.4   4.8   2.2   0.21  1.65  0.93
# m13m4.varintercept        70.1   5.0   2.9   0.15  2.52  1.11
# m13m4.fixed               79.2   4.3  12.0   0.00 11.07 11.45
plot(cmp)
# Summary: the best model is GP model, it has the smallest WAIC and the smallest number of effective parameters
# It is exciting that effective number of parameters for m13m4.gp is 3.9 that is less than 4.3 of "all fixed" model m13m4.fixed
# The actual number of parameters that are estimated for those models are 14 and 3 correspondingly.
# Model with pooled varying intercept and slope is ranked second by WAIC and has 4.8 effective params. 
# It has a weight equal to 0.21 in model scoring.

## 13H1 ####
data("bangladesh")
d <- bangladesh
# rename variable with dots to variables with underscores, convert factor distric id into continuous integer id
d$district_id <- as.integer(as.factor(d$district))
d$use_contraception <- d$use.contraception
d$age_centered <- d$age.centered
d <- dplyr::select(d,-use.contraception, -age.centered, -district)
head(d)
# fit model from homework 12h1 to compare with
m13h1.pooled.intercept <- map2stan( 
  alist(
    use_contraception ~ dbinom(1, p),
    logit(p) <- a_district[district_id] ,
    a_district[district_id] ~ dnorm(a, sigma),
    a ~ dnorm(0, 10) ,
    sigma ~ dcauchy(0,1)
  ), 
  data=d,
  warmup=2000, iter=5000, chains=3, cores=3
)
precis(m13h1.pooled.intercept, depth=2)
# ...
# a              -0.54   0.09      -0.69      -0.40  6455    1
# sigma           0.52   0.09       0.38       0.65  2739    1
plot(m13h1.pooled.intercept)

m13h1 <- map2stan(
  alist(
    use_contraception ~ dbinom(1, p),
    logit(p) <- a_district[district_id] + b_urban_district[district_id]*urban,
    c(a_district, b_urban_district)[district_id] ~ dmvnorm2( c(a,b_urban), sigma, Rho),
    a ~ dnorm(0, 10),
    b_urban ~ dnorm(0, 1),
    sigma ~ dcauchy(0, 2),
    Rho ~ dlkjcorr(2)
  ),
  data=d,
  warmup=2000, iter=5000, chains=3, cores=3
)
(cmp <- compare(m13h1.pooled.intercept, m13h1))
plot(cmp)
# model                    WAIC pWAIC dWAIC weight    SE   dSE
# m13h1                  2468.6  52.7   0.0      1 28.13    NA
# m13h1.pooled.intercept 2514.4  35.7  45.9      0 25.07 13.88
# m13h1 definetely wins

precis(m13h1, depth=2)
# ...
#a                    -0.71   0.10      -0.87      -0.54  9000    1
#b_urban               0.70   0.17       0.43       0.98  5191    1
#sigma[1]              0.59   0.10       0.42       0.74  2527    1
#sigma[2]              0.81   0.21       0.51       1.16   892    1
# ...
#Rho[1,2]             -0.66   0.17      -0.91      -0.44  1453    1
# Rho[1,2] is a correlation coefficient between slope and intercept. Mean value of its posterior distribution is -0.66. 
# It tells us that intercept per district and a slope for the "urban" variable are negatively correlated. 
# In other words, if a probability of using contraception is high in the district then "urban" adds a small portion to an overall score. 
# And on the contrary, if the intercept is low (e.g. close to pooled mean -0.71 due to lack of data) 
# then "urban" variable has higher influence on resulting probability (bigger slope)
# Let's illustrate it with the intercept vs. slope plot

# plot 13h1-1: slope vs. intercept per district
post <- extract.samples(m13h1)
a <- apply( post$a_district , 2 , mean )
b_urban <- apply( post$b_urban_district , 2 , mean )
plot(a, b_urban, xlab="intercept per district" , ylab="slope for urban variable", pch=16 , col=rangi2)
# contour plot of estimated multivariateive distribution of intercept and slope
Mu_est <- c( mean(post$a_district) , mean(post$b_urban_district) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma[,1] )
sb_est <- mean( post$sigma[,2] )
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
# draw contours
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) {
  lines(ellipse(Sigma_est,centre=Mu_est,level=l), col=col.alpha("black",0.2))
}

# plot 13h1-2: probability of dependency between use contraception rate and urban variable
# Let's create counterfactual plots: for each district predict probability of using contraception for urban and rural women
district_ids <- min(d$district_id):max(d$district_id)
# create data.frame with data for prediction
new_data_rural <- data.frame(district_id=district_ids, urban=0)
new_data_urban <- data.frame(district_id=district_ids, urban=1)
# calculate posterior distribution of predictions
use_cont_rural <- link(m13h1, data=new_data_rural, n=3000)
use_cont_urban <- link(m13h1, data=new_data_urban, n=3000)
# use mean of distribution per example to plot
use_cont_rural_mean <- apply(use_cont_rural, 2, mean)
use_cont_urban_mean <- apply(use_cont_urban, 2, mean)

plot(use_cont_rural_mean, use_cont_urban_mean, pch=16, col='blue',
     xlab='rate in rural', 
     ylab='rate in urban',
     main='Estimated rate of using contraception in rural vs. urban area',
     xlim=c(0.1, 0.8), ylim=c(0.1, 0.8))
abline(a=0, b=1, lty=2)
# The plot shows that on average use of contraception is higher in the urban area (most of the points lie above line x=y)

# Next few lines of code add ellipses around each point with radiuses along each axis equal to the standard error of the value.
# Resulting plot is messy but gives some intuition of estimates uncertainty.
plot_ellipse <- function(xc, yc, r1, r2){
  theta <- seq(0, 2 * pi, length = 200)
  lines(x = xc + r1 * cos(theta), y = yc + r2 * sin(theta), lty=2)
}
use_cont_rural_sd <- apply(use_cont_rural, 2, sd)
use_cont_urban_sd <- apply(use_cont_urban, 2, sd)
for(i in district_ids){
  plot_ellipse(use_cont_rural_mean[i], use_cont_urban_mean[i], use_cont_rural_sd[i], use_cont_urban_sd[i])
}

# plot 13h1-3:
plot(use_cont_rural_mean, use_cont_urban_mean-use_cont_rural_mean, pch=16, col='blue',
     xlab='rate in rural', 
     ylab='difference between rate in urban and rate in rural',
     xlim=c(0.1, 0.8))
# Plot illustrates the fact that with increasing rate of contraception use in the rural are, difference between urban and rural rates declines.
# This plot is similar to the first one, but it uses natural probability scale instead of logit scale of model coefficients

## 13H2 ####
data("Oxboys")
d <- Oxboys

d %>% ggplot(aes(age, height, group=Subject, color=as.factor(Subject))) + geom_line() + geom_point(size=2) + 
      ggtitle("Boy's height as function of age", subtitle = "Each line corresponds to a particuar boy(subject). Age is centered and standardized.")

check_index(d$Subject)

d$height_normalised <- (d$height - mean(d$height))/sd(d$height)
m13h2.centered <- map2stan(
  alist(
    height_normalised ~ dnorm(mu, sigma), #predict normalised height to be able to compare intercept and slope on the same scale
    mu <- a + a_individual[Subject] + (b_age + b_age_individual[Subject])*age,
    c(a_individual, b_age_individual)[Subject] ~ dmvnorm2(0, sigma_ind, Rho),
    a ~ dnorm(0, 100),
    b_age ~ dnorm(0, 1),
    sigma_ind ~ dcauchy(0, 2),
    Rho ~ dlkjcorr(2),
    sigma ~ dcauchy(0, 2)
  ),
  data=d,
  warmup=1000, iter=3000, chains=2, cores=2
)
# Number of effective samples for a_individual and b_age_individual are really small (~200)
# I've found that moving parameter 'a' and 'b_age' into the dmvnorm2 makes inference faster and lead to much better n_eff.
# The problem with it, that is no longer possible to compare a_individual to b_age_individual to determine which one has greater influence on the result.
precis(m13h2.centered, depth=2)

post <- extract.samples(m13h2.centered)
a <- apply( post$a_individual , 2 , mean )
b_age <- apply( post$b_age_individual , 2 , mean )
plot(a, b_age, xlab="intercept per subject" , ylab="slope for age per subject", pch=16 , col=rangi2, xlim=c(-2.2,1), ylim=c(-2.2,1))
abline(a=0, b=1, lty=2)
# Roughly half of the observations lie above line x=y and another part - below. 
# From this fact it's hard to say which of varying parts (per intercept or slope) has higher influence on the height estimates.
# Deviance of the a_individual is bigger (sigma_ind[1]) than deviance of the slope b_age_individual(sigma_ind[2]). 
# I assume, that it can be interpreted as a variance of height across individuals is bigger than a variance of growth speed.

## 13H3 ####
m13h2.centered2 <- map2stan(
  alist(
    height_normalised ~ dnorm(mu, sigma), #predict normalised height to be able to compare intercept and slope on the same scale
    mu <- a_individual[Subject] + b_age_individual[Subject]*age,
    c(a_individual, b_age_individual)[Subject] ~ dmvnorm2(c(a, b_age), sigma_ind, Rho), 
    a ~ dnorm(0, 100),
    b_age ~ dnorm(0, 1),
    sigma_ind ~ dcauchy(0, 2),
    Rho ~ dlkjcorr(2),
    sigma ~ dcauchy(0, 2)
  ),
  data=d,
  warmup=1000, iter=3000, chains=2, cores=2
)
precis(m13h2.centered2, depth=2)

post <- extract.samples(m13h2.centered2)
a <- apply( post$a_individual, 2 , mean )
b_age <- apply( post$b_age_individual, 2 , mean )
plot(a, b_age, xlab="intercept per subject" , ylab="slope for age per subject", pch=16 , col=rangi2, xlim=c(-2.2,1), ylim=c(-2.2,1))
abline(a=0, b=1, lty=2)
# contour plot of estimated multivariateive distribution of intercept and slope
Mu_est <- c( mean(post$a_individual) , mean(post$b_age_individual) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma_ind[,1] )
sb_est <- mean( post$sigma_ind[,2] )
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
# draw contours
for ( l in c(0.1,0.3,0.5,0.8,0.99) ) {
  lines(ellipse(Sigma_est,centre=Mu_est,level=l), col=col.alpha("black",0.2))
}
# Model suggests that intercept and slope are correlated - for big intercepts(boys who are higher on average) the speed of growth are larger(slope for age is bigger)

## 13H4 ####
m13h4 <- map2stan(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a_individual[Subject] + b_age_individual[Subject]*age,
    c(a_individual, b_age_individual)[Subject] ~ dmvnorm2(c(a, b_age), sigma_ind, Rho), 
    a ~ dnorm(0, 100),
    b_age ~ dnorm(0, 1),
    sigma_ind ~ dcauchy(0, 2),
    Rho ~ dlkjcorr(2),
    sigma ~ dcauchy(0, 2)
  ),
  data=d,
  warmup=1000, iter=3000, chains=2, cores=2
)
precis(m13h4, depth=2)

post <- extract.samples(m13h4)
Mu_est <- c( mean(post$a_individual) , mean(post$b_age_individual) )
rho_est <- mean( post$Rho[,1,2] )
sa_est <- mean( post$sigma_ind[,1] )
sb_age_est <- mean( post$sigma_ind[,2] )
cov_ab <- sa_est*sb_age_est*rho_est
Sigma_est <- matrix( c(sa_est^2, cov_ab, cov_ab, sb_age_est^2) , ncol=2)

N_examples <- 100
new_individual_params <- mvrnorm(N_examples, Mu_est, Sigma_est)
plot(1, 1,
     xlim=c(-10, 10), #such boundaries are unrealistic but give better visualisation of the trend
     ylim=c(0, 200),  #c(h_mean-3*h_sd, h_mean+3*h_sd),
     type='n', xlab='age', ylab='height')
age_seq = seq(-10, 10, by = 0.1)
h_mean = mean(d$height)
for(i in 1:N_examples){
  intercept <- new_individual_params[i, 1]
  slope <- new_individual_params[i, 2]
  height <- intercept + age_seq*slope
  lines(age_seq, height, col=col.alpha('blue'))
}
#abline(h=h_mean, lty=2)
abline(a=Mu_est[1], b=Mu_est[2], lty=2, lwd=3)

