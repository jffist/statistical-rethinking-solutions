rm(list=ls())
library(ggplot2)
library(dplyr)

## Medium ####
# 2M1, 2M2
data <- c('W','W', 'W')
data <- c('W','W', 'W', 'L')
data <- c('L','W', 'W', 'L', 'W', 'W', 'W')

n = length(data)
w = sum(data=='W')

p_grid <- seq(0, 1, length.out = 25)
ng = length(p_grid)

#prior <- rep(1, ng)
prior <- ifelse(p_grid<0.5, 0, 1)

likelihood <- dbinom(w, n, p_grid)

posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

df <- data.frame(param=c(p_grid,p_grid),
                 prob=c(prior,posterior),
                 ptype=c(rep('prior',ng), rep('posterior',ng)))
x <- ggplot(df %>% filter(ptype=='posterior'), 
            aes(x=param, y=prob, group=ptype, color=ptype)) +
  geom_line() + geom_point()
print(x)


# 2M3
# p(Earth|Land) = P(Land|Earth)*P(Earth)/ (P(Land|Earth)*P(Earth)+P(Land|Mars)*P(Mars))
p_e_l = (1-0.7)*0.5 / ((1-0.7)*0.5 + 1.*0.5)
print(p_e_l)

# 2M4
# BB WB WW
# By probability formula
# P(Down=B|Upper=B) = P(Down=B,Upper=B) / P(Upper=B) = P(BB) / (P(BB) + P(WB)*P(Upper=B|WB)) =  
#    = (1/3) / (1/3 + 1/3 * 1/2) = 2/3
# by counting
card.bb.likelihood <- 2
card.wb.likelihood <- 1
card.ww.likelihood <- 0
likelihood <- c(card.bb.likelihood, card.wb.likelihood, card.ww.likelihood)
prior <- c(1,1,1)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)
posterior[1] == 2/3

# 2M5
# by counting
card.bb.likelihood <- 2
card.wb.likelihood <- 1
card.ww.likelihood <- 0
card.bb.likelihood <- 2

likelihood <- c(card.bb.likelihood, card.wb.likelihood, card.ww.likelihood, card.bb.likelihood)
prior <- c(1,1,1,1)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)

# the result is sum of probabilities to draw card 1 or 4
posterior[1]+posterior[4]

#by formula
# P(Down=B|Upper=B) = P(Down=B,Upper=B) / P(Upper=B) = 
#    = P(BB) / (P(BB) + P(WB)*P(Upper=B|WB)) =  
#    = ( P(BB_1)+P(BB_4) ) / (P(BB_1) + P(BB_4) + P(WB)*P(Upper=B|WB))
#    =  (2/4) / (2/4 + 1/4 * 1/2) = 4/5 = 0.8

# 2M6
# BB(1) WB(2) WW(3) --> priors (1/6, 2/6, 3/6)
# By probability formula
# P(Down=B|Upper=B) = P(Down=B,Upper=B) / P(Upper=B) = P(BB) / (P(BB) + P(WB)*P(Upper=B|WB)) =  
#    = (1/6) / (1/6 + 2/6 * 1/2) = 1/(1+1) = 0.5
# by counting
card.bb.likelihood <- 2
card.wb.likelihood <- 1
card.ww.likelihood <- 0
likelihood <- c(card.bb.likelihood, card.wb.likelihood, card.ww.likelihood)
prior <- c(1,2,3)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)
posterior[1] == 0.5

# 2M7
card.bb.likelihood <- 2*3 #2 blacks choices, for each there are 3 option = ww->2 + wb->1
card.wb.likelihood <- 1*2 
card.ww.likelihood <- 0
likelihood <- c(card.bb.likelihood, card.wb.likelihood, card.ww.likelihood)
prior <- c(1,1,1)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)
posterior[1] == 0.75

## Hard ####
# 2H1
# p(twins) = p(species=A)*p(twins|A)+p(species=B)*p(twins|B)
# before first breeding P(species=A)=p(species=B)=0.5
# to get answer for the second round we need to calculate new p(A) and p(B)
# it's posterior of first breeding that becomes prior for the second one
# p(tweens_2|tweens_1) = p(species=A|tweens_1)*p(twins|A)+p(species=B|tweens_1)*p(twins|B)
p_twins_A <- 0.1
p_twins_B <- 0.2
likelihood <- c(p_twins_A, p_twins_B)
prior <- c(1, 1)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)

# result
sum(posterior*likelihood)

# 2H2
# vanila Bayes rule
p_twins_A <- 0.1
p_twins_B <- 0.2
likelihood <- c(p_twins_A, p_twins_B)
prior <- c(1, 1)
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)

posterior[1] #0.33

# 2H3
p_twins_A <- 0.1
p_twins_B <- 0.2
# first Bayesian update
likelihood_twins <- c(p_twins_A, p_twins_B)
prior <- c(1, 1)
posterior <- prior * likelihood_twins
posterior <- posterior/sum(posterior)

# second Bayesian update
likelihood_single <- c(1-p_twins_A, 1-p_twins_B)
prior <- posterior
posterior <- prior * likelihood_single
posterior <- posterior/sum(posterior)

posterior[1] #0.36

# alternatively it could be implemented with single update
p_twins_A <- 0.1
p_twins_B <- 0.2
# likelihood of two events (p(twins_step1 & single_step2|species=X))
likelihood_twins_single <- c(p_twins_A*(1-p_twins_A), 
                             p_twins_B*(1-p_twins_B))
prior <- c(1, 1)
posterior <- prior * likelihood_twins_single
posterior <- posterior/sum(posterior)

posterior[1] #0.36

# 2H4
#p1
likelihood_test <- c(0.8, 1-0.65)
prior <- c(1, 1)
posterior_vet <- prior * likelihood_test
posterior_vet <- posterior_vet/sum(posterior_vet)

posterior_vet[1] #0.6956522
#p2
p_twins_A <- 0.1
p_twins_B <- 0.2
likelihood_twins <- c(p_twins_A, p_twins_B)
prior <- posterior_vet
posterior <- prior * likelihood_twins
posterior <- posterior/sum(posterior)

posterior[1] #0.5333333
