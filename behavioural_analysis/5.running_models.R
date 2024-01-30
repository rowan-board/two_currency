# running directed action


# load the packages
library(rstan)
library(tidyverse)
library(R.matlab)
library(bayesplot)
library(wesanderson)
library(waterfalls)

# load the data
load('modelling/datalists/online_datalist.RData')

# running stan
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

modelFile <- 'modelling/models/full_model.stan'

nIter     <- 2000
nChains   <- 4 
nWarmup   <- floor(nIter/2)
nThin     <- 1

fit_rl <- stan(modelFile, 
               data    = datalist, 
               chains  = nChains,
               iter    = nIter,
               warmup  = nWarmup,
               thin    = nThin,
               init    = "random",
               seed    = 3210
)


# check values rhats etc
print(fit_rl)

rl_sum <- summary(fit_rl)
rl_sum <- rl_sum$summary
write.csv(rl_sum, 'modelling/summaries/full_model_online.csv')
rl_ppc <- summary(fit_rl, pars = c('pred'))
rl_ppc <- rl_ppc$summary
write.csv(rl_ppc, 'modelling/summaries/full_model_online_ppc.csv')
rl_q <- summary(fit_rl, pars = c('ar', 'ag', 'br','bg'))
rl_q <- rl_q$summary
write.csv(rl_q, 'modelling/summaries/q.csv')
rl_chosen_net <- summary(fit_rl, pars = c('chosen_net'))
rl_chosen_net <- rl_chosen_net$summary
write.csv(rl_chosen_net, 'modelling/summaries/full_model_online_chosen_net.csv')
rl_chosen_devalued <- summary(fit_rl, pars = c('chosen_devalued'))
rl_chosen_devalued <- rl_chosen_devalued$summary
write.csv(rl_chosen_devalued, 'modelling/summaries/chosen_devalued.csv')
rl_colour_diff <- summary(fit_rl, pars = c('chosen_colour_diff'))
rl_colour_diff <- rl_colour_diff$summary
write.csv(rl_colour_diff, 'modelling/summaries/colour_diff.csv')
rl_chosen_diff <- summary(fit_rl, pars = c('chosen_diff'))
rl_chosen_diff <- rl_chosen_diff$summary
write.csv(rl_chosen_diff, 'modelling/summaries/full_model_online_chosen_diff.csv')
rl_unchosen_diff <- summary(fit_rl, pars = c('unchosen_diff'))
rl_unchosen_diff <- rl_unchosen_diff$summary
write.csv(rl_unchosen_diff, 'modelling/summaries/unchosen_diff.csv')

save(fit_rl, file='modelling/fitted_models/full_model.RData')

# checking diagnostics with plots
plot_dens_tau <- stan_plot(fit_rl, pars=c('mu_tau','sigma_tau','tau'), show_density=T, fill_color = 'skyblue')
plot_dens_lr <- stan_plot(fit_rl, pars=c('a_lr', 'b_lr', 'lr'), show_density=T, fill_color = 'skyblue')
plot_trace <- stan_trace(fit_rl_4, pars=c('k_tau','theta_tau', 'a_lr', 'b_lr'), inc_warmup = F)


### model comparison ### 

model_1 <- load('modelling/fitted_models/model_1.RData')
model_2 <- load('modelling/fitted_models/model_2.RData')
model_3 <- load('modelling/fitted_models/model_3.RData')
model_4 <- load('modelling/fitted_models/model_4.RData')
model_5 <- load('modelling/fitted_models/model_5.Rdata')
model_6 <- load('modelling/fitted_models/model_6.RData')
model_7 <- load('modelling/fitted_models/model_7.RData')
model_8 <- load('modelling/fitted_models/model_9.RData')
model_6 <- load('modelling/fitted_models/model_9.RData')
model_10 <- load('modelling/fitted_models/model_10.RData')

# get the log likelihoods
LL1 <- loo::extract_log_lik(fit_rl_1)
LL2 <- loo::extract_log_lik(fit_rl_2)
LL3 <- loo::extract_log_lik(fit_rl_3)
LL4 <- loo::extract_log_lik(fit_rl_4)
LL5 <- loo::extract_log_lik(fit_rl_5)
LL6 <- loo::extract_log_lik(fit_rl_9)
LL7 <- loo::extract_log_lik(fit_rl_7)
LL8 <- loo::extract_log_lik(fit_rl_9)
LL8 <- loo::extract_log_lik(fit_rl_8)
LL8 <- loo::extract_log_lik(fit_rl_10)

LL1_mean <- mean(LL1)
LL2_mean <- mean(LL2)
LL3_mean <- mean(LL3)
LL4_mean <- mean(LL4)
LL5_mean <- mean(LL5)
LL6_mean <- mean(LL6)
LL7_mean <- mean(LL7)
LL8_mean <- mean(LL8)
LL9_mean <- mean(LL9)

# calculate and compare looic
loo1 <- loo::loo(LL1)
loo2 <- loo::loo(LL2)
loo3 <- loo::loo(LL3)
loo4 <- loo::loo(LL4)
loo5 <- loo::loo(LL5)
loo6 <- loo::loo(LL6)
loo7 <- loo::loo(LL7)
loo8 <- loo::loo(LL8)
loo9 <- loo::loo(LL9)

loo::loo_compare(loo1,loo2,loo3,loo4,loo5,loo6)

# plot the looic
looic = rbind(loo1[[7]]-loo5[[7]], loo2[[7]]-loo5[[7]], loo3[[7]]-loo5[[7]], loo4[[7]]-loo5[[7]], loo5[[7]]-loo5[[7]], loo6[[7]]-loo5[[7]])
looic_se = rbind(loo1[[10]], loo2[[10]], loo3[[10]], loo4[[10]],loo5[[10]], loo6[[10]])
looic_df <- cbind(looic, looic_se) %>%
  as.data.frame()
looic_df$model <- c("Null_model","simple_RW","two_channel","two_channel_motivation","motivation_counterfactual","free_motivation")



ggplot(looic_df, aes(fct_rev(fct_reorder(model, V1)),
                     V1)) +
  geom_bar(stat = "identity", fill='#FF9B5E') +
  theme_classic()+
  theme(axis.text=element_text(size=12),
           axis.title=element_text(size=16),
           plot.title=element_text(size=18)) +
  labs(y ="LOOIC difference",x='Model') +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  ggtitle('Two currency task models')

## participant wise looics
LL5_ppt <- loo5$pointwise[,4]
x=as.character(1:length(LL5_ppt))
df <- as.data.frame(cbind(x,LL5_ppt))

# subtract min from all values
min_looic <- as.numeric(min(df$LL5_ppt))
df$LL5_ppt <- as.numeric(df$LL5_ppt)-min_looic


ggplot(df, aes(fct_rev(fct_reorder(x, LL5_ppt)),
                     LL5_ppt)) +
  geom_bar(stat = "identity", fill='#FF9B5E') +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title=element_text(size=18)) +
  labs(y ="LOOIC difference",x='subject') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle('Individual fits')




# manual calculations for BIC

nt = 500 # number of trials

BIC_1 = sum(1 * log(nt) - 2*colMeans(LL1))
BIC_2 = sum(2 * log(nt) - 2*colMeans(LL2))
BIC_3 = sum(2 * log(nt) - 2*colMeans(LL3))
BIC_4 = sum(2 * log(nt) - 2*colMeans(LL4))
BIC_5 = sum(2 * log(nt) - 2*colMeans(LL5))
BIC_6 = sum(3 * log(nt) - 2*colMeans(LL6))
BIC_7 = sum(3 * log(nt) - 2*colMeans(LL7))
BIC_8 = sum(4 * log(nt) - 2*colMeans(LL8))
