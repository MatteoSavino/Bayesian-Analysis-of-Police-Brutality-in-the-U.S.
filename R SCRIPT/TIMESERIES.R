#########################################
####           TIME SERIES           #### 
#########################################

{
library(rstan)
library(rethinking)
library(coda)
# for plots 
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)
# for time series
library(TTR)
library(tseries)
}

###############################

data <- read.csv('/Users/Francesco/Documents/POLICE_PROJECT/DATA/fatal_police_WP2.csv')    
data[,3] <- as.Date(data[,3], format="%Y-%m-%d")

# Data are already ordered!
#names(data)
#head(data)

# --------------------------------------------------------------------------------------

#cali <- data[which(data$state=='CA'),]
#head(cali)
#pop_cali <- data.frame(W=14206900,B=2074900,H=15261300) # 2018
victims <- data %>%
  mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
  group_by(year) %>%
  summarise(H = length(which(race=='H')),
            B = length(which(race=='B')),
            W = length(which(race=='W')))

victims %>%
  ggplot(aes(year)) + 
  geom_path(aes(y = B, group=1),col ="purple4") + 
  geom_path(aes(y = H,group=1),col ="gold") +
  geom_path(aes(y = W,group=1),col="grey3") +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

library("reshape2")
test_data_long <- melt(victims, id="year")  # convert to long format

ggplot(data=test_data_long,
       aes(x=year, y=value, colour=variable, group=2)) +
  geom_path()





victims <- data %>%
  mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
  group_by(year, month) %>%
  summarise(H = length(which(race=='H')),
            B = length(which(race=='B')),
            W = length(which(race=='W')))
victims <- victims[-71,]
matplot(victims[,c(3,4,5)],type='l',lty=1,col=c('orange','black','grey'),ylab='victims')
par(mfcol=c(1,3))
for(i in 3:5) acf(victims[,i])
for(i in 3:5) pacf(victims[,i])
adf.test(ts(victims[,5]))
dev.off()

perc_armed <- data %>%
  mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
  group_by(year, month) %>%
  summarise(armed = length(which(armed!='unarmed'))/length(armed)*100)

matplot(perc_armed[,3],type='l',lty=1,col='red',ylab='% armed')
acf(perc_armed$armed)
adf.test(perc_armed$armed)

# We create a dummy variable for the festivities:
dummy <- ifelse(victims$month=='12',1,0)





########################################################################
#################### AR1 model #########################################
########################################################################

Y <- (victims$B/victims$W)[1:65]
plot(Y,type='l', ylab='black/white')
abline(h=mean(Y),col='grey')
par(mfrow=c(1,2))
acf(Y)
pacf(Y)
adf.test(Y) ### It is stationary!!!


N <- length(Y)

data_StSp <-list(N = N, 
                 Y = Y,
                 sigma2phi = 10, 
                 a_sigma2 = 2,
                 b_sigma2 = 10,
                 sigma2m0 = 10)

inits <- function() 
{
  list(phi = 1,
       m0 = 0, 
       sigma2 = var(Y))
}

################### FIT THE MODEL
file_name = '/Users/Francesco/Documents/POLICE_PROJECT/STAN/StSp_model.stan'
StSp_model <- stan(file = file_name, 
                   data = data_StSp,
                   chains = 2, 
                   iter = 100000, 
                   warmup = 10000, 
                   thin= 10, 
                   seed = 42, 
                   init = inits,
                   algorithm = 'NUTS')

rstan::traceplot(StSp_model, pars = c("phi", "m0", "sigma2"), inc_warmup = TRUE)

# Diagnostic ------------------------------------------------------------

coda_chain <- As.mcmc.list(StSp_model, pars = c("phi", "m0", "sigma2"))
summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
#geweke.plot(coda_chain, frac1=0.1, frac2=0.5)

stan_ac(StSp_model)

# Posterior distributions -----------------------------------------------

# beta
plot_post <- StSp_model %>% 
  rstan::extract(c("phi", "m0", "sigma2")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# PLOT ------------------------------------------------------------------

params <- data.frame(rstan::extract(StSp_model, c("phi", "m0", "sigma2"), perm = T))

resultStSp <- matrix(NA, nrow = nrow(params), ncol = 5) # we simulate 5 months...
for(i in 1:5){
  if(i == 1){
    resultStSp[,i] <- apply(params, 1, function(x) rnorm(1, x[2] + x[1] * Y[length(Y)], sqrt(x[3])))
  } else {
    resultStSp[,i] <- apply(cbind(params, resultStSp[,i-1]), 1, function(x) rnorm(1, x[2] + x[1] * x[4], sqrt(x[3])))
  }
}

#y = colMeans(resultStSp)
#plot(y,type='l')

ggplot(data.frame(x = (length(Y)+1):(length(Y) + 5), 
                  y = colMeans(resultStSp),
                  ylow = apply(resultStSp, 2, quantile, p = 0.05),
                  yup = apply(resultStSp, 2, quantile, p = 0.95))) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y, x = 1:length(Y)), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +           # how many uncertainties!!!
  geom_vline(aes(xintercept = length(Y)), col = 2, lty = 2) +
  ggtitle("AR1 model - COVID data Lombardia")

# a disaster!!!

################################### -------------------------------------------------------------
############# ARX(1) ##############
###################################

Y <- victims$B/victims$W
Y <- Y[1:(length(Y)-5)]
X <- perc_armed$armed
N <- length(Y)

data_DRM <-list(N = N, 
                Y = Y,
                X = X[1:length(Y)], ### we ignore the last 5 days, which will be used for the prediction
                sigma_coef = rep(10, 3), 
                sigma2X_par = c(2, 10), 
                sigma2Y_par = c(2, 10))

inits <- function() 
{
  list(sigma2Y = var(Y),
       m0Y = 0,
       beta = 0, 
       gamma = 0)
}

# FIT THE MODEL

file_name = '/Users/Francesco/Documents/POLICE_PROJECT/STAN/DLR_model.stan'
DRM_model <- stan(file = file_name, 
                  data = data_DRM,
                  chains = 2, 
                  iter = 100000, 
                  warmup = 10000, 
                  thin= 10, 
                  seed = 42, 
                  init = inits,
                  algorithm = 'NUTS')
rstan::traceplot(DRM_model, pars = c("beta", "gamma"), inc_warmup = TRUE)
rstan::traceplot(DRM_model, pars = c("m0Y", "sigma2Y"), inc_warmup = TRUE)

# Diagnostic 

coda_chain <- As.mcmc.list(DRM_model, pars = c("beta", "gamma"))
summary(coda_chain)

# Compute the gelman and rubin's convergence diagnostic
# ‘potential scale reduction factor’
# And quoting Gelman:
# “Inference is normal science. Model-checking is revolutionary science.”
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
#geweke.plot(coda_chain, frac1=0.1, frac2=0.5)

# autocorrelation and plot
stan_ac(DRM_model)

# Posterior distributions 

# coefs
plot_post <- DRM_model %>% 
  rstan::extract(c("beta", "gamma")) %>% 
  as_tibble() %>% 
  map_df(as_tibble, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# intercepts
plot_post <- DRM_model %>% 
  rstan::extract(c("m0Y", "sigma2Y")) %>% 
  as_tibble() %>% 
  map_df(as_tibble, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# PREDICTION PLOT 
params <- data.frame(rstan::extract(DRM_model, c("beta", "m0Y", "gamma", "sigma2Y"), perm = T))
Xtemp <- X[(length(Y) + 1):(length(Y) + 5)]  # we consider the months we have initially discarded (66:70)
resultDLR <- matrix(NA, nrow = nrow(params), ncol = 5 )
for(i in 1:5){
  if(i == 1){
    resultDLR[,i] <- apply(params, 1, function(x) rnorm(1, x[2] + x[1] * Xtemp[i] + x[3] * Y[length(Y)], sqrt(x[4])))
  } else {
    resultDLR[,i] <- apply(cbind(params, resultDLR[,i-1]), 1, function(x) rnorm(1,  x[2] + x[1] * Xtemp[i] + x[3] * x[5], sqrt(x[4])))
  }
}
#resultDLR <- resultDLR[,1:5]
#resultStSp <- resultStSp[,1:5]
Y_true <- (victims$B/victims$W)[66:70]

ggplot(data.frame(x = (length(Y)+1):(length(Y) + 5), 
                  y = colMeans(resultDLR),
                  ylow = apply(resultDLR, 2, quantile, p = 0.05),
                  yup = apply(resultDLR, 2, quantile, p = 0.95))) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5),     # AR(1) prediction
                              y = colMeans(resultStSp),
                              ylow = apply(resultStSp, 2, quantile, p = 0.05),
                              yup = apply(resultStSp, 2, quantile, p = 0.95)), 
            aes(x = x, y = y), lty = 2, col = 3) + 
  geom_ribbon(data = data.frame(x = (length(Y)+1):(length(Y) + 5), 
                                y = colMeans(resultStSp),
                                ylow = apply(resultStSp, 2, quantile, p = 0.05),
                                yup = apply(resultStSp, 2, quantile, p = 0.95)),
              aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2, fill = 3) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y, x = 1:length(Y)), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = length(Y)), col = 2, lty = 2) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5),     # True values
                              y = Y_true),
              aes(x = x, y = y), lty = 2, col = 2) +
  ggtitle("AR1 model & ARX1 - USA")


##################################### -------------------------------------------------------------
############# ARX(2,1) ##############
#####################################

Y <- victims$B/victims$W
Y <- Y[1:(length(Y)-5)]
X <- perc_armed$armed
Z <- dummy[1:65]

N <- length(Y)

data_DRM <-list(N = N, 
                Y = Y,
                X = X[1:length(Y)],
                Z = Z,
                sigma_coef = rep(10, 5), 
                sigma2X_par = c(2, 10), 
                sigma2Y_par = c(2, 10))

inits <- function() 
{
  list(sigma2Y = var(Y),
       m0Y = 0,
       beta1 = 0,
       beta2 = 0,
       gamma1 = 0,
       gamma2 = 0
  )
}

# FIT THE MODEL

file_name = '/Users/Francesco/Documents/POLICE_PROJECT/STAN/ARX2_model.stan'
ARX_model <- stan(file = file_name, 
                  data = data_DRM,
                  chains = 2, 
                  iter = 200000, 
                  warmup = 10000, 
                  thin= 10, 
                  seed = 42, 
                  init = inits,
                  algorithm = 'NUTS')

rstan::traceplot(ARX_model, pars = c("beta1", "beta2","gamma1","gamma2"), inc_warmup = TRUE)
rstan::traceplot(ARX_model, pars = c("m0Y", "sigma2Y"), inc_warmup = TRUE)

coda_chain <- As.mcmc.list(ARX_model, pars = c("beta1", "beta2","gamma1","gamma2"))
summary(coda_chain)

# Gelman's convergence diagnostic
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
#geweke.plot(coda_chain, frac1=0.1, frac2=0.5)

# autocorrelation and plot
stan_ac(ARX_model)

# coefs
plot_post <- ARX_model %>% 
  rstan::extract(c("beta1", "beta2")) %>% 
  as_tibble() %>% 
  map_df(as_tibble, .id = 'param')
plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# coefs
plot_post <- ARX_model %>% 
  rstan::extract(c("gamma1", "gamma2")) %>% 
  as_tibble() %>% 
  map_df(as_tibble, .id = 'param')
plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# PREDICTION PLOT ------------------------
params <- data.frame(rstan::extract(ARX_model, c("m0Y","beta1","beta2","gamma1","gamma2", "sigma2Y"), perm = T))
Xtemp <- X[(length(Y) + 1):(length(Y) + 5)] # the values we discarded (66:70)
Ztemp <- dummy[(length(Y) + 1):(length(Y) + 5)]
resultARX <- matrix(NA, nrow = nrow(params), ncol = 5)
for(i in 1:5){
  if(i == 1){
    resultARX[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[2]*Xtemp[i] + x[3]*Z[i] + x[4]*Y[length(Y)] + x[5]*Y[length(Y)-1] , sqrt(x[6])))
  } 
  if(i==2){
    resultARX[,i] <- apply(cbind(params, resultARX[,i-1]), 1, function(x) rnorm(1,  x[1] + x[2]*Xtemp[i] + x[3]*Z[i] + x[4]*x[7] + x[5]*Y[length(Y)], sqrt(x[6])))
  }
  if(i>2){
    resultARX[,i] <- apply(cbind(params,resultARX[,i-1],resultARX[,i-2]), 1, function(x) rnorm(1,  x[1] + x[2]*Xtemp[i] + x[3]*Z[i] + x[4]*x[8] + x[5]*x[7], sqrt(x[6])))
  }
}

ggplot(data.frame(x = (length(Y)+1):(length(Y) + 5), 
                  y = colMeans(resultARX),
                  ylow = apply(resultARX, 2, quantile, p = 0.05),
                  yup = apply(resultARX, 2, quantile, p = 0.95))) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5), 
                              y = colMeans(resultStSp),
                              ylow = apply(resultStSp, 2, quantile, p = 0.05),
                              yup = apply(resultStSp, 2, quantile, p = 0.95)), 
            aes(x = x, y = y), lty = 2, col = 3) + 
  geom_ribbon(data = data.frame(x = (length(Y)+1):(length(Y) + 5), 
                                y = colMeans(resultStSp),
                                ylow = apply(resultStSp, 2, quantile, p = 0.05),
                                yup = apply(resultStSp, 2, quantile, p = 0.95)),
              aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2, fill = 3) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = Y, x = 1:length(Y)), aes(x = x, y = Y)) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2) +
  geom_vline(aes(xintercept = length(Y)), col = 2, lty = 2) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5),     # True values
                              y = Y_true),
            aes(x = x, y = y), lty = 2, col = 2) +
  ggtitle("ARX model & ARX(2,1) - USA")





#########################
###      ARMA         ###
#########################

Y <- (victims$B/victims$W)[1:65]
plot(Y,type='l', ylab='black/white')
abline(h=mean(Y),col='grey')

N <- length(Y)

data_ARMA <-list(N = N, 
                 Y = Y,
                 #sigma2phi = 10, 
                 a_sigma2 = 2,
                 b_sigma2 = 10,
                 sigma2m0 = 10)

inits <- function() 
{
  list(#phi = 1,
       theta = 1,
       m0 = 0, 
       sigma2 = var(Y))
}

################### FIT THE MODEL
file_name = '/Users/Francesco/Documents/POLICE_PROJECT/STAN/ARMA_model.stan'
ARMA_model <- stan(file = file_name, 
                   data = data_ARMA,
                   chains = 2, 
                   iter = 100000, 
                   warmup = 10000, 
                   thin= 10, 
                   seed = 42, 
                   init = inits,
                   algorithm = 'NUTS')

rstan::traceplot(ARMA_model, pars = c("phi","theta", "m0", "sigma2"), inc_warmup = TRUE)
rstan::traceplot(ARMA_model, pars = c("theta", "m0", "sigma2"), inc_warmup = TRUE)

# Diagnostic 

coda_chain <- As.mcmc.list(ARMA_model, pars = c("phi","theta", "m0", "sigma2"))
coda_chain <- As.mcmc.list(ARMA_model, pars = c("theta", "m0", "sigma2"))

summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
#geweke.plot(coda_chain, frac1=0.1, frac2=0.5)

stan_ac(ARMA_model)

# Posterior distributions 

# beta
plot_post <- ARMA_model %>% 
  rstan::extract(c("phi","theta", "m0", "sigma2")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# PLOT 
# TODO
params <- data.frame(rstan::extract(ARMA_model, c("m0","phi","theta","sigma2","epsilon"), perm = T))
params <- cbind(params,as.data.frame(matrix(NA, nrow = nrow(params), ncol = 4))) # We introduce 4 more columns for the new epsilon
resultARMA <- matrix(NA, nrow = nrow(params), ncol = 5) # we simulate 5 months...
for(i in 1:5){
  if(i == 1){
    # params[69] = epsilon(65) !!!
    resultARMA[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[2]*Y[length(Y)] + x[3]*x[69], sqrt(x[4])))
    params[69+i] <-  resultARMA[,i] - (params[1]+params[2]*Y[length(Y)]+params[3]*params[69])
  } else {
    resultARMA[,i] <- apply(cbind(params, resultARMA[,i-1]), 1, function(x) rnorm(1, x[1] + x[2]*x[74] + x[3]*x[69+i-1], sqrt(x[4])))
    params[69+i] <-  resultARMA[,i] - (params[1]+params[2]*resultARMA[,i-1]+params[3]*params[69+i-1])
  }
}

params <- data.frame(rstan::extract(ARMA_model, c("m0","theta","sigma2","epsilon"), perm = T))
params <- cbind(params,as.data.frame(matrix(NA, nrow = nrow(params), ncol = 4))) # We introduce 4 more columns for the new epsilon
resultARMA <- matrix(NA, nrow = nrow(params), ncol = 5) # we simulate 5 months...
for(i in 1:5){
  if(i == 1){
    # params[69] = epsilon(65) !!!
    resultARMA[,i] <- apply(params, 1, function(x) rnorm(1, x[1] + x[2]*x[68], sqrt(x[3])))
    params[68+i] <-  resultARMA[,i] - (params[1]+params[2]*params[68])
  } else {
    resultARMA[,i] <- apply(cbind(params, resultARMA[,i-1]), 1, function(x) rnorm(1, x[1] + x[2]*x[68+i-1], sqrt(x[3])))
    params[68+i] <-  resultARMA[,i] - (params[1]+params[2]*params[68+i-1])
  }
}
Y_true <- (victims$B/victims$W)[66:70]
ggplot(data.frame(x = (length(Y)+1):(length(Y) + 5), 
                  y = colMeans(resultARMA),
                  ylow = apply(resultARMA, 2, quantile, p = 0.05),
                  yup = apply(resultARMA, 2, quantile, p = 0.95))) +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2,fill=3) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5), 
                              y = colMeans(resultStSp),
                              ylow = apply(resultStSp, 2, quantile, p = 0.05),
                              yup = apply(resultStSp, 2, quantile, p = 0.95)), 
            aes(x = x, y = y), lty = 1, col = 3) + 
  geom_ribbon(data = data.frame(x = (length(Y)+1):(length(Y) + 5), 
                                y = colMeans(resultStSp),
                                ylow = apply(resultStSp, 2, quantile, p = 0.05),
                                yup = apply(resultStSp, 2, quantile, p = 0.95)),
              aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2, fill = 1) +
  geom_line(aes(x = x, y = y), lty = 2 ) + theme_bw() +
  geom_line(data = data.frame(Y = Y, x = 1:length(Y)), aes(x = x, y = Y)) +
  #geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.2,fill=3) +
  geom_vline(aes(xintercept = length(Y)), col = 2, lty = 2) +
  geom_line(data = data.frame(x = (length(Y)+1):(length(Y) + 5),     # True values
                              y = Y_true),
            aes(x = x, y = y), lty = 2, col = 2) +
  ggtitle("AR(1) VS MA(1) - USA")
