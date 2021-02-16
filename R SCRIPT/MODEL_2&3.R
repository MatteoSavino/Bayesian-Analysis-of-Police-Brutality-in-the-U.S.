{
library(rstan) 
library(coda)
library(rethinking)
# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
require(gplots)
require(ggpubr)
}

my_WAIC <- function(fit, param){
  llik   <- rstan::extract(fit, param)[[1]]
  p_WAIC <- sum(apply(llik, 2, var))
  lppd   <- sum(apply(llik, 2, function(x) log(mean(exp(x)))))
  WAIC   <- - 2 * lppd + 2 * p_WAIC
  return(WAIC)
}

# ------------------------------------------------------------------------#
############################### MODEL 2 ###################################
# ------------------------------------------------------------------------#

setwd('/Users/Francesco/Documents/POLICE_PROJECT')
data <- read.csv('DATA/fatal_police_WP2.csv')      # fatal_police_kaggle/PoliceKillingsUS.csv


#                  OR
#        fatal_police_WP.csv        
# Data are already ordered!
names(data)
head(data)

####################################
####       STAN PAPER           ####
####################################

### We erase incomplete cases:
d <- data[which((data[,'armed'])!=''),]

Results<-c()
Results <- ifelse(d$armed !="unarmed" & d$race=="B", "BlackArmed",NA)
Results <- ifelse(d$armed =="unarmed" & d$race=="B", "BlackUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="H", "HispanicArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="H", "HispanicUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="W", "WhiteArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="W", "WhiteUnarmed",Results)

table <- table(d$state,Results)
dim(table)
table <- table[which(complete.cases(table[,'BlackArmed'])),] 
dim(table)
head(table)


# Extract population size data for each state (REM: Data from 2018)
usa <- read.csv('DATA/usa_population.csv',header=T) # open USA_populaton.csv 
#load('geodata.RData') 
#dim(usa)
#dim(k2) 
#usa <- usa[-1,]
#data <- data[-1,]
#table <- table[-1,]

Nwhite <- usa$White
Nblack <- usa$Black
Nhispanic <- usa$Hispanic
N_race <- c(Nblack,Nblack,Nhispanic,Nhispanic,Nwhite,Nwhite) # !!!

# Extract shooting data for each state
UnarmedBlack <- as.numeric(table[,'BlackUnarmed'])
ArmedBlack <- as.numeric(table[,'BlackArmed'])
UnarmedHispanic <- as.numeric(table[,'HispanicUnarmed'])
ArmedHispanic <- as.numeric(table[,'HispanicArmed'])
UnarmedWhite <- as.numeric(table[,'WhiteUnarmed'])
ArmedWhite <- as.numeric(table[,'WhiteArmed'])
ArmedWhite[which(ArmedWhite==0)] <- 1
UnarmedWhite[which(UnarmedWhite==0)] <- 1
Data <- c(UnarmedBlack,ArmedBlack,UnarmedHispanic,ArmedHispanic,UnarmedWhite,ArmedWhite) # !!!

offset <- log(c(UnarmedWhite/Nwhite*Nblack, ArmedWhite/Nwhite*Nblack,
              UnarmedWhite/Nwhite*Nhispanic, ArmedWhite/Nwhite*Nhispanic,
             UnarmedWhite/Nwhite*Nwhite, ArmedWhite/Nwhite*Nwhite))



N <- length(Data) # !!!
# We create a dummy variable W:  (we consider all the 6 groups)
W <- cbind(rep(c(1,0),c(51,N-51)),
           rep(c(0,1,0),c(51,51,N-102)),
           rep(c(0,1,0),c(102,51,N-153)),
           rep(c(0,1,0),c(153,51,N-204)),
           rep(c(0,1,0),c(204,51,N-255)),
           rep(c(0,1),c(N-51,51)))

# This variable will be useful for u[i] !!!
#group <- rep(1:6, each = 51)
group <- rep(1:51, 6)

#------------------------------- We upload the Geo Data:
#load('DATA/DATA.RData')
#load('CODICE MATTEO/dati_usa.RData')

#X <- X[,-c(1,5,6,7)]
#X[,c(2,3)] <- X[,c(2,3)]/100
#X <- cbind(X,k2[,c(16,17,18)]/100,k2[,14],k2[,19],k2[,c(22,23,24)])
#DATA_MAX <- X
#save(DATA_MAX,file='DATA/DATA_MAX.RData')
load('DATA/DATA_MAX.RData')
head(DATA_MAX)
DATA_MAX <- DATA_MAX[,-c(4,5,6)]
B_W <- data.frame(B_W=Nblack/Nwhite*100)
H_W <- data.frame(H_W=Nhispanic/Nwhite*100)
DATA_MAX <- cbind(DATA_MAX[,1:3],B_W,H_W,DATA_MAX[,4:8])
DATA_MAX$income <- scale(DATA_MAX$income)
DATA_MAX$`Median Age` <- scale(DATA_MAX$`Median Age`)
DATA_MAX$Violent.crime <- scale(DATA_MAX$Violent.crime)
DATA_MAX <- DATA_MAX[,-10]
head(DATA_MAX)

# Which variables do we want???
X <- DATA_MAX[,1:6]
head(X)
# We compute the covariance matrix:
X <- as.matrix(X)
#B0 <- solve(t(X)%*%X)
# Now we multiply the X matrix in order to fit all the data:
X <- rbind(X,X,X,X,X,X)
dim(X)
p_fix <- ncol(X)
p_ran <- ncol(W)



model_dat  <- list(N=N,
                   p_fix=p_fix,
                   p_ran=p_ran,
                   group=group,
                   X=X,         
                   W=W,
                   #--------- normal-inv-gamma model parameter
                   mu_beta=0,
                   a_sigma2_beta=2.,
                   b_sigma2_beta=10.,
                   #--------- random effect parameters:
                   mu_gamma=-5,
                   a_sigma2_gam=2.,#0.01,
                   b_sigma2_gam=10.,#0.01,
                   #N_race=N_race,
                   Y=Data,
                   offset=offset
                   )   

inits <- function() 
{
  list( beta = rep(0.1,p_fix),
        sigma2_beta = rep(10,p_fix),
        gamma = rep(0.1,p_ran),
        sigma2_gamma = rep(10,p_ran)
        #u = rep(0.1,p_ran)
        )
}

### Here we select one of our STAN files
file_stan <- '/Users/Francesco/Documents/POLICE_PROJECT/STAN_MODELS/POISSON.stan'

fitKilling <- stan(file=file_stan,
                   data = model_dat,
                   init=inits,
                   seed=2613,
                   thin=1,
                   iter = 30000,
                   warmup=5000,
                   chains=1,
                   refresh=100,
                   algorithm='NUTS',
                   control=list(max_treedepth=15))

#save(fitKilling,file='STAN_MODELS/poisson_model_true.RData')
load('STAN_MODELS/poisson_model_true.RData')

# We immediately test our model
my_WAIC(fitKilling, "log_lik") 


# R = exp { gamma + beta * x }

rstan::traceplot(fitKilling, pars = "sigma2_beta", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "sigma2_gamma", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "beta", inc_warmup = FALSE)
rstan::traceplot(fitKilling, pars = "gamma", inc_warmup = FALSE)


# Diagnostic ------------------------------------------------------------

coda_chain <- As.mcmc.list(fitKilling, pars = c('beta','sigma2_beta','gamma','sigma2_gamma'))
summary(coda_chain)
# Compute the gelman and rubin's convergence diagnostic
# ‘potential scale reduction factor’
# And quoting Gelman:
# “Inference is normal science. Model-checking is revolutionary science.”
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain, frac1=0.1, frac2=0.6)

# autocorrelation and plot
acfplot(coda_chain, lag.max = 100)

# Posterior distributions -----------------------------------------------

plot_post <- fitKilling %>% 
  rstan::extract("beta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none") +
  geom_vline(aes(xintercept=0, col='red'))

# gamma
plot_post <- fitKilling %>% 
  rstan::extract("gamma") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none") +
  geom_vline(aes(xintercept=0, col='red'))

# sigma
plot_post <- fitKilling %>% 
  rstan::extract("sigma2_gamma") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")


#### PLOT #### ----------------------------------------------------

# We extract the values of mu:
mu <- fitKilling %>% 
         rstan::extract("mu")
mu <- as.data.frame(mu)


#mu_BU <- mu[,1:51] 
#mu_BA <- mu[,(1+51):(51+51)]
#mu_HU <- mu[,(1+51+51):(51+51+51)]
#mu_HA <- mu[,(1+51+51+51):(51+51+51+51)]


quantile(mu[,1],0.95)
q0.05 <- apply(mu, 2, quantile, p=0.05)
q0.95 <- apply(mu, 2, quantile, p=0.95)
mean_mu <- colMeans(mu)

Risk_mu <- mean_mu / exp(offset)
Risk_0.05 <- q0.05 / exp(offset)
Risk_0.95 <- q0.95 / exp(offset)

Risk <- as.data.frame(cbind(Risk_mu,Risk_0.05,Risk_0.95))
# UnarmedBlack, ArmedBlack, UnarmedHispanic, ArmedHispanic, UnarmedWhite, ArmedWhite
Risk_BU <- Risk[1:51,] 
Risk_BA <- Risk[(1+51):(51+51),]
Risk_HU <- Risk[(1+51+51):(51+51+51),]
Risk_HA <- Risk[(1+51+51+51):(51+51+51+51),]

Risk_BU <- Risk_BU[order(Risk_BU$Risk_mu),]
Risk_BA <- Risk_BA[order(Risk_BA$Risk_mu),]
Risk_HU <- Risk_HU[order(Risk_HU$Risk_mu),]
Risk_HA <- Risk_HA[order(Risk_HA$Risk_mu),]



plot(Risk_BU$Risk_mu,main='Risk_BU',xlab='',ylab='',ylim=c(0,8),pch=19)
abline(h=mean(Risk_BU$Risk_mu),col='purple4',lty=3)
#abline(h=mean(Risk_BU$Risk_0.05),col='purple4',lty=3)
#abline(h=mean(Risk_BU$Risk_0.95),col='purple4',lty=3)
abline(h=quantile(Risk_BU$Risk_mu,0.05),col='purple4',lty=3)
abline(h=quantile(Risk_BU$Risk_mu,0.95),col='purple4',lty=3)

for(i in 1:51)
  segments(x0=i,y0=Risk_BU$Risk_0.05[i], y1=Risk_BU$Risk_0.95[i], col='grey4')
points(Risk_BU$Risk_mu,pch=19)








