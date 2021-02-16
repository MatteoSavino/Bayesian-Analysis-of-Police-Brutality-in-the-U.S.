library(dplyr) 
library(rstan) 
library(rethinking)

my_WAIC <- function(fit, param){
  llik   <- rstan::extract(fit, param)[[1]]
  p_WAIC <- sum(apply(llik, 2, var))
  lppd   <- sum(apply(llik, 2, function(x) log(mean(exp(x)))))
  WAIC   <- - 2 * lppd + 2 * p_WAIC
  return(WAIC)
}


# ------------------------------------------------------------------------#
############################### MODEL 1 ###################################
# ------------------------------------------------------------------------#
setwd('/Users/Francesco/Documents/POLICE_PROJECT')
data <- read.csv('DATA/fatal_police_WP.csv')      # fatal_police_kaggle/PoliceKillingsUS.csv   
                                     #                  OR
                                     #        fatal_police_WP.csv        
# Data are already ordered!
names(data)
head(data)
n <- dim(data)[1]

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

# or 
table(d$state,d$race) 

# We consider the states, not the counties!!!
table <- table(d$state,Results) 
dim(table)
table <- table[which(complete.cases(table[,'BlackArmed'])),] 
dim(table)

# Extract population size data for each state
# REM: Data from 2018:
usa <- read.csv('DATA/usa_population.csv',header=T) # open USA_populaton.csv 
dim(usa)

Nwhite <- usa$White
Nblack <- usa$Black
Nhispanic <- usa$Hispanic

# Extract shooting data for each state
UnarmedBlack <- as.numeric(table[,'BlackUnarmed'])
ArmedBlack <- as.numeric(table[,'BlackArmed'])
UnarmedHispanic <- as.numeric(table[,'HispanicUnarmed']	)
ArmedHispanic <- as.numeric(table[,'HispanicArmed'])
UnarmedWhite <- as.numeric(table[,'WhiteUnarmed']		)
ArmedWhite <- as.numeric(table[,'WhiteArmed'])

N<-length(Nblack)

model_dat  <- list(N=N,
                   mu_Mu=-15,
                   sigma2_Mu=10,
                   Nblack=Nblack,
                   Nwhite=Nwhite,
                   Nhispanic=Nhispanic,
                   UnarmedBlack=UnarmedBlack,	
                   ArmedBlack=ArmedBlack,		
                   UnarmedHispanic=UnarmedHispanic,		
                   ArmedHispanic=ArmedHispanic,		
                   UnarmedWhite=UnarmedWhite,		
                   ArmedWhite=ArmedWhite
)   


file_stan <- '/Users/Francesco/Documents/POLICE_PROJECT/STAN_MODELS/PAPER_MODEL.stan'
fitKilling <- stan(file=file_stan, data = model_dat, init=0,
                   thin=1, iter=1000, warmup=500,chains = 1,refresh=10)

my_WAIC(fitKilling, "log_lik_ab") 
my_WAIC(fitKilling, "log_lik_aw") 
my_WAIC(fitKilling, "log_lik_ah") 
my_WAIC(fitKilling, "log_lik_ub") 
my_WAIC(fitKilling, "log_lik_uw") 
my_WAIC(fitKilling, "log_lik_uh") 


### DIAGNOSTICS ###

#save(fitKilling, file="fitKilling_2chains.dat")
#load("STAN_MODELS/fitKilling_2chains.dat")

print(fitKilling, 
      probs = c(0.025, 0.5, 0.975), 
      par = c('Mu','RR_Black_Unarmed_Versus_White_Armed', 'RR_Black_Armed_Versus_Unarmed'))


library(coda)

coda_chain <- As.mcmc.list(fitKilling, pars=c(  "Mu_RR_Black_Armed_Versus_Unarmed",
                                                "Mu_RR_White_Armed_Versus_Unarmed",           
                                                "Mu_RR_Hispanic_Armed_Versus_Unarmed",        
                                                
                                                "Mu_RR_Black_Armed_Versus_White_Armed",       
                                                "Mu_RR_Hispanic_Armed_Versus_White_Armed",    
                                                "Mu_RR_Hispanic_Armed_Versus_Black_Armed",    
                                                
                                                "Mu_RR_Black_Unarmed_Versus_White_Unarmed",   
                                                "Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed",
                                                "Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed",
                                                
                                                "Mu_RR_Black_Unarmed_Versus_White_Armed",     
                                                "Mu_RR_Hispanic_Unarmed_Versus_White_Armed"))

coda_chain <- As.mcmc.list(fitKilling, pars='Theta')
                                                
summary(coda_chain)

# We need 2 chains!!!
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
gelman.plot(coda_chain)

# Geweke's convergence diagnostic
geweke.diag(coda_chain, frac1=0.1, frac2=0.5)
geweke.plot(coda_chain, frac1 = 0.1, frac2 = 0.5, nbins = 20) 

# autocorrelation and plot
quartz()
acfplot(coda_chain, lag.max = 300)

## Traceplot
rstan::traceplot(fitKilling, pars = "Mu_RR_Black_Armed_Versus_Unarmed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_White_Armed_Versus_Unarmed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Armed_Versus_Unarmed", inc_warmup = TRUE)

rstan::traceplot(fitKilling, pars = "Mu_RR_Black_Armed_Versus_White_Armed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Armed_Versus_White_Armed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Armed_Versus_Black_Armed", inc_warmup = TRUE)

rstan::traceplot(fitKilling, pars = "Mu_RR_Black_Unarmed_Versus_White_Unarmed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed", inc_warmup = TRUE)

rstan::traceplot(fitKilling, pars = "Mu_RR_Black_Unarmed_Versus_White_Armed", inc_warmup = TRUE)
rstan::traceplot(fitKilling, pars = "Mu_RR_Hispanic_Unarmed_Versus_White_Armed", inc_warmup = TRUE)

print(fitKilling, 
      probs = c(0.025, 0.5, 0.975), 
      par = "Mu_RR_Black_Unarmed_Versus_White_Armed")


### PLOT ###

#library(rgdal)
#usa_spdf <- readOGR('DATA/us_state/tl_2017_us_state.shp')
# Ran this code in order to plot good map
#library(sp)
#usa_spdf <- spTransform(gadm,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#library(broom)
#usa_fortified <- tidy(usa_spdf)

library(maps)
# We plot the map of the US!!!
usa <- map_data("state")
ggplot() +
  geom_map(data = usa, map = usa, aes(x = long, y = lat, map_id=region), fill='blue', col='black')
unique(usa$region)
# We retrieve the data we are interested in:
RR_BU_WU <- fitKilling %>% 
  rstan::extract("RR_Black_Unarmed_Versus_White_Unarmed")
RR_HU_WU <- fitKilling %>% 
  rstan::extract("RR_Hispanic_Unarmed_Versus_White_Unarmed")

RR_BU_WU <- colMeans(as.data.frame(RR_BU_WU))
RR_HU_WU <- colMeans(as.data.frame(RR_HU_WU))

# We merge the data!
load('DATA/DATA.RData')
X <- cbind(X,RR_BU_WU=RR_BU_WU,RR_HU_WU=RR_HU_WU)
X <- X %>%
  mutate(state_long=c('alaska','alabama','arkansas','arizona','california','colorado','connecticut',
                      'district of columbia','delaware','florida','georgia','hawaii',"iowa","idaho","illinois",
                      "indiana","kansas","kentucky","louisiana","massachusetts","maryland","maine",
                      "michigan","minnesota","missouri","mississippi",
                      "montana","nebraska","nevada", "new hampshire","new jersey","new mexico",
                      "new york","north carolina","north dakota","ohio","oklahoma","oregon",
                      "pennsylvania","rhode island","south carolina","south dakota","tennessee","texas",
                      "utah","vermont","virginia","washington","west virginia","wisconsin","wyoming"  ))

usa = usa %>%
  left_join(. , X, by=c("region"="state_long"))

# Finally we plot!!!

ggplot() +
  geom_polygon(data = usa, aes(fill = RR_HU_WU, x = long, y = lat, group = group)) +
  scale_fill_gradient(name = "H/W", low = "plum1", high = "magenta4", guide = "colorbar", na.value="#eeeeee") 
ggplot() +
  geom_polygon(data = usa, aes(fill = RR_BU_WU, x = long, y = lat, group = group)) +
  scale_fill_gradient(name = "H/W", low = "plum1", high = "magenta4", guide = "colorbar", na.value="#eeeeee") 







