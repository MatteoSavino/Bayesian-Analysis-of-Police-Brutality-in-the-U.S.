{
  #library(rstan) 
  #library(coda)
  #library(rethinking)
  # for plots:
  library(RColorBrewer)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(purrr)
  library(ggsci)
  require(gplots)
  require(ggpubr)
}

# Libraries for CAR models:
{
  library(glmnet)
  library(CARBayes)
  library(CARBayesdata)
  
  library(spdep)
  library(maptools) 
  library(sp)
  library(rgdal)
}

# ----------------------------------------------------------------------------#
############################### MODEL 3 CAR ###################################
# ----------------------------------------------------------------------------#

setwd('/Users/Francesco/Documents/POLICE_PROJECT')
data <- read.csv('DATA/fatal_police_WP2.csv')      # fatal_police_kaggle/PoliceKillingsUS.csv   

####################################
####        CALIFORNIA          ####
####################################

cali <- data %>%
  filter(state=='CA')
# We erase the outlier:
cali <- cali[-which.max(cali$longitude),]

#ca_spdf <- readOGR('DATA/CA_Counties/CA_Counties_TIGER2016.shp')
CA_shp <- readOGR('/Users/Francesco/Documents/POLICE_PROJECT/DATA/CALIFORNIA/tl_2016_06_place.shp')
names(CA_shp@data)
head(CA_shp@data)

# We can create a link between our data!!!
unique(cali$city)
unique(CA_shp@data$NAME)
length(which(cali$city %in% CA_shp@data$NAME))

`%notin%` <- Negate(`%in%`)
missing <- setdiff(CA_shp@data$ST, geodata$ST)
good_row <- which(usashp@data$ST %notin% missing)
USA_spatial <- usashp[good_row, ]
USA_spatial@data <- data.frame(USA_spatial@data, geodata, victims, Nwhite=usa$White, Nblack=usa$Black, Nhispanic=usa$Hispanic,
                               Black.rate=usa$Black/usa$White*100, Hispanic.rate=usa$Hispanic/usa$White*100)  
USA_spatial@data$ST == geodata$ST
is(USA_spatial)




# ______________________________________________________________________________
# ______________________________________________________________________________

  ####################
  #### CAR MODELS ####
  ####################

setwd('/Users/Francesco/Documents/POLICE_PROJECT')
data <- read.csv('DATA/fatal_police_WP2.csv')    
d <- data[which((data[,'armed'])!=''),]

Results<-c()
Results <- ifelse(d$armed !="unarmed" & d$race=="B", "BlackArmed",NA)
Results <- ifelse(d$armed =="unarmed" & d$race=="B", "BlackUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="H", "HispanicArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="H", "HispanicUnarmed",Results)
Results <- ifelse(d$armed !="unarmed" & d$race=="W", "WhiteArmed",Results)
Results <- ifelse(d$armed =="unarmed" & d$race=="W", "WhiteUnarmed",Results)

# We consider the states!!!
table <- table(d$state,Results) 
dim(table)
table <- table[which(complete.cases(table[,'BlackArmed'])),] 
dim(table)
head(table)
victims <- data.frame(BA=table[,1],
                      BU=table[,2],
                      HA=table[,3],
                      HU=table[,4],
                      WA=table[,5],
                      WU=table[,6])
head(victims)

usa <- read.csv('DATA/usa_population.csv',header=T) # open USA_populaton.csv 

# Now we download the shape files:

# '/Users/Francesco/Documents/POLICE_PROJECT/DATA/us_state/tl_2017_us_state.shp'
# '/Users/Francesco/Documents/POLICE_PROJECT/DATA/usa_shapefiles/cb_2018_us_nation_20m.shp'
usashp <- readOGR('/Users/Francesco/Documents/POLICE_PROJECT/DATA/usa_shapefile/cb_2018_us_state_500k.shp')
names(usashp)
unique(usashp$STUSPS)

# We order our shape data and change name to variable STUSPS:
usashp <- usashp[order(usashp$STUSPS),]
names(usashp)[5] <- 'ST'
head(usashp@data)

# Now we load our geo-data:

load('DATA/DATA_MAX.RData')
state_names <- unique(data$state)[order(unique(data$state))]
DATA_MAX$income <- scale(DATA_MAX$income)
DATA_MAX$`Median Age` <- scale(DATA_MAX$`Median Age`)
DATA_MAX$Violent.crime <- scale(DATA_MAX$Violent.crime)

geodata <- data.frame(ST=state_names,DATA_MAX)
rm(DATA_MAX)
head(geodata)

# We erase the superfluous observations and merge usashp with geodata:
# usashp + geodata = USA_spatial

`%notin%` <- Negate(`%in%`)
missing <- setdiff(usashp@data$ST, geodata$ST)
good_row <- which(usashp@data$ST %notin% missing)
USA_spatial <- usashp[good_row, ]
USA_spatial@data <- data.frame(USA_spatial@data, geodata, victims, Nwhite=usa$White, Nblack=usa$Black, Nhispanic=usa$Hispanic,
                               Black.rate=usa$Black/usa$White*100, Hispanic.rate=usa$Hispanic/usa$White*100)  
USA_spatial@data$ST == geodata$ST
is(USA_spatial)

head(USA_spatial@data)
# spplot(USA_spatial)

# We compute the similarity matrix:
W.nb <- poly2nb(USA_spatial, row.names = USA_spatial@data$ST)
# ATTENTION!!! Obviously Alaska and Hawaii have no connection!!!
# We better erase them:
rows <- which(USA_spatial$ST %in% c('AK','HI'))
USA_spatial <- USA_spatial[-rows, ]

W.nb <- poly2nb(USA_spatial, row.names = USA_spatial@data$ST)
W <- nb2mat(W.nb, style="B")

W[1:9,1:9]

### We study some scatter plots: ------------------------------------------------------

names(USA_spatial@data)

detach(USA_spatial@data)
attach(USA_spatial@data)

# We have a zero value for WA!!!
which(USA_spatial@data$WA==0)
USA_spatial@data$WA[which(USA_spatial@data$WA==0)] <- 1

which(USA_spatial@data$WU==0)
USA_spatial@data$WU[which(USA_spatial@data$WU==0)] <- 1

logRisk_BA <- log(USA_spatial@data$BA/USA_spatial@data$Nblack*USA_spatial@data$Nwhite/USA_spatial@data$WA)
logRisk_HA <- log(USA_spatial@data$HA/USA_spatial@data$Nhispanic*USA_spatial@data$Nwhite/USA_spatial@data$WA)

pairs(cbind(logRisk_BA, income, poverty, education, White, Black, Hispanic), 
      labels = c("R_BA", "income", "poverty", "education", "White", "Black", 'Hispanic'))
pairs(cbind(logRisk_BA, Gini.Coefficient, Median.Age, Unemployment, Violent.crime, Homicide.rate), 
      labels = c("R_BA", "Gini", "Median Age", "Unemployment", "Violence", "Homicide"))
pairs(cbind(logRisk_BA, Gini.Coefficient, Median.Age, Unemployment, Violent.crime, Homicide.rate), 
      labels = c("R_BA", "Gini", "Median Age", "Unemployment", "Violence", "Homicide"))
pairs(cbind(logRisk_BA, Nwhite, Nblack, Nhispanic, Black.rate, Hispanic.rate), 
      labels = c("R_BA", "Nwhite", "Nblack", "Nhispanic",'Black.rate','Hispanic.rate'))



#### HOW MANY ZEROS ??? ####

par(mfrow=c(3,2))
hist(victims$WA, col='darkgrey', breaks=50, main=length(which(victims$WA==0)))
hist(victims$WU, col='darkgrey', breaks=50, main=length(which(victims$WU==0)))
hist(victims$BA, col='purple4', breaks=50, main=length(which(victims$BA==0)))
hist(victims$BU, col='purple4', breaks=50, main=length(which(victims$BU==0)))
hist(victims$HA, col='gold', breaks=50, main=length(which(victims$HA==0)))
hist(victims$HU, col='gold', breaks=50, main=length(which(victims$HU==0)))

dev.off()

#### MODEL 0 #### ---------------------------------------------------------------


# BA: education, unemployment
# HA: "                     ", poverty
# HU: unemployment
# BU: poverty

names(USA_spatial)

formula = BA ~ income + poverty + education + Gini.Coefficient + Unemployment + # *
  Violent.crime  + offset(log(WA)-log(Nwhite)+log(Nblack))

formula = BA ~ income + poverty + education + Gini.Coefficient + # *
  Violent.crime  + offset(log(WA)-log(Nwhite)+log(Nblack))

formula = BU ~ income + poverty + education + Gini.Coefficient +  # /
  Violent.crime + offset(log(WU)-log(Nwhite)+log(Nblack))

formula = HA ~ income + poverty + education + Gini.Coefficient + # ***
  Violent.crime + offset(log(WA)-log(Nwhite)+log(Nhispanic))

formula = HU ~ income + poverty + education + Gini.Coefficient + # /
  Violent.crime + offset(log(WU)-log(Nwhite)+log(Nhispanic))


model_POI <- glm(formula = formula, family="poisson", data=USA_spatial@data)
summary(model_POI)
# We plot the Pearson residuals:
plot(residuals(model_POI, type='pearson'))


W.list <- nb2listw(W.nb, style = "B")

# H_0: no spatial autocorrelation, H_1: positive spatial autocorrelation
moran.mc(x = residuals(model_POI,type='pearson'), listw = W.list, nsim=10000)



#### MODEL 1 SPATIAL MODEL #### ---------------------------------------------------

names(USA_spatial)

#formula = BA ~ income + poverty + education + Gini.Coefficient + Violent.crime
  
# Given the high numbers of zero we used a ZIP distribution!
model_CAR_leroux <- S.CARleroux(formula = formula, 
                                          data=USA_spatial@data, 
                                          family="zip", W=W,
                                          formula.omega=~1,
                                          burnin=50000, n.sample=300000, 
                                          thin=10)
### Diagnostics

print(model_CAR_leroux)
model_CAR_leroux$model
names(model_CAR_leroux)
model_CAR_leroux$modelfit  # Goodness of fit indices
names(model_CAR_leroux$samples)

# Traceplots: 
plot(model_CAR_leroux$samples$beta)
plot(model_CAR_leroux$samples$rho) 
# rho ~ 1 ==> strong spatial correlation

summarise.samples(model_CAR_leroux$samples$beta, quantiles=c(0.5, 0.025, 0.975))

plot(model_CAR_leroux$samples$nu2)
plot(model_CAR_leroux$samples$tau2)
plot(model_CAR_leroux$samples$rho)

plot(model_CAR_leroux$samples$phi[,1:6], type='l') # spatial component phi


# Residuals (Pearson: i.e. standardized residuals for Poisson distribution)
plot(model_CAR_leroux$residuals[,2])

# We plot the fitted values:
plot(model_CAR_leroux$fitted.values, type = 'l', ylab='Black Unarmed Victims', xlab ='States') 
points(USA_spatial$BA, pch = 19, col = "purple4")


### We plot the USA

USA_spatial@data$risk <- model_CAR_leroux$fitted.values /
  ((USA_spatial@data$WA) * USA_spatial@data$Nblack / USA_spatial@data$Nwhite)

USA_spatial@data$risk <- model_CAR_leroux$fitted.values /
  ((USA_spatial@data$WA) * USA_spatial@data$Nhispanic / USA_spatial@data$Nwhite)


hist(USA_spatial@data$risk, col='purple4')

spplot(USA_spatial, c("risk"), main='Hispanic Armed Risk Rates',
       #sp.layout=list( boundary.final, col ='red', pch=20), #at=breakpoints, 
       par.settings=list(axis.line = list(col =  'transparent')),
       # col.regions=cm.colors(n=20),
       col="white")

graphics.off()



#### MODEL 2 DISSIMILARITY MODEL #### -----------------------------------------------

# MODEL fitted: Y_black ~ Poisson(E R), E_k = Y_white * Nblack / Nwhite
#               where Rk is the risk of being shot if you are black in state k and
#               ln(Rk) = beta_0 + \Phi_k, phi_k spatial random effect

# We compute some dissimilarity matrices:
income <- USA_spatial@data$income
poverty <- USA_spatial@data$poverty
education <- USA_spatial@data$education
Gini <- USA_spatial@data$Gini.Coefficient
crime <- USA_spatial@data$Violent.crime

Z.income <- as.matrix(dist(cbind(income, income), method="euclidean", diag=TRUE, upper = TRUE))
Z.education <- as.matrix(dist(cbind(education, education), method="euclidean", diag=TRUE, upper = TRUE))
Z.poverty <- as.matrix(dist(cbind(poverty, poverty), method="euclidean", diag=TRUE, upper = TRUE))
Z.Gini <- as.matrix(dist(cbind(Gini, Gini), method="euclidean", diag=TRUE, upper = TRUE))
Z.crime <- as.matrix(dist(cbind(crime, crime), method="euclidean", diag=TRUE, upper = TRUE))

names(USA_spatial)


formula = BA ~ offset(log(WA)-log(Nwhite)+log(Nblack))
formula = HA ~ offset(log(WA)-log(Nwhite)+log(Nhispanic))


model.dissimilarity <- S.CARdissimilarity(formula = formula, 
                       data=USA_spatial@data, 
                       family="poisson", W=W,
                       Z=list(Z.income=Z.income,
                              Z.poverty=Z.poverty,
                              Z.education=Z.education,
                              Z.Gini=Z.Gini,
                              Z.crime=Z.crime),
                       burnin=10000, n.sample=50000, thin=10)


print(model.dissimilarity)

# Diagnostics

plot(model.dissimilarity$samples$phi[,1:6], type='l')

# Residuals (Pearson: i.e. standardized residuals for Poisson distribution)
plot(model.dissimilarity$residuals[,2])

# We plot the fitted values:
plot(model.dissimilarity$fitted.values, type = 'l', ylab='Black Armed Victims', xlab ='States') 
points(USA_spatial$HA, pch = 19, col = "purple4")

# Estimated probability of having a border ----> exp(-z_kj \alpha) 
# (NA if they are not adjacent)
model.dissimilarity$localised.structure$W.border.prob[1:9, 1:9]
# If the probability is small, than a border is detected, 
# otherwise no
# i.e. if exp(-\alpha z_{ji}) > 0.5, there is NO border
model.dissimilarity$localised.structure$W.posterior[1:9, 1:9]

border.locations <- model.dissimilarity$localised.structure$W.posterior
boundary.final <- highlight.borders(border.locations=border.locations,
                                    spdata=USA_spatial)
USA_spatial@data$risk <- model.dissimilarity$fitted.values /
                                     ((USA_spatial@data$WA) * USA_spatial@data$Nblack / USA_spatial@data$Nwhite)

USA_spatial@data$risk <- model.dissimilarity$fitted.values /
                         ((USA_spatial@data$WA) * USA_spatial@data$Nhispanic / USA_spatial@data$Nwhite)

hist(USA_spatial@data$risk, col='gold')

spplot(USA_spatial, c("risk"), main='Hispanic Armed Risk Rates',
       sp.layout=list( boundary.final, col ='red', pch=20), #at=breakpoints, 
       par.settings=list(axis.line = list(col =  'transparent')),
       # col.regions=cm.colors(n=20),
       col="white")

# Pay attention to Illinois!


#### _________________________________________________________________________________

#### We focus just on the most interesting states:

# Using the DATA_PLOT script we know that the north east states are:

north_east <- c('Connecticut', 'Delaware','District Of Columbia',
                'Illinois','Indiana' ,             'Iowa' ,               
                'Kansas'   ,            'Kentucky'    ,         'Maine'  ,             
                'Maryland'  ,           'Massachusetts'   ,     'Michigan'  ,          
                'Minnesota' ,          'Missouri'   ,          'Nebraska'  ,          
                'New Hampshire' ,       'New Jersey'  ,         'New York'  ,          
                'North Dakota'  ,       'Ohio'      ,                    
                'Pennsylvania' ,        'Rhode Island',         'South Dakota' ,       
                'Vermont'       ,       'Virginia'       ,      'West Virginia' ,      
                'Wisconsin')
north_east_rows <- which(USA_spatial$NAME %in% north_east)
USA_spatial_NE@data <- USA_spatial@data[north_east_rows, ]

spplot(USA_spatial_NE, c("risk"),
       sp.layout=list( boundary.final, col ='red', pch=20), #at=breakpoints, 
       par.settings=list(axis.line = list(col =  'transparent')),
       # col.regions=cm.colors(n=20),
       col="white"
       )






