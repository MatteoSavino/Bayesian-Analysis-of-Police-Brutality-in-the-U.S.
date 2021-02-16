#########################################
####         DATA ANALYSIS           #### 
#########################################

{
  library(tidyverse)
  # for STAN
  library(rstan)
  library(rethinking)
  library(coda)
  # for regression
  library(car)
  library(rgl)
  # for excel files
  library(readxl)
}

###############################

setwd('/Users/Francesco/Documents/POLICE_PROJECT')

data <- read_csv('fatal_police_WP.csv')    
names(data)

# We have data by cities
geo_data1 <- read_csv('fatal_police_kaggle/MedianHouseholdIncome2015.csv')
geo_data2 <- read_csv('fatal_police_kaggle/PercentagePeopleBelowPovertyLevel.csv')
geo_data3 <- read_csv('fatal_police_kaggle/PercentOver25CompletedHighSchool.csv')
geo_data4 <- read_csv('fatal_police_kaggle/ShareRaceByCity.csv')


# Now we compute data by states:
Income <- geo_data1 %>%
  mutate(`Median Income`= as.numeric(`Median Income`)) %>%
  filter(is.na(`Median Income`)==FALSE) %>%
   group_by( `Geographic Area`) %>%
  summarize(income=sum(`Median Income`)/length(`Median Income`))

Poverty <- geo_data2 %>%
  mutate(poverty_rate= as.numeric(poverty_rate)) %>%
  filter(is.na(poverty_rate)==FALSE) %>%
  group_by( `Geographic Area`) %>%
  summarize(poverty=sum(poverty_rate)/length(poverty_rate))

Education <- geo_data3 %>%
  mutate(percent_completed_hs= as.numeric(percent_completed_hs)) %>%
  filter(is.na(percent_completed_hs)==FALSE) %>%
  group_by( `Geographic Area`) %>%
  summarize(education=sum(percent_completed_hs)/length(percent_completed_hs))

Ethnicity <- geo_data4 %>%
  mutate(share_black = as.numeric(share_black),
         share_white = as.numeric(share_white)) %>%
  filter(is.na(share_black)==FALSE) %>%
  filter(is.na(share_white)==FALSE) %>%
  group_by(`Geographic area`) %>%
  summarize(black=sum(share_black)/length(share_black),
            white=sum(share_white)/length(share_white))

Datageo <- cbind(Income,poverty=Poverty$poverty,education=Education$education,
                 black=Ethnicity$black,white=Ethnicity$white)

colnames(Datageo)[1] <- 'state'
head(Datageo)

rm(geo_data1,geo_data2,geo_data3,geo_data4,Income,Poverty,Education,Ethnicity)

### Let us plot something
Datageo %>% 
  ggplot(aes(x=1:51,y=income)) +
  geom_col(col="cyan")

Datageo %>% 
  ggplot(aes(x=1:51,y=poverty)) +
  geom_col(col="magenta")

Datageo %>% 
  ggplot(aes(x=1:51,y=education)) +
  geom_col(col="yellow")

Datageo %>% 
  ggplot(aes(x=1:51,y=black)) +
  geom_col(col="grey")


### Let us investigate:
Datageo$State[which.max(Datageo$income)]
Datageo$State[which.min(Datageo$poverty)]
Datageo$State[which.max(Datageo$education)]
Datageo$State[which.max(Datageo$black)]

# ---
Datageo$State[which.min(Datageo$income)]
Datageo$State[which.max(Datageo$poverty)]
Datageo$State[which.min(Datageo$education)]
Datageo$State[which.min(Datageo$black)]



#### LINEAR MODELS #### --------------------------------------------------------------------

# black/white for every state
perc <- data %>%
  group_by(state) %>%
  summarise(percentage = length(which(race=='B'))/length(which(race=='W')))

### There is an outlier! It is DC:
plot(perc$percentage)
data %>% 
  filter(state=='DC') %>%
  summarise(black = length(which(race=='B')),
            white = length(which(race=='W')))
# We discard DC

Z <- Datageo[-8,]
Y <- perc[-8,]

plot(Y$percentage ~ Z$income)
plot(Y$percentage ~ Z$poverty)
plot(Y$percentage ~ Z$education)
plot(Y$percentage ~ Z$black)

# ------------------------------------------------
### We adopt a new approach: we are more selective:
perc <- data %>%
  group_by(state) %>%
  filter(length(which(race=='B'))>10) %>%
  filter(length(which(race=='W'))>10) %>%
  summarise(percentage = length(which(race=='B'))/length(which(race=='W')))

Z <- NULL
for(i in 1:28){
  for(j in 1:51){
    if(Datageo$state[j]==perc$state[i]){
      Z <- rbind(Z,Datageo[j,])
    } 
  }
}
Y <- perc
Z$state == Y$state

plot(Y$percentage ~ Z$income)
plot(Y$percentage ~ Z$poverty)
plot(Y$percentage ~ Z$education)
plot(Y$percentage ~ Z$black)


### LM ---------------------------
fit <- lm(Y$percentage ~ income + education + black, data=Z)
summary(fit)
fit <- lm(Y$percentage ~ income + poverty + education + black, data=Z)
summary(fit)
fit <- lm(Y$percentage ~ income + poverty + education + I(log(black)), data=Z)
summary(fit)
fit <- lm(Y$percentage ~ income +  education + I(log(black)), data=Z)
summary(fit)
fit <- lm(Y$percentage ~ income + I(log(black)), data=Z)
summary(fit)

vif(fit)
hist(fit$residuals)
shapiro.test(fit$residuals)
qqPlot(fit$residuals)


########################## -------------------------------------------------------
#### ELECTION DATASET ####
##########################

# xls files
election <- read_excel("DATA/usa_election.xlsx")
names(election)
#summary(election)
election <- election[,c(1,2,4,19,20,21,22,23,31,32,33,36,52,60:70)]
names(election)
colnames(election)[11] <- 'population'
election$population <- as.numeric(election$population)

#election2 <- election %>%
#  filter(is.na(Population)==0) %>%
#  filter(is.na(Black)==0) %>%
#  filter(Black <= 100) %>%
#  mutate(tot_black = Black / 100 * Population)

#plot(election2$Population) +  ### maybe some outliers!
#abline(h=1e+06, col='red') +
#points(which(election2$Population >= 1e+06),
#       election2$Population[which(election2$Population >= 1e+06)], pch=20, col='red')
#election2$County[which(election2$Population >= 1e+06)]

#plot(election2$tot_black)   ### maybe some outliers!
#which.max(election2$tot_black)
#which.max(election2$Population)
#rbind(election2$County[which.max(election2$tot_black)],
#      election2$Population[which.max(election2$tot_black)]
#      )

#election %>%
#  filter(is.na(Population)==0) %>%
# filter(is.na(Black)==0) %>%
#  mutate(tot_black = Black / 100 * Population) %>%
#  group_by(ST) %>%
#  summarise(Population=sum(Population),
#            tot_black=sum(tot_black),
#            perc = tot_black/Population*100) #%>%
  #summarise(tot_black=sum(tot_black)) %>%
  #summarise(perc = sum(tot_black) / sum(Population) * 100) #%>%
  #ggplot(aes(x=ST,y=perc_)) +
  #geom_col()


### We create some variables: ###
names(election)

Gini <- election %>%
  filter(is.na(Gini.Coefficient)==0) %>%
  mutate(Gini.Coefficient = as.numeric(Gini.Coefficient)) %>%
  mutate(gini_tot = Gini.Coefficient/100*population) %>%
  group_by(ST) %>%
  summarise(Gini = sum(gini_tot)/sum(population)*100) 
dim(Gini)
Gini <- Gini[-52,2]


# We create a unique dataset:
dim(Gini)
dim(Datageo)
X <- cbind(Datageo,Gini)
save(X,file = 'DATA.RDAta')
rm(X)
load('DATA.RData')



